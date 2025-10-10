import json
import click
import os
from osgeo import gdal
import numpy as np

import subprocess
import geopandas as gpd
from affine import Affine
import rasterio
from rasterio.features import rasterize

from mosaic import apply_glt_noClick
from spec_io import load_data, write_cog, open_tif


###
@click.command()
@click.argument('rfl_file', type=str, required=True)
@click.argument('l2a_mask_file', type=str, required=True)
@click.argument('glt_file', type=str, required=True)
@click.argument('frcov_mask', type=str, required=True)
@click.option('--urban_data', type=click.Path(exists=True), default="/store/shared/landcover/complete_landcover.vrt")
@click.option('--coastal_data', type=click.Path(exists=True), default="/store/shared/landcover/GSHHS_f_L1.shp")
@click.option('--glt_nodata_value', type=int, default = 0)
def process_files(rfl_file, l2a_mask_file, glt_file, frcov_mask, urban_data, coastal_data, glt_nodata_value):
    """
    Generate QC product for EMIT fractional cover

    Writes single band COG with following values:
        Cloud (EMIT cloud + cirrus flag)    = 1
        Urban                               = 2
        Water (EMIT water + coastal mask)   = 3
        Snow/Ice                            = 4

    Args:
        acq_id (str):
        input_loc (path str): path to EMIT reflectance, mask, GLT files
        output_loc (path str): path to save generated mask output files
        urban_data (path str): path to ESA worldcover dataset (.vrt/tif)
        coastal_data (path str): path to GSHHS coastal dataset (.shp)
        glt_nodata_value (int): defaults to 0 (nodata for .envi files)
    """

    output_directory = os.path.dirname(frcov_mask)
    glt_nodata_value = 0

    acq_id = os.path.basename(l2a_mask_file).split('_')[0]

    os.makedirs(output_directory, exist_ok=True)

    ############ Generate QC and save to COG ############

    # Orthorectify EMIT mask file
    ortho_mask_file = os.path.join(output_directory, acq_id + 'l2amask_ortho.tif')
    apply_glt_noClick(glt_file, l2a_mask_file, ortho_mask_file, nodata_value=-9999,
                      bands=None, output_format='tif', glt_nodata_value=glt_nodata_value)

    # Write ortho'ed extent to json file
    json_filename = os.path.join(output_directory, acq_id + '_extent.json')
    geotiff_extent_to_geojson(ortho_mask_file, json_filename)

    # Urban mask and orth
    urban_out_file = os.path.join(output_directory, acq_id + '_ortho_urban.tif')
    meta = urban_mask_cog(ortho_mask_file, urban_out_file, json_filename, urban_data, ortho_mask_file)

    # Coastal mask and ortho
    coastal_out_file = os.path.join(output_directory, acq_id + '_ortho_coastal.tif')
    coastal_mask_cog(ortho_mask_file, json_filename, coastal_out_file, coastal_data, meta)

    # NDSI (generate and then ortho)
    ndsi_file = os.path.join(output_directory, acq_id + '_ndsi.tif')
    ndsi_cog(rfl_file, ndsi_file)

    ndsi_ortho_file = os.path.join(output_directory, acq_id + '_ortho_ndsi.tif')
    apply_glt_noClick(glt_file, ndsi_file, ndsi_ortho_file, nodata_value=-9999, bands=None, output_format='tif', glt_nodata_value=glt_nodata_value)

    ## Convert to singleband COG
    _, urban_mask = open_tif(urban_out_file)
    _, coastal_mask = open_tif(coastal_out_file)
    _, ndsi_mask = open_tif(ndsi_ortho_file)

    emit_meta, emit_mask = open_tif(ortho_mask_file)
    emit_cloud = emit_mask[:,:,9] # SpecTf cloud flag
    emit_cirrus = emit_mask[:,:,1]
    emit_water = emit_mask[:,:,2]

    ## Convert to singleband
    single_band_stack = os.path.join(frcov_mask)
    singleband_raster_hierarchy(emit_cloud, emit_cirrus, emit_water,
                                urban_mask[:,:,0], ndsi_mask[:,:,0], coastal_mask[:,:,0],
                                single_band_stack, emit_meta)

    ## Clean up and remove intermediary files
    os.remove(ndsi_file)
    os.remove(ndsi_ortho_file)
    os.remove(ortho_mask_file)
    os.remove(coastal_out_file)
    os.remove(urban_out_file)
    os.remove(json_filename)

####

def geotiff_extent_to_geojson(tiff_path, geojson_path):
    """
    Extracts the bounding box of a GeoTIFF file and saves it as a GeoJSON file.

    Args:
        tiff_path (str): Path to the input GeoTIFF file.
        geojson_path (str): Path to the output GeoJSON file.
    """

    # Open raster and get extent of valid data
    with rasterio.open(tiff_path) as src:
        bounds = src.bounds
        crs = src.crs

    # Extract bounding box
    polygon = {
        "type": "Polygon",
        "coordinates": [[
            [bounds.left, bounds.bottom],
            [bounds.left, bounds.top],
            [bounds.right, bounds.top],
            [bounds.right, bounds.bottom],
            [bounds.left, bounds.bottom]  # close the polygon
        ]]
    }

    # Write GeoJSON to file
    feature = {
        "type": "Feature",
        "geometry": polygon,
        "properties": {
            "crs": str(crs)
        }
    }
    geojson = {
        "type": "FeatureCollection",
        "features": [feature]
    }
    with open(geojson_path, 'w') as f:
        json.dump(geojson, f, indent=2)


def singleband_raster_hierarchy(cloud, cirrus, water, urban, snow_ice, coastal, out_file, meta):
    """
    Condense multiple row x col arrays into a single band COG with a hierarchical classification process

    Writes single band COG with following values:
        Cloud (EMIT cloud + cirrus flag)    = 1
        Urban                               = 2
        Water (EMIT water + coastal mask)   = 3
        Snow/Ice                            = 4

    # --- hierarchy order --- #
    if cloud or cirrus or (cloud + cirrus):
        QC = 1
    if urban:
        QC = 2
    if water or coastal or (water + coastal):
        QC = 3
    if snow/ice:
        QC = 4

    Args:
        cloud, cirrus, water, urban, snow_ice, coastal (arr): 6 row x col arrays, where 1 = value to be masked out for that variable
        out_file (str): path to save singleband raster output
        meta (GenericGeoMetadata): An object containing the wavelengths and FWHM.
    """

    # apply hierarchical categorization logic
    result = np.zeros((cloud.shape[0], cloud.shape[1]), dtype=np.uint8)
    result[(cloud == 1) | (cirrus == 1)] = 1
    result[(urban == 1) & (result == 0)] = 2
    result[((water == 1) | (coastal == 1)) & (result == 0)] = 3
    result[(snow_ice == 1) & (result == 0)] = 4

    result = result.reshape((result.shape[0], result.shape[1], 1))
    write_cog(out_file, result, meta)

def warp_array_to_ref(array, source_ds, ref_path, nodata_value=0):
    """
    Warp input array to match projection/resolution of reference gtiff

    Args:
        array (np arr): input array to reproject
        source_ds (gdal dataset): input dataset of array
        ref_path (str): path to reference tif
        nodata_value (int): nodata value to for gdal dataset

    Out:
        out_arr (np arr): gdalwarped input array
    """

    ref_ds = gdal.Open(ref_path)
    if ref_ds is None:
        raise FileNotFoundError(f"Could not open {ref_path}")

    # Save in memory
    mem_ds = gdal.GetDriverByName("MEM").Create('', source_ds.RasterXSize, source_ds.RasterYSize, 1, gdal.GDT_UInt16)
    mem_ds.SetGeoTransform(source_ds.GetGeoTransform())
    mem_ds.SetProjection(source_ds.GetProjection())
    mem_ds.GetRasterBand(1).WriteArray(array[:, :, 0])
    mem_ds.GetRasterBand(1).SetNoDataValue(nodata_value)

    # Extract spatial transform
    gt = ref_ds.GetGeoTransform()
    bounds = (gt[0], gt[3] + gt[5]*ref_ds.RasterYSize, gt[0] + gt[1]*ref_ds.RasterXSize, gt[3])

    # apply gdal transform to desired reference dataset
    warp_options = gdal.WarpOptions(format='MEM',
                            dstSRS=ref_ds.GetProjection(),
                            outputBounds=bounds,
                            width=ref_ds.RasterXSize,
                            height=ref_ds.RasterYSize,
                            srcNodata=nodata_value,
                            dstNodata=nodata_value,
                            resampleAlg='bilinear')
    warped = gdal.Warp('', mem_ds, options=warp_options)
    out_arr = warped.ReadAsArray().reshape(ref_ds.RasterYSize, ref_ds.RasterXSize,1)

    return out_arr


def urban_mask_cog(ortho_file, out_file, json_file, urban_data, ref_path, output_res = 0.000542232520256, nodata_value = 0):
    """
    Generate mask of urban/built-up areas and save as COG

    Args:
        ortho_file (str): path to orthorectified EMIT mask file
        out_file (str): path to save urban area COG
        json_file (str): path to json of EMIT tile extent
        urban_data (str): path to ESA worldcover dataset (.vrt/tif)
        ref_path (str): path to reference tif to align data with
        output_res (float): default to EMIT res
        nodata_value (int): nodata value for gdal dataset

    Out:
        meta (GenericGeoMetadata): An object containing the wavelengths and FWHM.
    """

    print(f"Running Urban Masking on {json_file}")

    # Get SRS info from orthoed file
    ds_mask = gdal.Open(ortho_file)
    if ds_mask is None:
        raise FileNotFoundError(f"Could not open {ortho_file}")
    wkt = ds_mask.GetProjection()

    # Build warp options -- coarse clipping to bounding box
    temp_file= os.path.join(os.path.dirname(out_file), os.path.splitext(os.path.basename(out_file))[0]) + '_TEMPclipped.tif'
    warp_options = gdal.WarpOptions(
        cutlineDSName=json_file,
        cropToCutline=True,
        dstNodata=nodata_value,
        xRes=output_res,
        yRes=-output_res,
        dstSRS=wkt
    )
    gdal.Warp(destNameOrDestDS=temp_file, srcDSOrSrcDSTab=urban_data, options=warp_options)

    # Generate geotiff mask of urban areas (50 in ESA worldcover)
    meta, _ = open_tif(temp_file)
    ds = gdal.Open(temp_file)
    band = ds.GetRasterBand(1)
    urban_array = band.ReadAsArray()
    result = np.logical_and(urban_array >= 0, urban_array == 50).astype(np.uint8)
    result = result.reshape((result.shape[0], result.shape[1], 1))

    # Exact clipping to valid data points in EMIT data mask
    emit_mask = (ds_mask.GetRasterBand(1).ReadAsArray() != -9999)
    emit_mask = emit_mask.reshape((emit_mask.shape[0], emit_mask.shape[1], 1))
    result_clip = np.where(emit_mask, result, 0)
    result_warp = warp_array_to_ref(result_clip, ds, ref_path)

    # Write to COG
    write_cog(out_file, result_warp, meta)

    os.remove(temp_file)
    return meta


def coastal_mask_cog(ortho_file, json_file, out_file, coastal_data, meta, output_res = 0.000542232520256):
    """
    Generate mask of coastal water features and save as COG

    Args:
        ortho_file (str): path to orthorectified EMIT mask file
        json_file (str): path to json of EMIT tile extent
        out_file (str): path to save coastal area COG
        coastal_data (str): path to GSHHS coastal dataset (.shp)
        meta (GenericGeoMetadata):  An object containing the wavelengths and FWHM.
        output_res (float): default to EMIT res
    """

    print(f"Running Coastal Masking on {json_file}")

    # Clip large coastal data to approx. tile extent
    tile_extent = gpd.read_file(json_file)
    coastal = gpd.read_file(coastal_data)
    clipped = gpd.overlay(coastal, tile_extent, how="intersection")

    # Get extent from json
    minx, miny, maxx, maxy = tile_extent.total_bounds
    width, height = int((maxx - minx) / output_res), int((maxy - miny) / output_res)
    transform = Affine.translation(minx, maxy) * Affine.scale(output_res, -output_res)

    if clipped.empty:  # No intersecting coastal features -- return mask of 0
        raster = np.zeros((height, width), dtype=np.uint8)

    else:
        #  Mask for inside EMIT tile = 1, outside tile = 0 --> needed to prevent classification as water in corners
        ds_mask = gdal.Open(ortho_file)
        if ds_mask is None:
            raise FileNotFoundError(f"Could not open {ortho_file}")
        tile_mask = (ds_mask.GetRasterBand(1).ReadAsArray() != -9999)

        # Land = 0,  Water = 1
        coastal_raster = rasterize(
            [(geom, 0) for geom in clipped.geometry if not geom.is_empty],
            out_shape=(height, width),
            transform=transform,
            fill=1,
            dtype=np.uint8
        )
        raster = coastal_raster * tile_mask

    # Write coastal mask to COG
    result = raster.reshape((height, width, 1))
    write_cog(out_file, result, meta)


def ndsi_cog(input_file, output_file, green_wl = 560, swir_wl = 1600, green_width = 0, swir_width = 0, threshold = 0.4, ortho=True):
    """
    Calculate NDSI (normalized difference snow index) and save as cog

    Args:
        input_file (str): Path to the EMIT reflectance
        output_file (str): Path to save NDSI COG
        green_wl (int): Green band wavelength [nm].
        swir_wl (int): SWIR1 band wavelength [nm].
        green_width (int): Green band width [nm]; 0 = single wavelength.
        swir_width (int): SWIR1 band width [nm]; 0 = single wavelength.
    """

    print(f"Running NDSI Calculation on {input_file}")

    meta, rfl = load_data(input_file, lazy=True, load_glt=ortho)

    green = rfl[..., meta.wl_index(green_wl, green_width)]
    swir = rfl[..., meta.wl_index(swir_wl, swir_width)]

    ndsi = (green - swir) / (green + swir)
    ndsi = ndsi.squeeze()
    ndsi[np.isfinite(ndsi) == False] = -9999
    ndsi = ndsi.reshape((ndsi.shape[0], ndsi.shape[1], 1))

    ndsi[ndsi > threshold] = 1
    ndsi[ndsi <= threshold] = 0

    write_cog(output_file, ndsi, meta, ortho=ortho)


## NOT CURRENTLY USED
def singleband_raster_unique(raster_stack, out_file):
    """
    Leverage distinct sums (aka 2^n Sidon set) to condense a multiband raster into a single band without losing
    information about pixels that are QC flagged for multiple reasons (e.g., both a cloud and an urban pixel).
    There are 63 distinct combinations for the current 6 band QC product.

    If a given image has N bands, then the values attributed to each band in the single-band raster are {2^0, 2^1, 2^2, ... 2^N}

    Band 1 - Cloud = 1
    Band 2 - Cirrus = 2
    Band 3 - Water = 4
    Band 4 - Urban = 8
    Band 5 - Snow/Ice = 16
    Band 6 - Coastal = 32

    Example --> pixel value of 52 = pixel flagged as water, snow/ice, coastal (4+16+32 = 52)
    """

    subprocess.run([
        'gdal_calc.py',
        '-A', raster_stack, '--A_band=1',
        '-B', raster_stack, '--B_band=2',
        '-C', raster_stack, '--C_band=3',
        '-D', raster_stack, '--D_band=4',
        '-E', raster_stack, '--E_band=5',
        '-F', raster_stack, '--F_band=6',
        '--calc', '(A*1)+(B*2)+(C*4)+(D*8)+(E*16)+(F*32)',
        '--outfile', out_file,
        '--NoDataValue=0',
        '--type=Int8',
        '--co', 'COMPRESS=LZW', ## remainder are from write_cog code
        '--co', 'BIGTIFF=YES',
        '--co', 'COPY_SRC_OVERVIEWS=YES',
        '--co', 'TILED=YES',
        '--co', 'BLOCKXSIZE=256',
        '--co', 'BLOCKYSIZE=256',
        '--overwrite'
    ], check=True)


##########


@click.group()
def cli():
    pass

cli.add_command(process_files)

if __name__ == '__main__':
    cli()