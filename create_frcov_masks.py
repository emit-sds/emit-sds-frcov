import json 
import click
import glob
import os 
from osgeo import gdal, ogr 
import numpy as np

import subprocess

from mosaic import apply_glt
from spec_io import load_data, write_cog, open_tif

# TODO - clip to EMIT tile extent or convert to NaN outside tile extent 

### 
@click.command()
@click.argument('fid', type=str, required=True)
@click.argument('input_loc', type=click.Path(exists=True), default="/store/emit/ops/data/acquisitions/")
@click.argument('output_loc', type=click.Path(exists=True), default="/home/colemanr/Unmixing/outputs/")
@click.argument('urban_data_loc', type=click.Path(exists=True), default="/store/shared/landcover/complete_landcover.vrt")
@click.argument('coastal_data_loc', type=click.Path(exists=True), default="/home/colemanr/Unmixing/coastal_mask/GSHHS_f_L1.shp")
@click.argument('json_file_loc', type=click.Path(exists=True), default="/store/brodrick/emit/emit-visuals/track_coverage_pub.json")
@click.argument('glt_nodata_value', type=int, default = 0)
def process_files(fid, input_loc, output_loc, urban_data_loc, coastal_data_loc, json_file_loc, glt_nodata_value):
    """
    Generate QC product for EMIT fractional cover 

    Writes single band COG with following values: 
        Cloud (EMIT cloud + cirrus flag)    = 1
        Urban                               = 2
        Water (EMIT water + coastal mask)   = 3
        Snow/Ice                            = 4

    Args: 
        fid (str):
        input_loc (path str): path to EMIT reflectance, mask, GLT files
        output_loc (path str): path to save generated mask output files 
        urban_data_loc (path str): path to ESA worldcover dataset (.vrt/tif)
        coastal_data_loc (path str): path to GSHHS coastal dataset (.shp)
        json_file_loc (path str): path to EMIT info json file 
        glt_nodata_value (int): defaults to 0 (nodata for .envi files) 
    """

    # Write specific fid to geojson string
    with open(json_file_loc) as json_data:
        coverage = json.load(json_data)
    for feature in coverage["features"]:
        if feature["properties"].get("fid") == fid:
            geojson_str = json.dumps(feature["geometry"])
    json_filename = os.path.join(output_loc, 'json_tile', fid + '_extent.json')
    geojson_str_to_feature_file(geojson_str, fid, json_filename)

    # Get file names
    try: 
        rfl_file = glob.glob(os.path.join(input_loc, fid[4:12], fid.split('_')[0], "l2a", "*l2a_rfl_*.hdr"))[0]
        mask_file = glob.glob(os.path.join(input_loc, fid[4:12], fid.split('_')[0], "l2a", "*l2a_mask_*.hdr"))[0]
        glt_file = glob.glob(os.path.join(input_loc, fid[4:12], fid.split('_')[0], "l1b", "*l1b_glt_*.hdr"))[0]
    except Exception as e: 
        print('No path found for FID')

    ############ Generate QC and save to COG ############

    # Orthorectify EMIT mask file 
    ortho_mask_file = os.path.join(output_loc, 'masks', fid + '_ortho_mask.tif')
    apply_glt.main([glt_file, mask_file, ortho_mask_file, '--glt_nodata_value', glt_nodata_value], standalone_mode=False)

    # Urban mask and orth
    urban_out_file = os.path.join(output_loc, 'masks', fid + '_ortho_urban.tif')
    ref_path = ortho_mask_file
    meta = urban_mask_cog(ortho_mask_file, urban_out_file, json_filename, urban_data_loc, ref_path)

    # Coastal mask and ortho
    coastal_out_file = os.path.join(output_loc, 'masks', fid + '_ortho_coastal.tif')
    coastal_mask_cog(json_filename, coastal_out_file, coastal_data_loc, meta, ref_path)
    
    # NDSI (generate and then ortho)
    ndsi_file = os.path.join(output_loc, 'masks', fid + '_ndsi.tif')
    ndsi(rfl_file, ndsi_file)

    ndsi_ortho_file = os.path.join(output_loc, 'masks', fid + '_ortho_ndsi.tif')
    apply_glt.main([glt_file, ndsi_file, ndsi_ortho_file, '--glt_nodata_value', glt_nodata_value], standalone_mode=False)

    ## Convert to singleband COG 
    _, urban_mask = open_tif(urban_out_file)
    _, coastal_mask = open_tif(coastal_out_file)
    _, ndsi_mask = open_tif(ndsi_ortho_file)

    emit_meta, emit_mask = open_tif(ortho_mask_file)
    emit_cloud = emit_mask[:,:,0]
    emit_cirrus = emit_mask[:,:,1]
    emit_water = emit_mask[:,:,2]

    ## Convert to singleband 
    single_band_stack = os.path.join(output_loc, 'masks', fid + '_ortho_hierarchy.tif')
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
    write_cog(out_file, result, meta, nodata_value=0) # nodata value of 0 b/c of np array initialization

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
    mem_ds = gdal.GetDriverByName("MEM").Create('', source_ds.RasterXSize, source_ds.RasterYSize, 1, gdal.GDT_Byte)
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
    ds = gdal.Open(ortho_file)
    if ds is None:
        raise FileNotFoundError(f"Could not open {ortho_file}")
    wkt = ds.GetProjection()
    ds = None
    temp_file = "clipped.tif"

    # Build warp options
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

    # Align to ref_path (EMIT mask file) and write to COG
    result_warp = warp_array_to_ref(result, ds, ref_path)
    write_cog(out_file, result_warp, meta)

    os.remove(temp_file)
    return meta 


def coastal_mask_cog(json_file, out_file, coastal_data, meta, ref_path, output_res = 0.000542232520256): 
    """
    Generate mask of coastal water features and save as COG
    
    Args:
        json_file (str): path to json of EMIT tile extent
        out_file (str): path to save coastal area COG
        coastal_data (str): path to GSHHS coastal dataset (.shp)
        meta (GenericGeoMetadata):  An object containing the wavelengths and FWHM.
        ref_path (str): path to reference tif to align data with
        output_res (float): default to EMIT res 
    """

    print(f"Running Coastal Masking on {json_file}")

    # Clip large coastal shapefile to EMIT extent (too slow) 
    gdal.VectorTranslate(
        "temp.shp",                              # Output file
        coastal_data,                            # Input file (coastal data)
        options=gdal.VectorTranslateOptions(
            format="ESRI Shapefile",
            clipSrc=json_file                    # Clip to tile extent 
        )
    )

    # Get extent from json 
    json_ds = ogr.Open(json_file)
    json_layer = json_ds.GetLayer()
    minx, maxx, miny, maxy = json_layer.GetExtent()
    x_res = int((maxx - minx) / output_res)
    y_res = int((maxy - miny) / output_res)
    
    # Create in-memory raster
    mem_ds = gdal.GetDriverByName("MEM").Create("", x_res, y_res, 1, gdal.GDT_Byte)
    geotransform = (minx, output_res, 0, maxy, 0, -output_res)
    mem_ds.SetGeoTransform(geotransform)

    # Set projection from shapefile
    shp_ds = ogr.Open("temp.shp")
    layer = shp_ds.GetLayer()   
    srs = layer.GetSpatialRef()
    if srs:
        mem_ds.SetProjection(srs.ExportToWkt())
    else:
        raise ValueError("Shapefile has no spatial reference.")

    # Initialize with land (1), then burn coastal (0)
    mem_ds.GetRasterBand(1).Fill(1)
    gdal.RasterizeLayer(mem_ds, [1], layer, burn_values=[0], options=["INVERT=FALSE"])

    # Read result as NumPy array
    result = mem_ds.GetRasterBand(1).ReadAsArray()
    result = result.reshape((result.shape[0], result.shape[1], 1))

    # Align to ref_path (EMIT mask file) and write to COG
    result_warp = warp_array_to_ref(result, mem_ds, ref_path) 
    write_cog(out_file, result_warp, meta)

    # Clean up temp files
    shp_ds = None
    mem_ds = None
    for ext in [".shp", ".shx", ".dbf", ".prj"]:
        try:
            os.remove(f"temp{ext}")
        except FileNotFoundError:
            pass


def ndsi(input_file, output_file, green_wl = 560, swir_wl = 1600, green_width = 0, swir_width = 0, threshold = 0.4, ortho=True):
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

def geojson_str_to_feature_file(geojson_str, fid, output_filename):
    """
    Convert a GeoJSON geometry string to save as a json

    Args:
        geojson_str (str): geojson string describing EMIT tile extent
        fid (str): EMIT fid
        output_filename (str): output json file 
    """
    geometry = json.loads(geojson_str)
    
    feature = {
        "type": "Feature",
        "geometry": geometry,
        "properties": {
            "FID": fid
        }
    }

    feature_collection = {
        "type": "FeatureCollection",
        "features": [feature]
    }

    with open(output_filename, 'w') as f:
        json.dump(feature_collection, f, indent=2)

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