import json 
import click
import glob
import os 
from osgeo import gdal, ogr 
import numpy as np

from mosaic import apply_glt
from spec_io import load_data, write_cog, open_tif

### 
@click.command()
@click.argument('fid', type=str, required=True)
@click.argument('input_loc', type=click.Path(exists=True), default="/store/emit/ops/data/acquisitions/")
@click.argument('output_loc', type=click.Path(exists=True), default="/home/colemanr/Unmixing/outputs/")
@click.argument('urban_data_loc', type=click.Path(exists=True), default="/store/shared/landcover/complete_landcover.vrt")
@click.argument('coastal_data_loc', type=click.Path(exists=True), default="/home/colemanr/Unmixing/coastal_mask/GSHHS_f_L1.shp")
@click.argument('json_file', type=click.Path(exists=True), default="/store/brodrick/emit/emit-visuals/track_coverage_pub.json")
def process_files(fid, input_loc, output_loc, urban_data_loc, coastal_data_loc, json_file):

    # Write specific fid to geojson string
    with open(json_file) as json_data:
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

    ############ Generate masks ############

    # # Orthorectify EMIT mask file 
    ortho_mask_file = os.path.join(output_loc, 'masks', fid + '_ortho_mask.tif')
    glt_nodata_value = 0
    apply_glt.main([glt_file, mask_file, ortho_mask_file, '--glt_nodata_value', glt_nodata_value], standalone_mode=False)

    # Urban 
    urban_out_file = os.path.join(output_loc, 'masks', fid + '_ortho_urban.tif')
    meta = urban_mask_cog(ortho_mask_file, urban_out_file, json_filename, urban_data_loc)

    # Coastal 
    coastal_out_file = os.path.join(output_loc, 'masks', fid + '_ortho_coastal.tif')
    coastal_mask_cog(json_filename, coastal_out_file, coastal_data_loc, meta)
    
    # NDSI (generate + ortho)
    ndsi_file = os.path.join(output_loc, 'masks', fid + '_ndsi.tif')
    ndsi(rfl_file, ndsi_file)

    ndsi_ortho_file = os.path.join(output_loc, 'masks', fid + '_ortho_ndsi.tif')
    apply_glt.main([glt_file, ndsi_file, ndsi_ortho_file, '--glt_nodata_value', glt_nodata_value], standalone_mode=False)

    ## Stack raster masks and remove individual files?

    os.remove(ndsi_file)


#### 

def urban_mask_cog(ortho_file, out_file, json_file, urban_data, output_res = 0.000542232520256, nodata_value = 0):

    print(f"Running Urban Masking on {json_file}")

    # Get SRS info from orthoed file 
    ds = gdal.Open(ortho_file)
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

    write_cog(out_file, result, meta)
    os.remove(temp_file)
    return meta 


def coastal_mask_cog(json_file, out_file, coastal_data, meta, output_res = 0.000542232520256): 

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

    # Get extent from shapefile for new raster 
    shp_ds = ogr.Open("temp.shp")
    layer = shp_ds.GetLayer()   
    minx, maxx, miny, maxy = layer.GetExtent()
    x_res = int((maxx - minx) / output_res)
    y_res = int((maxy - miny) / output_res)

    # Create in-memory raster
    mem_ds = gdal.GetDriverByName("MEM").Create("", x_res, y_res, 1, gdal.GDT_Byte)
    geotransform = (minx, output_res, 0, maxy, 0, -output_res)
    mem_ds.SetGeoTransform(geotransform)

    # Set projection from shapefile
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

    # Write COG
    write_cog(out_file, result, meta)

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
    Calculate NDSI (normalized difference snow index) -- not orthorectified yet 

    Args:
        input_file (str): Path to the input file.
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
    Convert a GeoJSON geometry string and save as a json in correct format 
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


##########



@click.group()
def cli():
    pass

cli.add_command(process_files)

if __name__ == '__main__':
    cli()