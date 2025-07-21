import json 
import click
import glob
import os 
from osgeo import gdal, ogr 
import numpy as np

import spectral_util 
import mosaic

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

    # Orthorectify EMIT mask file 
    ortho_mask_file = os.path.join(output_loc, 'masks', fid + '_ortho_mask.tif')
    glt_nodata_value = 0
    mosaic.apply_glt.main([glt_file, mask_file, ortho_mask_file, '--glt_nodata_value', glt_nodata_value], standalone_mode=False)

    # Urban 
    urban_out_file = os.path.join(output_loc, 'masks', fid + '_urban.tif')
    urban_mask(ortho_mask_file, json_filename, urban_data_loc, urban_out_file)

    # Coastal 
    coastal_out_file = os.path.join(output_loc, 'masks', fid + '_coastal.tif')
    coastal_mask(json_filename, coastal_data_loc, coastal_out_file)
    
    # NDSI + RGB (generate + ortho)
    ndsi_file = os.path.join(output_loc, 'masks', fid + '_ndsi.tif')
    rgb_file = os.path.join(output_loc, 'masks', fid + '_rgb.tif')
    spectral_util.ndsi.main([rfl_file, ndsi_file], standalone_mode=False)
    spectral_util.rgb.main([rfl_file, rgb_file], standalone_mode=False)

    ndsi_ortho_file = os.path.join(output_loc, 'masks', fid + '_ortho_ndsi.tif')
    rgb_ortho_file = os.path.join(output_loc, 'masks', fid + '_ortho_rgb.tif')
    mosaic.apply_glt.main([glt_file, ndsi_file, ndsi_ortho_file, '--glt_nodata_value', glt_nodata_value], standalone_mode=False)
    mosaic.apply_glt.main([glt_file, rgb_file, rgb_ortho_file, '--glt_nodata_value', glt_nodata_value], standalone_mode=False)


#### #### 

def urban_mask(ortho_file, json_file, urban_data, output_file, output_res = 0.000542232520256, nodata_value = 0):

    # Get SRS info from orthoed file 
    ds = gdal.Open(ortho_file)
    wkt = ds.GetProjection()
    ds = None

    # Build warp options
    warp_options = gdal.WarpOptions(
        cutlineDSName=json_file,
        cropToCutline=True,
        dstNodata=nodata_value,
        xRes=output_res,
        yRes=-output_res,
        dstSRS=wkt
    )
    gdal.Warp(destNameOrDestDS="clipped.tif", srcDSOrSrcDSTab=urban_data, options=warp_options)

    # Generate geotiff mask of urban areas (50 in ESA worldcover)
    ds = gdal.Open("clipped.tif") # temporary file 
    band = ds.GetRasterBand(1)
    array = band.ReadAsArray()
    result = np.logical_and(array >= 0, array == 50).astype(np.uint8)

    # Create output file with the same georeference and projection
    driver = gdal.GetDriverByName("GTiff")
    out_ds = driver.Create(output_file, ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Byte)
    out_ds.SetGeoTransform(ds.GetGeoTransform())
    out_ds.SetProjection(ds.GetProjection())

    # Write result and set NoData value
    out_band = out_ds.GetRasterBand(1)
    out_band.WriteArray(result)
    out_band.SetNoDataValue(nodata_value)

    os.remove("clipped.tif") # del temporary file 

    print('...Writing urban mask')

def coastal_mask(json_file, coastal_data, output_file, output_res = 0.000542232520256): 
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

    # Create output raster
    out_ds = gdal.GetDriverByName("GTiff").Create(output_file, x_res, y_res, 1, gdal.GDT_Byte)
    out_ds.SetGeoTransform((minx, output_res, 0, maxy, 0, -output_res))
    out_band = out_ds.GetRasterBand(1)
    out_band.Fill(1)

    # Set projection and rasterize 
    srs = layer.GetSpatialRef()
    if srs:
        out_ds.SetProjection(srs.ExportToWkt())

    gdal.RasterizeLayer(
        out_ds,
        [1],  # band index
        layer,
        burn_values=[0],
        options=["INVERT=FALSE"]
    )

    os.remove("temp.shp") # del temporary file 

    print('...Writing coastal mask')

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