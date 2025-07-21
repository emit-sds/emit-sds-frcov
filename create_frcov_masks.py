import json 
import click
import glob
import os 
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
    mosaic.urban_mask.main([ortho_mask_file, json_filename, urban_data_loc, urban_out_file], standalone_mode=False)

    # Coastal 
    coastal_out_file = os.path.join(output_loc, 'masks', fid + '_coastal.tif')
    mosaic.coastal_mask.main([json_filename, coastal_data_loc, coastal_out_file], standalone_mode=False)
    
    # NDSI + RGB (generate + ortho)
    ndsi_file = os.path.join(output_loc, 'masks', fid + '_ndsi.tif')
    rgb_file = os.path.join(output_loc, 'masks', fid + '_rgb.tif')
    spectral_util.ndsi.main([rfl_file, ndsi_file], standalone_mode=False)
    spectral_util.rgb.main([rfl_file, rgb_file], standalone_mode=False)

    ndsi_ortho_file = os.path.join(output_loc, 'masks', fid + '_ortho_ndsi.tif')
    rgb_ortho_file = os.path.join(output_loc, 'masks', fid + '_ortho_rgb.tif')
    mosaic.apply_glt.main([glt_file, ndsi_file, ndsi_ortho_file, '--glt_nodata_value', glt_nodata_value], standalone_mode=False)
    mosaic.apply_glt.main([glt_file, rgb_file, rgb_ortho_file, '--glt_nodata_value', glt_nodata_value], standalone_mode=False)


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