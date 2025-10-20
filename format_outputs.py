import argparse
import datetime
import os
from PIL import Image
import numpy as np
import tempfile
from osgeo import gdal
import shutil

from mosaic import apply_glt_noClick
from spec_io import write_cog, open_tif

product_metadata = {
    'pv': {
        'cov': {
            'name': 'EMIT_L2B_FRCOVPV',
            'description': 'Photosynthetic Vegetation Fractional Cover Values',
            'units': 'percent'
        },
        'unc': {
            'name': 'EMIT_L2B_FRCOVPVUNC',
            'description': 'Photosynthetic Vegetation Fractional Cover Uncertainty Values',
            'units': 'percent'
        }
    },
    'npv': {
        'cov': {
            'name': 'EMIT_L2B_FRCOVNPV',
            'description': 'Non-photosynthetic Vegetation Fractional Cover Values',
            'units': 'percent'
        },
        'unc': {
            'name': 'EMIT_L2B_FRCOVNPVUNC',
            'description': 'Non-photosynthetic Vegetation Fractional Cover Uncertainty Values',
            'units':  'percent'
        }
    },
    'bare': {
        'cov': {
            'name': 'EMIT_L2B_FRCOVBARE',
            'description': 'Bare Soil Fractional Cover Values',
            'units': 'percent'
        },
        'unc': {
            'name': 'EMIT_L2B_FRCOVBAREUNC',
            'description': 'Bare Soil Fractional Cover Uncertainty Values',
            'units':  'percent'
        }
    },
    'qc': {
            'name': 'EMIT_L2B_FRCOVQC',
            'description': 'Fractional Cover Quality Flag',
        }
}

def add_metadata_to_cog(input_file, product_metadata, software_build_version, product_version):
    metadata = {
        "keywords": "Imaging Spectroscopy, minerals, EMIT, dust, radiative forcing",
        "sensor": "EMIT (Earth Surface Mineral Dust Source Investigation)",
        "instrument": "EMIT",
        "platform": "ISS",
        "Conventions": "CF-1.63",
        "institution": "NASA Jet Propulsion Laboratory/California Institute of Technology",
        "license": "https://creativecommons.org/publicdomain/zero/1.0/",
        "naming_authority": "LPDAAC",
        "date_created": datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ"),
        "keywords_vocabulary": "NASA Global Change Master Directory (GCMD) Science Keywords",
        "stdname_vocabulary": "NetCDF Climate and Forecast (CF) Metadata Convention",
        "creator_name": "Jet Propulsion Laboratory/California Institute of Technology",
        "creator_url": "https://earth.jpl.nasa.gov/emit/",
        "project": "Earth Surface Mineral Dust Source Investigation",
        "project_url": "https://earth.jpl.nasa.gov/emit/",
        "publisher_name": "NASA LPDAAC",
        "publisher_url": "https://lpdaac.usgs.gov",
        "publisher_email": "lpdaac@usgs.gov",
        "identifier_product_doi_authority": "https://doi.org",
        "software_build_version": software_build_version,
        "product_version": product_version,
        "description": product_metadata.get("description", "")
    }

    ds_mem = gdal.Translate("", input_file, format="MEM")
    ds_mem.SetMetadata(metadata)

    band = ds_mem.GetRasterBand(1)
    band.SetDescription(product_metadata.get("name", ""))
    if "units" in product_metadata:
        band.SetMetadataItem("UNITS", product_metadata["units"])
    band.FlushCache()

    translate_options = gdal.TranslateOptions(format="COG")
    with tempfile.NamedTemporaryFile(suffix=".tif", delete=False) as tmp:
        gdal.Translate(tmp.name, ds_mem, options=translate_options)
        shutil.move(tmp.name, input_file)
        os.chmod(input_file, 0o664)

    ds_mem = None

def apply_mask(frcov_file, frcov_unc_file, mask_file, output_base, glt_file, software_version, product_version, glt_nodata_value):
    
    """
    Apply GLT and mask to fractional cover and uncertainty files.

    Args:
        frcov_file (str): path to input fractional cover file
        frcov_unc_file (str): path to input fractional cover uncertainty file
        mask_file (str): path to mask file
        output_base (str): base path for saving output files
        glt_file (str): path to EMIT GLT file
        glt_nodata_value (int): nodata value for GLT file (default=0)
    """

    output_directory = os.path.dirname(output_base)
    os.makedirs(output_directory, exist_ok=True)

    add_metadata_to_cog(mask_file, product_metadata['qc'], software_version, product_version)  

    ortho_frcov_file = output_base + '_frcov_ort.tif'
    apply_glt_noClick(glt_file, frcov_file, ortho_frcov_file, nodata_value=-9999,
                      bands=None, output_format='tif', glt_nodata_value=glt_nodata_value)

    ortho_frcov_unc_file = output_base + '_frcov_unc_ort.tif'
    apply_glt_noClick(glt_file, frcov_unc_file, ortho_frcov_unc_file, nodata_value=-9999,
                      bands=None, output_format='tif', glt_nodata_value=glt_nodata_value)

    _, mask = open_tif(mask_file)

    frcov_meta, frcov = open_tif(ortho_frcov_file)
    frcov[mask[:,:,0] > 0] = -9999

    frcov_unc_meta, frcov_unc = open_tif(ortho_frcov_unc_file)
    frcov_unc[mask[:,:,0] > 0] = -9999

    cover_types = ['npv', 'pv', 'bare']
    
    for band in range(3):
        masked_ortho_frcov_file = os.path.join(output_directory, output_base + f'_frcov_{cover_types[band]}.tif')
        masked_ortho_frcov_unc_file = os.path.join(output_directory, output_base  + f'_frcovunc_{cover_types[band]}.tif')
        
        write_cog(masked_ortho_frcov_file, frcov[:,:,[band]], frcov_meta)
        write_cog(masked_ortho_frcov_unc_file, frcov_unc[:,:,[band]], frcov_unc_meta)
        
        add_metadata_to_cog(masked_ortho_frcov_file, product_metadata[cover_types[band]]['cov'], software_version, product_version)  
        add_metadata_to_cog(masked_ortho_frcov_unc_file, product_metadata[cover_types[band]]['unc'], software_version, product_version)  

    bare = frcov[:,:,2]
    pv = frcov[:,:,1]
    npv = frcov[:,:,0]
    rgb = np.dstack([bare, pv, npv])
    rgb[rgb == -9999] = 0
    rgb = np.clip(rgb, 0, 1)
    rgb = (rgb * 255).astype(np.uint8)
    alpha = np.where(mask[:,:,0] == -9999, 0, 255).astype(np.uint8)
    rgba = np.dstack([rgb, alpha])
    png_path = os.path.join(output_directory, output_base + '_frcov.png')
    Image.fromarray(rgba, mode='RGBA').save(png_path)
    
    bare_color = np.array([210, 180, 140], dtype=float)
    pv_color = np.array([34, 139, 34], dtype=float)
    npv_color = np.array([218, 165, 32], dtype=float)

    rgb = (bare[:,:,None] * bare_color +
           pv[:,:,None] * pv_color +
           npv[:,:,None] * npv_color)

    nodata_mask = (mask[:,:,0] == -9999)
    rgb[nodata_mask] = 0
    rgb = np.clip(rgb, 0, 255).astype(np.uint8)
    alpha = np.where(nodata_mask, 0, 255).astype(np.uint8)
    rgba = np.dstack([rgb, alpha])
    
    png_path = os.path.join(output_directory, output_base + '_frcov_col.png')
    Image.fromarray(rgba, mode='RGBA').save(png_path)
 
 
def main():
    parser = argparse.ArgumentParser(description='Apply GLT and mask to fractional cover data')
    parser.add_argument('frcov_file', type=str)
    parser.add_argument('frcov_unc_file', type=str)
    parser.add_argument('mask_file', type=str)
    parser.add_argument('glt_file', type=str)
    parser.add_argument('output_base', type=str)
    parser.add_argument('--software_version', type=str)
    parser.add_argument('--product_version', type=str)
    parser.add_argument('--glt_nodata_value', type=int, default=0)
    args = parser.parse_args()

    apply_mask(args.frcov_file, 
               args.frcov_unc_file, 
               args.mask_file, 
               args.output_base, 
               args.glt_file,
               args.software_version,
               args.product_version, 
               args.glt_nodata_value)


if __name__ == '__main__':
    main()