import argparse
import os

from mosaic import apply_glt_noClick
from spec_io import write_cog, open_tif


def apply_mask(frcov_file, frcov_unc_file, mask_file, output_base, glt_file, glt_nodata_value):

    output_directory = os.path.dirname(output_base)
    os.makedirs(output_directory, exist_ok=True)


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
        masked_ortho_frcov_file = os.path.join(output_directory, output_base + f'_frcovunc_{cover_types[band]}.tif')
        masked_ortho_frcov_unc_file = os.path.join(output_directory, output_base  + f'_frcov_{cover_types[band]}.tif')
        write_cog(masked_ortho_frcov_file, frcov[:,:,[band]], frcov_meta)
        write_cog(masked_ortho_frcov_unc_file, frcov_unc[:,:,[band]], frcov_unc_meta)
        
    os.remove(ortho_frcov_file)
    os.remove(ortho_frcov_unc_file)

def main():
    parser = argparse.ArgumentParser(description='Apply GLT and mask to fractional cover data')
    parser.add_argument('frcov_file', type=str)
    parser.add_argument('frcov_unc_file', type=str)
    parser.add_argument('mask_file', type=str)
    parser.add_argument('glt_file', type=str)
    parser.add_argument('output_base', type=str)
    parser.add_argument('--glt_nodata_value', type=int, default=0)
    args = parser.parse_args()

    apply_mask(args.frcov_file, args.frcov_unc_file, args.mask_file, args.output_base, args.glt_file, args.glt_nodata_value)


if __name__ == '__main__':
    main()