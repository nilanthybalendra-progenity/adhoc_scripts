'''This script facilitates TSV diffs between the RND Prototype and Bifrost'''
# test commands
# python diff_tool.py bridge-lot /mnt/ruo_rw/rnd/staff/nilanthy.balendra/bifrost/diff_test/b_expected/ /mnt/ruo_rw/rnd/staff/nilanthy.balendra/bifrost/diff_test/b_results/
# python diff_tool.py calibrator-assign /mnt/ruo_rw/rnd/staff/nilanthy.balendra/bifrost/diff_test/c_expected/ /mnt/ruo_rw/rnd/staff/nilanthy.balendra/bifrost/diff_test/c_results/

import argparse
import os
import toml
import pandas as pd

from pathlib import Path
from pandas.testing import assert_frame_equal


def special(filename, command, bifrost_data, prototype_data):
    """Perform file specific format manipulation"""

    if (filename == 'results.tsv') & (command == 'calibrator-assign'):
        bifrost_data['conc_filename']   = bifrost_data['conc_filename'].str[:-4] #remove .tsv
        prototype_data['conc_filename'] = prototype_data['conc_filename'].str[:-10] # remove _trans.tsv
    
    elif filename == 'selected_samples.tsv':
        bifrost_data['conc_trans_filename']   = bifrost_data['conc_trans_filename'].str[:-4]
        prototype_data['conc_trans_filename'] = prototype_data['conc_trans_filename'].str[:-10]
        bifrost_data['Dataset'] = bifrost_data['Dataset'].str.capitalize()

    elif (filename in ['inter_controls_cv.tsv', 'intra_controls_cv.tsv']):
        bifrost_data['Dataset'] = bifrost_data['Dataset'].str.capitalize()
    
    elif (filename == 'results.tsv') & (command == 'bridge-lot'):
        bifrost_data['conc filename']   = bifrost_data['conc filename'].str[:-4]
        prototype_data['conc filename'] = prototype_data['conc filename'].str[:-10]
        bifrost_data['Dataset'] = bifrost_data['Dataset'].str.capitalize() # capitalize Comparator and Incoming

    return bifrost_data, prototype_data


def check_exists(command, prototype_dir, bifrost_dir, prototype_files, bifrost_files):
    """Checks if files exists in the prototype and bifrost directories"""

    check_file = pd.DataFrame()

    check_file['file_paths'] = [(prototype_dir / f) for f in prototype_files] + [(bifrost_dir / f) for f in bifrost_files] 
    check_file['files'] = prototype_files + bifrost_files
    check_file['exists'] = [os.path.isfile(f) for f in check_file['file_paths']]

    missing = check_file.loc[check_file['exists'] == False]

    if not missing.empty:
        raise ValueError(f'Files are missing: {missing["files"].tolist()}')


def check_df_equal(params, command, prototype_dir, bifrost_dir):
    """Perform operations to address known differences, then check if equal"""
    prototype = pd.read_csv(prototype_dir/ params['prototype_file'])
    bifrost = pd.read_csv(bifrost_dir / params['bifrost_file'], sep='\t')

    # make prototype column names match bifrost
    col_map = {params['prototype_col'][i]: params['bifrost_col'][i] for i in range(len(params['prototype_col']))} 
    prototype = prototype.rename(columns=col_map)

    # make percents in prototype into ratios 
    to_ratio = params['to_ratio']
    prototype[to_ratio] = prototype[to_ratio] / 100 

    # drop extra bifrost columns
    if params['drop']:
        bifrost = bifrost.drop(params['drop_bcols'], axis=1)
        prototype = prototype.drop(params['drop_pcols'], axis=1)

    # make column order match bifrost
    col_order = bifrost.columns.values.tolist()
    prototype = prototype[col_order]

    # perform file specific operations
    if params['special']:
        bifrost, prototype = special(params['bifrost_file'], command, bifrost, prototype)
    
    # sort rows to match
    if params['sort']:
        bifrost = bifrost.sort_values(by=params["by"], ignore_index=True)
        prototype  = prototype.sort_values(by=params["by"], ignore_index=True)

    assert_frame_equal(bifrost, prototype, check_like=True, atol=1e-3)


def arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()

    parser.add_argument('command', type=str, help='Outputs being compared (bridge-lot or calibrator-assign')
    parser.add_argument('prototype_dir', type=Path, help='Path to prototype result directory')
    parser.add_argument('bifrost_dir', type=Path, help='Path to bifrost result directory')

    return parser


def run():
    args = arg_parser().parse_args()

    config = toml.load('diff_config.toml')
    
    if args.command == 'bridge-lot': 
        params = config['bridge-lot']
        check_exists(args.command, args.prototype_dir, args.bifrost_dir, params['prototype_files'], params['bifrost_files'])

        check_df_equal(params['inter_controls_cv'], args.command, args.prototype_dir, args.bifrost_dir)
        check_df_equal(params['intra_controls_cv'], args.command, args.prototype_dir, args.bifrost_dir)
        check_df_equal(params['od_means'], args.command, args.prototype_dir, args.bifrost_dir)
        check_df_equal(params['conc_difference_stats'], args.command, args.prototype_dir, args.bifrost_dir)
        check_df_equal(params['conc_difference'], args.command, args.prototype_dir, args.bifrost_dir)
        check_df_equal(params['selected_samples'], args.command, args.prototype_dir, args.bifrost_dir)
        check_df_equal(params['results'], args.command, args.prototype_dir, args.bifrost_dir)

    elif args.command == 'calibrator-assign':
        params = config['calibrator-assign']
        check_exists(args.command, args.prototype_dir, args.bifrost_dir, params['prototype_files'], params['bifrost_files'])

        check_df_equal(params['concentration_per_exp'], args.command, args.prototype_dir, args.bifrost_dir)
        check_df_equal(params['intravials'], args.command, args.prototype_dir, args.bifrost_dir)
        check_df_equal(params['vials'], args.command, args.prototype_dir, args.bifrost_dir)
        check_df_equal(params['mean'], args.command, args.prototype_dir, args.bifrost_dir)
        check_df_equal(params['results'], args.command, args.prototype_dir, args.bifrost_dir)
        
    else:  
        raise ValueError('Invalid command, must be either bridge-lot or calibrator-assign')

if __name__ == "__main__":
    run()



