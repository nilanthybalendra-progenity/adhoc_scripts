# this is a script to aggregate all nipt counts for a given flowcell
# python agg_counts.py /mnt/prod_ro/sequencing/ops/production/nipt9002/results/ HHNF7DMXX /mnt/ruo_rw/rnd/staff/nilanthy.balendra/tools/adhoc_scripts/
# neb_train_wrapper wraps this function to allow it to submit jobs to the cluster

import argparse
import os
import pandas as pd
from pathlib import Path


def aggregate_counts(args):
    fc = args.fc
    print(fc)
    out_dir = args.out_dir

    main_dir = Path('/mnt/prod_ro/sequencing/ops/production/nipt9002/results/')
    samp_dir = main_dir / f'{fc}.v1' / 'Project_Poly_9002'

    if os.path.isdir(samp_dir): # then the data is in the old location:
        samples = os.listdir(samp_dir) #all samples
        place = 'old'

    else: #it is in the new spot
        main_dir = Path('/mnt/prod_bfx_analysis_ro/analysis_data/progenity_workflow_polyphemus/sample_result')
        files = pd.DataFrame()
        files['BFX_ID'] = [s for s in os.listdir(main_dir) if s.split('_')[0] in fc]
        files['SAMPLE_ID'] = files['BFX_ID'].str[:-4]  # remove version tag

        files['TIME'] = [os.stat(main_dir / s).st_ctime for s in files['BFX_ID'].tolist()]
        files.sort_values(by='TIME', inplace=True)
        files.drop_duplicates(subset='SAMPLE_ID', keep='last', inplace=True) #take the latest run of the samples

        samples = files['BFX_ID'].tolist()
        place = 'new'

    appended_count_files = []
    ind = 0
    for sid in samples:
        if place == 'old':
            counts_file = samp_dir / sid / f'{sid[7:]}_counts.tsv'
        else:
            counts_file = main_dir / sid / 'counts' / 'counts.tsv'

        if os.path.isfile(counts_file):
            counts = pd.read_csv(counts_file, sep='\t', index_col=0, header=0)
            if ind > 0:
                counts.drop(labels='ROPType', axis=1, inplace=True)
            ind +=1
            appended_count_files.append(counts)

        # if ind ==5:
        #     break

    all_counts = pd.concat(appended_count_files, axis=1)

    #all_counts.to_csv(out_dir / f'{fc}.tsv',sep='\t')
    all_counts.to_parquet(out_dir / f'{fc}.pq', index=True)


def arg_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument('fc', type=str, help='flowcell ID')
    parser.add_argument('out_dir', type=Path, help='path to output directory')

    return parser


def cli():
    parser = arg_parser()
    args = parser.parse_args()

    aggregate_counts(args)


if __name__ == '__main__':
    cli()