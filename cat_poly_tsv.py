#this script aggregates together run.tsv and sample.tsv files from Polyphemus. Pandas is used as columns are not
#always in the same order

import argparse
import os
import pandas as pd

from pathlib import Path


def aggregate_tsv(suffix, out_dir):

    main_path = Path('/mnt/prod_bfx_analysis_ro/analysis_data/progenity_workflow_polyphemus/flowcell_result')
    #main_path = Path('/mnt/prod_ro/sequencing/ops/production/nipt9002/results')  # old prod location, pre 2.0

    run_filename = out_dir / f'poly_run_{suffix}.tsv'
    sample_filename = out_dir / f'poly_sample_{suffix}.tsv'

    fc = os.listdir(main_path)

    all_run = []
    all_sample = []

    print('Reading in files...')

    for f in fc:

        run_path = main_path / f / 'tsv' / 'run.tsv'
        sample_path = main_path / f / 'tsv' / 'samples.tsv'
        #run_path = main_path / f / 'Bioinformatics_Project_Poly_9002' /'tsv' / 'Run_Project_Poly_9002_run.tsv' # old prod location, pre 2.0

        if os.path.isfile(run_path):
            run_tsv = pd.read_csv(run_path, sep='\t', header = 0)
            all_run.append(run_tsv)

        if os.path.isfile(sample_path):
            sample_tsv = pd.read_csv(sample_path, sep='\t', header=0)
            all_sample.append(sample_tsv)

    print('Writing out files...')
    poly_run = pd.concat(all_run, axis=0, sort=False, join='inner')
    poly_run.to_csv(run_filename, sep='\t', index=None)

    poly_run = pd.concat(all_sample, axis=0, sort=False, join='inner')
    poly_run.to_csv(sample_filename, sep='\t', index=None)

    print('Complete!')
    print(run_filename)
    print(sample_filename)


def arg_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument('suffix', type=str, help='output file suffix (ex. date like 1009)')
    parser.add_argument('output_path', type=Path, help='output path')

    return parser


def cli():
    parser = arg_parser()
    args = parser.parse_args()

    aggregate_tsv(args.suffix, args.output_path)


if __name__ == "__main__":
    cli()
