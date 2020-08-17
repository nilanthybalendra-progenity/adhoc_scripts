#this script aggregates together run.tsv and sample.tsv files from Polyphemus. Pandas is used as columns are not
#always in the same order
from pathlib import Path

import pandas as pd
import os

main_path = Path('/mnt/prod_bfx_analysis_ro/analysis_data/progenity_workflow_polyphemus/flowcell_result')
#main_path = Path('/mnt/prod_ro/sequencing/ops/production/nipt9002/results')

run_filename = 'poly_run_0814.tsv'
sample_filename = 'poly_sample_0814.tsv'

fc = os.listdir(main_path)

all_run = []
all_sample = []
for f in fc:
    print(f)

    run_path = main_path / f / 'tsv' / 'run.tsv'
    sample_path = main_path / f / 'tsv' / 'samples.tsv'
    #run_path = main_path / f / 'Bioinformatics_Project_Poly_9002' /'tsv' / 'Run_Project_Poly_9002_run.tsv'

    if os.path.isfile(run_path):
        run_tsv = pd.read_csv(run_path, sep='\t', header = 0)
        all_run.append(run_tsv)

    if os.path.isfile(sample_path):
        sample_tsv = pd.read_csv(sample_path, sep='\t', header=0)
        all_sample.append(sample_tsv)


poly_run = pd.concat(all_run, axis=0, sort=False, join='inner')
poly_run.to_csv(run_filename, sep='\t', index=None)

poly_run = pd.concat(all_sample, axis=0, sort=False, join='inner')
poly_run.to_csv(sample_filename, sep='\t', index=None)