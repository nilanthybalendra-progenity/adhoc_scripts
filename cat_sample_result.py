# this script cats together the sample result files for comparison in the V&V report
import os
import pandas as pd

from pathlib import Path

def get_file_list(main_path, fc):
    """flowcell result, sample alignment and sample result directories contain directories for each flowcell/sample
    appended by "v##". The latest output does not necessarily have the greatest ## value. This function puts together a
    list of the latest run directories based on creation time stamp.

    Args:
        main_path (Path): Path to directory to search
        fc (string): Flowcell ID

    Returns: List

    """
    files = pd.DataFrame()
    files['BFX_ID'] = [s for s in os.listdir(main_path) if s.split('_')[0] in fc]
    files['SAMPLE_ID'] = files['BFX_ID'].str[:-4] # remove version tag

    files['TIME'] = [os.stat(main_path / s).st_ctime for s in files['BFX_ID'].tolist()]

    files.sort_values(by='TIME', inplace=True)
    files.drop_duplicates(subset='SAMPLE_ID', keep='last', inplace=True)

    return files['BFX_ID'].tolist()


main_path = Path('/mnt/bfx_analysis_ro/analysis_data/progenity_workflow_levitate/sample_result')
fc = ['HCNJLDRXX', 'HCWM2DRXX', 'HCNLGDRXX']



general = []
cnv = []
rt = []

output_dir = Path('/mnt/ruo_rw/rnd/SQA_Data/output/progenity_workflow_levitate/vv_report_v1.0.6/objective_4')

for f in fc:
    dirs = get_file_list(main_path, f)

    for d in dirs:
        print(d)
        tmp_general = pd.read_csv(main_path / d / 'results' / 'general.tsv', sep='\t', header=0)
        tmp_rt = pd.read_csv(main_path / d / 'results' / 'rt.tsv', sep='\t', header=0)
        tmp_cnv = pd.read_csv(main_path / d / 'results' / 'cnv.tsv', sep='\t', header=0)

        general.append(tmp_general)
        rt.append(tmp_rt)
        cnv.append(tmp_cnv)


    all_gen = pd.concat(general, axis=0, sort=True)
    all_rt = pd.concat(rt, axis=0, sort=True)
    all_cnv = pd.concat(cnv, axis=0, sort=True)

    all_gen.to_csv(output_dir / f'{f}_v1.0.6_candidate_general.tsv', sep='\t', index=None)
    all_rt.to_csv(output_dir / f'{f}_v1.0.6_candidate_rt.tsv', sep='\t', index=None)
    all_cnv.to_csv(output_dir / f'{f}_v1.0.6_candidate_cnv.tsv', sep='\t', index=None)
    print(f'done_{f}')




