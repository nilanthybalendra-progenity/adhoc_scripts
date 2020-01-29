import pandas as pd
import os
from pathlib import Path
import argparse


def get_file_list(main_path, fc):
    """flowcell result, sample alignment and sample result directories contain directories for each flowcell/sample
    appended by "v##". The latest output does not necessarily have the greatest ## value. This function puts together a
    list of the latest and second latest run directories based on creation time stamp.
    Args:
        main_path (Path): Path to directory to search
        fc (string): Flowcell ID
    Returns: List List
    """
    files = pd.DataFrame()
    files['BFX_ID'] = [s for s in os.listdir(main_path) if s.split('_')[0] in fc]
    files['SAMPLE_ID'] = files['BFX_ID'].str[:-4] # remove version tag

    files['TIME'] = [os.stat(main_path / s).st_ctime for s in files['BFX_ID'].tolist()]

    files.sort_values(by='TIME', inplace=True)
    latest = files.drop_duplicates(subset='SAMPLE_ID', keep='last')

    left_over = files.loc[~files['BFX_ID'].isin(latest['BFX_ID'].tolist())]

    second_latest = left_over.drop_duplicates(subset='SAMPLE_ID', keep='last')

    return str(main_path) +'/'+ latest['BFX_ID'], str(main_path) +'/'+second_latest['BFX_ID']


def aggregate(fc, file_suff, output_path):
    """Generates the production output file list and candidate output file list given flowcell ID. Then the files
    indicated by file_suff are aggregated for each release and output as written out, as well as the two file lists.

    Args:
        fc (string): Flowceel ID
        file_suff(string): File suffix (general, rt or cnv)
        output_path(Path): Output path

    Returns:

    """
    main_path = Path('/mnt/bfx_analysis_ro/analysis_data/progenity_workflow_levitate/sample_result')

    file = f'{file_suff}.tsv'

    sorter = {'general.tsv': 'LIS_RESULT_NAME',
              'rt.tsv':      ['BFX_RT_BFXID', 'BFX_RT_ALLELEID'],
              'cnv.tsv':     ['BFX_NEBULA_BFXID', 'BFX_NEBULA_ALLELEID']}


    cand_path, prod_path = get_file_list(main_path, fc)

    prod_path.to_csv(output_path / 'production_dirs.txt', index=None, header=False)
    cand_path.to_csv(output_path / 'candidate_dirs.txt', index=None, header=False)

    p_tmp = []
    c_tmp = []
    print('aggregating files')
    for p_file, c_file in zip(prod_path.to_list(), cand_path.to_list()):
        p_tmp.append(pd.read_csv(Path(p_file) / 'results' / file, sep='\t', header=0))
        c_tmp.append(pd.read_csv(Path(c_file) / 'results' / file, sep='\t', header=0))

    p_all = pd.concat(p_tmp, axis=0, join='outer', sort=True).sort_values(by=sorter[file])
    c_all = pd.concat(c_tmp, axis=0, join='outer', sort=True).sort_values(by=sorter[file])

    if file == 'general.tsv':
        print('dropping known different columns')
        cols = ['BFX_RESULT_CREATEDDATETIME', 'BFX_RESULT_SAMPLEALIGNMENTRESULTURI', 'LIS_PIPELINE_VERSION', 'LIS_REQUEST_CREATEDDATETIME', 'LIS_RESULT_REQUESTGUID', 'LIS_RESULT_RESULTGUID']
        p_all = p_all.drop(labels=cols, axis=1)
        c_all = c_all.drop(labels=cols, axis=1)

    if file == 'cnv.tsv':
        p_all['BFX_NEBULA_GOF'] = p_all['BFX_NEBULA_GOF'].round(7)
        c_all['BFX_NEBULA_GOF'] = c_all['BFX_NEBULA_GOF'].round(7)

    print('outputting files')
    p_all.to_csv(output_path / f'{fc}_production_{file}', sep='\t',index=None)
    c_all.to_csv(output_path / f'{fc}_candidate_{file}', sep='\t',index=None)


def arg_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument('fcid', type=str, help='flowcell id')
    parser.add_argument('file', type=str, help='file type to aggregate, "general", "rt" or "cnv"')
    parser.add_argument('output_path', type=Path, help='Path to output directory')

    return parser


def cli():
    parser = arg_parser()
    args = parser.parse_args()

    aggregate(args.fcid, args.file, args.output_path)


if __name__ == '__main__':
    cli()