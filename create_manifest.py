import numpy as np
import os
import pandas as pd

from pathlib import Path


def control_list(data):
    tmp = data.loc[data['CONTROL_SAMPLE'] != 'Test', ['PROPS_ID', 'KNOWN_PLOIDY', 'HOST_SEX', 'FETAL_SEX', 'WELL', 'RUN', 'CONTROL_SAMPLE']]

    tmp.loc[tmp['WELL'] == 'H12', 'CONTROL_SAMPLE'] = 'NTC'
    tmp.loc[tmp['WELL'] == 'H12', 'KNOWN_PLOIDY'] = np.nan
    tmp.loc[tmp['WELL'] == 'H12', 'SET'] = True

    tmp.loc[tmp['WELL'] == 'G12', 'HOST_SEX'] = 'FEMALE'
    tmp.loc[tmp['WELL'] == 'G12', 'KNOWN_PLOIDY'] = 'FETAL EUPLOIDY'
    tmp.loc[tmp['WELL'] == 'G12', 'FETAL_SEX'] = 'FEMALE'
    tmp.loc[tmp['WELL'] == 'G12', 'SET'] = True

    tmp.loc[tmp['WELL'] == 'E12', 'HOST_SEX'] = 'MALE'
    tmp.loc[tmp['WELL'] == 'E12', 'KNOWN_PLOIDY'] = 'FETAL EUPLOIDY'
    tmp.loc[tmp['WELL'] == 'E12', 'SET'] = True

    with_run = tmp.loc[tmp['RUN'].isnull() == False]
    rest = tmp.loc[~tmp['RUN'].isnull() == False]

    with_run['RUN_END'] = with_run['RUN'].str[-1]
    with_run['TMP'] = with_run.apply(lambda row: True if (int(row['RUN_END']) % 2 == 0) else False, axis=1)

    with_run.loc[(with_run['WELL'] == 'F12') & (with_run['TMP'] == True), 'HOST_SEX'] = 'FEMALE'
    with_run.loc[(with_run['WELL'] == 'F12') & (with_run['TMP'] == True), 'KNOWN_PLOIDY'] = 'FETAL EUPLOIDY'
    with_run.loc[(with_run['WELL'] == 'F12') & (with_run['TMP'] == True), 'FETAL_SEX'] = 'MALE'
    with_run.loc[(with_run['WELL'] == 'F12') & (with_run['TMP'] == True), 'SET'] = True

    control_info = pd.concat([with_run, rest], axis=0, sort=False)

    return control_info.loc[
        control_info['SET'] == True,
        ['PROPS_ID', 'KNOWN_PLOIDY', 'HOST_SEX', 'FETAL_SEX', 'CONTROL_SAMPLE']].drop_duplicates(subset='PROPS_ID', keep='first')

# paths where data files are
main_data_dir =   Path('/mnt/bfx_projects/nipt_lifecycle/data/wip/metadata')

meta_data =       main_data_dir / 'from_BI'
analytical_data = main_data_dir / 'analytical_data'
other_data =      main_data_dir / 'other_data'

# read in analytical output files
columns_run_file = ['SAMPLE_ID', 'PROPS_ID', 'FLOWCELL', 'PLATE', 'WELL', 'CONTROL_SAMPLE', 'ANALYSIS_DATETIME',
                    'DUPLICATION_RATE', 'READS_TOTAL_ALIGNED_PCT', 'SNP_FETAL_PCT']
early_prod_run = pd.read_csv(analytical_data / 'poly_sample_early_prod.tsv', sep='\t', header=0)
late_prod_run  = pd.read_csv(analytical_data /'poly_sample_late_prod.tsv', sep='\t', header=0)
recent_prod_run = pd.read_csv(analytical_data / 'poly_sample_02042020.tsv', sep='\t', header=0)

# change columns for the latest data
recent_prod_run['PROPS_ID'] = recent_prod_run['SAMPLE_ID']
recent_prod_run.drop(labels=['SAMPLE_ID'], axis=1, inplace=True)
recent_prod_run['SAMPLE_ID'] = recent_prod_run['BFX_RESULT_ID']

# concat together
prod_run_data  = pd.concat([early_prod_run[columns_run_file], late_prod_run[columns_run_file], recent_prod_run[columns_run_file]], axis=0)
prod_run_data.drop_duplicates(subset='SAMPLE_ID', keep='first', inplace=True)

# read in model output files
raw_model31_calls = pd.read_csv(analytical_data / 'calls_model3.1.tsv', sep='\t', header=0)
columns_mod_file  = ['SAMPLE_ID', 'GOF', 'READ_COUNT', 'CHR13_CALL', 'CHR18_CALL', 'CHR21_CALL',
                     'CHRXY_CALL', 'CHR13_PLOIDY', 'CHR18_PLOIDY', 'CHR21_PLOIDY', 'CHRX_PLOIDY', 'CHRY_PLOIDY',
                     'CHR13_TVALUE', 'CHR18_TVALUE', 'CHR21_TVALUE', 'CHRX_TVALUE', 'CHRY_TVALUE', 'CHRY_FETAL_PCT']

model31_calls = pd.concat([raw_model31_calls[columns_mod_file], recent_prod_run[columns_mod_file]], axis=0).drop_duplicates(subset='SAMPLE_ID', keep='first')

# metadata files
version = 'v04'
p_run_data   = pd.read_csv(meta_data / f'progenity_run_data_{version}.tsv', sep='\t', header=0)
a_run_data   = pd.read_csv(meta_data / f'avero_run_data_{version}.tsv', sep='\t', header=0)
run_metadata = pd.concat([p_run_data, a_run_data], axis=0)

p_sample_data   = pd.read_csv(meta_data / f'progenity_sample_data_{version}.tsv', sep='\t', header=0)
a_sample_data   = pd.read_csv(meta_data / f'avero_sample_data_{version}.tsv', sep='\t', header=0)
sample_metadata = pd.concat([p_sample_data, a_sample_data], axis=0)
sample_metadata.drop(labels=['COMPANY'], axis=1, inplace=True) #this is already in run_metadata
sample_metadata.drop_duplicates(subset='SAMPLEID', keep='first', inplace=True)

p_reported = pd.read_csv(meta_data / f'progenity_reported_data_{version}.tsv', sep='\t', header=0)
a_reported = pd.read_csv(meta_data / f'avero_reported_data_{version}.tsv', sep='\t', header=0)
reported   = pd.concat([p_reported, a_reported], axis=0)

exclude_samples = pd.read_csv(other_data / 'exclude_samples.txt', sep='\t', header=0)
outcome_data    = pd.read_csv(other_data / 'outcomes.txt', sep='\t', header=0, dtype={'SID': object})

#aggregate plate and sample counts per flowcell and add to prod data
plate_count = prod_run_data.groupby(['FLOWCELL'])['PLATE'].nunique()
plate_count.rename('PLATE_COUNT',inplace=True)
prod_run_data = prod_run_data.merge(plate_count, how='left', left_on='FLOWCELL', right_on='FLOWCELL')

sample_count = prod_run_data.groupby(['FLOWCELL'])['SAMPLE_ID'].nunique()
sample_count.rename('SAMPLE_COUNT',inplace=True)
prod_run_data = prod_run_data.merge(sample_count, how='left', left_on='FLOWCELL', right_on='FLOWCELL')

# dropping some plates and samples identified by Natalie (see exclude.txt)
exclude_plates = ['PPU70017-1A', 'PPU70018-1B', 'PPU70020-1D']
prod_run_data = prod_run_data.loc[~prod_run_data['PLATE'].isin(exclude_plates)]
prod_run_data = prod_run_data.loc[~prod_run_data['SAMPLE_ID'].isin(exclude_samples['exclude'].to_list())]

# dropping some failed flowcells
failed_fc = ['HMLGHDMXX', 'HMK27DMXX', 'HM3F2DMXX', 'HLV2MDMXX', 'HLV7LDMXX', 'HLVW7DMXX']
prod_run_data = prod_run_data.loc[~prod_run_data['FLOWCELL'].isin(failed_fc)]

# this is a flowcell that was run twice using the same fcid (before and after nipt 2.0 launch). Dropping the first run.
prod_run_data = prod_run_data.loc[~(prod_run_data['FLOWCELL'].isin(['HKF7TDMXX']) & prod_run_data['ANALYSIS_DATETIME'].str.contains('2019-09-08'))]

# join aggregated Run_Project_Poly_9002_run.tsv with the model calls
tmp = prod_run_data.join(model31_calls.set_index('SAMPLE_ID'), sort=False, how='left', on='SAMPLE_ID')

#join in run metadata for progenity and avero
tmp['join_helper'] = tmp['FLOWCELL'] + '_' + tmp['PLATE'] + '_' + tmp['WELL'] + tmp['PROPS_ID']
run_metadata['join_helper'] = run_metadata['FLOWCELL'] + '_' + run_metadata['PLATE'] + '_' + run_metadata['WELL'] + run_metadata['SAMPLEID']
run_metadata.drop(labels=['SAMPLEID', 'PLATE', 'WELL', 'FLOWCELL'], axis=1, inplace=True)
tmp = tmp.join(run_metadata.set_index('join_helper'), sort=False, how='left', on='join_helper')

# merge in sample metadata (demographic data) for progenity and avero
tmp2 = tmp.merge(sample_metadata, how='left', left_on='PROPS_ID', right_on='SAMPLEID')

# putting reported data in a friendly format

euploid_values = ['FETAL EUPLOIDY', 'FETAL EUPLOIDY, FEMALE', 'FETAL EUPLOIDY, MALE', 'CHRY ABSENT', 'CHRY PRESENT']
reported['ID'] = reported['FLOWCELL'] + '_' + reported['SAMPLEID'].apply(str)
unique_sample_ids = reported.drop_duplicates(subset='ID', keep='first')

aneuploid = reported.loc[~reported['RESULT'].isin(euploid_values)]

# remap aneuploid RESULT values
chr = ['13', '18', '21']

for c in chr:
    aneuploid.loc[((aneuploid['OUTCOMENAME'] == f'Chromosome {c}') & (aneuploid['FETALPLOIDY'] == 'Trisomy')), 'RESULT'] = f'TRISOMY {c}'
    aneuploid.loc[((aneuploid['OUTCOMENAME'] == f'Chromosome {c}') & (aneuploid['FETALPLOIDY'] == 'Monosomy')), 'RESULT'] = f'MONOSOMY {c}'

# create calls dataframe
aneuploid_calls = aneuploid.groupby('ID')['RESULT'].apply(', '.join) #join calls if sample has more than one call
euploid_id      = [sid for sid in unique_sample_ids['ID'] if sid not in aneuploid['ID'].tolist()]
euploid_calls   = pd.Series('FETAL EUPLOIDY', index=euploid_id)
calls = euploid_calls.append(aneuploid_calls)

#join calls to manifest
tmp2['join_helper2'] = tmp2['FLOWCELL'] + '_' + tmp2['PROPS_ID']
tmp3 = tmp2.join(calls.rename('REPORTED_PLOIDY'), sort=False, how='left', on='join_helper2')

# now add other columns that are useful

# FETAL SEX
xy = reported.loc[reported['OUTCOMENAME']=='Chromosome XY', ['SAMPLEID','RESULT']]
xy.loc[xy['RESULT'].str.contains('FETAL EUPLOIDY, FEMALE|FETAL X0|FETAL XXX'), 'FETAL_SEX']                   = 'FEMALE'
xy.loc[xy['RESULT'].str.contains('FETAL EUPLOIDY, MALE|FETAL XXY|FETAL XYY|CHRY INDETERMINATE'), 'FETAL_SEX'] = 'MALE'
xy['SAMPLEID'] = xy['SAMPLEID'].astype(str)
xy.drop_duplicates(subset='SAMPLEID', inplace=True, keep='first')
tmp4 = tmp3.merge(xy, how='left', left_on='PROPS_ID', right_on='SAMPLEID')

#HOST SEX
tmp4['HOST_SEX'] = 'FEMALE'
tmp4.loc[tmp4['CONTROL_SAMPLE'] == 'Control', 'HOST_SEX'] = 'MALE' #lab only runs male controls???

#RUN NUMBER
tmp4['RERUN'] = tmp4.duplicated(subset='PROPS_ID', keep=False)

rerun_samples_sorted = tmp4.loc[(tmp4['RERUN'] == True) & (tmp4['CONTROL_SAMPLE'] == 'Test')].sort_values(by=['PROPS_ID', 'ANALYSIS_DATETIME'])
rerun_samples_sorted.loc[:, 'RUN_NUM'] = rerun_samples_sorted.groupby('PROPS_ID').cumcount()

# make it so it isn't 0 indexed :)
rerun_samples_sorted['RUN_NUM'] += 1

# put it together with stuff that was only run once
one_run = tmp4.loc[~((tmp4['RERUN'] == True) & (tmp4['CONTROL_SAMPLE'] == 'Test'))]
one_run.loc[:, 'RUN_NUM'] = 1
tmp5 = pd.concat([one_run, rerun_samples_sorted], axis=0)

# OUTCOME PLOIDY and SEX
tmp6 = tmp5.merge(outcome_data, how='left', left_on='PROPS_ID', right_on='SID')
tmp6['FETAL_SEX'] = np.where(tmp6['OUTCOME_SEX'].isnull() == True, tmp6['FETAL_SEX'], tmp6['OUTCOME_SEX'])

# KNOWN PLOIDY
sercare = ['1902150150', '1902150160', '1902150162']

calls = calls.to_frame()
calls.reset_index(inplace=True)
calls.columns = ['SAMPLE_ID', 'KNOWN_PLOIDY']
calls['PROPS_ID'] = calls['SAMPLE_ID'].str.split('_').str[-1]
calls.drop_duplicates(subset='PROPS_ID', inplace=True, keep='first')
tmp7 = tmp6.merge(calls[['PROPS_ID', 'KNOWN_PLOIDY']], how='left', left_on='PROPS_ID', right_on='PROPS_ID')

tmp7['KNOWN_PLOIDY'] = np.where(tmp7['OUTCOME_PLOIDY'].isnull() == True, tmp7['KNOWN_PLOIDY'], tmp7['OUTCOME_PLOIDY'])
tmp7.loc[tmp7['PROPS_ID'].isin(sercare), 'KNOWN_PLOIDY'] = 'TRISOMY 13, TRISOMY 18, TRISOMY 21'
tmp7.loc[tmp7['CONTROL_SAMPLE'] == 'Control', 'KNOWN_PLOIDY'] = 'FETAL EUPLOIDY'
tmp7.loc[pd.isna(tmp7['FN']), 'FN'] = False

#replace KNOWN_PLOIDY with blank for dubious values/calls
positives = tmp7.loc[(tmp7['KNOWN_PLOIDY'] != 'FETAL EUPLOIDY') & (tmp7['KNOWN_PLOIDY'].isnull() == False)]
outlier_sid = []

for c in chr:
    specific = positives.loc[positives['KNOWN_PLOIDY'].str.contains(f'TRISOMY {c}|MONOSOMY {c}')]
    specific.loc[:, 'ABS'] = abs(specific[f'CHR{c}_TVALUE'])
    outlier = specific.loc[(specific[f'CHR{c}_TVALUE'] < 8) & (specific['SNP_FETAL_PCT'] > 8)]
    outlier_sid += outlier['SAMPLE_ID'].tolist()

print(f'Outliers: {len(outlier_sid)}')

tmp7['KNOWN_PLOIDY'].loc[tmp7['SAMPLE_ID'].isin(outlier_sid)] = ''

# dealing with controls
control_info = control_list(tmp7)

for i, row in control_info.iterrows():
    cond = tmp7['PROPS_ID'] == row['PROPS_ID']
    tmp7.loc[cond, 'HOST_SEX']       = row['HOST_SEX']
    tmp7.loc[cond, 'FETAL_SEX']      = row['FETAL_SEX']
    tmp7.loc[cond, 'CONTROL_SAMPLE'] = row['CONTROL_SAMPLE']
    tmp7.loc[cond, 'KNOWN_PLOIDY']   = row['KNOWN_PLOIDY']

# finally some clean up
tmp7.drop(labels=['join_helper2', 'join_helper', 'SampleId', 'SID', 'RESULT', 'SAMPLEID_x', 'SAMPLEID_y'], axis=1, inplace=True)
tmp7.rename(columns={'OrderId': 'ORDER_ID', 'State': 'STATE', 'PostalCode': 'POSTAL_CODE'}, inplace=True)
tmp7.loc[tmp7['EXTRACTIONINSTRUMENTNAME'] == 'L000461', 'EXTRACTIONINSTRUMENTNAME'] = 'L00461'
tmp7['BMIATTIMEOFDRAW'].replace(0,np.nan, inplace=True)
tmp7['BMIATTIMEOFDRAW'].replace(-1,np.nan, inplace=True)
tmp7.loc[(tmp7['PROPS_ID'].str[0] != 'A') & (tmp7['COMPANY'].isnull()), 'COMPANY'] = 'Progenity'
tmp7.loc[(tmp7['PROPS_ID'].str[0] == 'A') & (tmp7['COMPANY'].isnull()), 'COMPANY'] = 'Avero'
tmp7.rename(columns={'PROPS_ID': 'INDIVIDUAL_ID'}, inplace=True)

tmp7.to_csv('manifest.tsv', sep='\t', index=None)
