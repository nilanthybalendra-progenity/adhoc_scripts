import numpy as np
import os
import pandas as pd

from pathlib import Path

# paths where data files are

main_data_dir =   Path('/mnt/bfx_projects/nipt_lifecycle/data/wip/metadata')

meta_data =       main_data_dir / 'from_BI'
analytical_data = main_data_dir / 'analytical_data'
other_data =      main_data_dir / 'other_data'

# read in analytical output files
columns_run_file = ['SAMPLE_ID', 'PROPS_ID', 'FLOWCELL', 'PLATE', 'WELL', 'CONTROL_SAMPLE', 'ANALYSIS_DATETIME',
                    'DUPLICATION_RATE', 'READS_0MM', 'READS_1MM', 'READS_2MM', 'EFFICIENCY_0MM', 'EFFICIENCY_1MM',
                    'EFFICIENCY_2MM', 'READS_TOTAL_ALIGNED_PCT']
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
columns_mod_file  = ['SAMPLE_ID', 'SNP_FETAL_PCT', 'GOF', 'READ_COUNT', 'CHR13_CALL', 'CHR18_CALL', 'CHR21_CALL',
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

#FETAL SEX
xy = reported.loc[reported['OUTCOMENAME']=='Chromosome XY', ['ID','RESULT']]
xy.loc[xy['RESULT'].str.contains('FETAL EUPLOIDY, FEMALE|FETAL X0|FETAL XXX'), 'FETAL_SEX']                   = 'FEMALE'
xy.loc[xy['RESULT'].str.contains('FETAL EUPLOIDY, MALE|FETAL XXY|FETAL XYY|CHRY INDETERMINATE'), 'FETAL_SEX'] = 'MALE'
tmp4 = tmp3.join(xy.set_index('ID'), sort=False, how='left', on='join_helper2')

#HOST SEX
tmp4['HOST_SEX'] = 'FEMALE'
tmp4.loc[tmp4['CONTROL_SAMPLE'] == 'Control', 'HOST_SEX'] = 'MALE' #lab only runs male controls

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

# OUTCOME PLOIDY
tmp6 = tmp5.merge(outcome_data, how='left', left_on='PROPS_ID', right_on='SID')

# KNOWN PLOIDY
sercare = ['1902150150', '1902150160', '1902150162']
tmp6['KNOWN_PLOIDY'] = tmp6.apply(lambda row: row['OUTCOME_PLOIDY'] if not pd.isna(row['OUTCOME_PLOIDY']) else row['REPORTED_PLOIDY'], axis=1)
tmp6.loc[tmp6['PROPS_ID'].isin(sercare), 'KNOWN_PLOIDY'] = 'TRISOMY 13, TRISOMY 18, TRISOMY 21'
tmp6.loc[tmp6['CONTROL_SAMPLE'] == 'Control', 'KNOWN_PLOIDY'] = 'FETAL EUPLOIDY'
tmp6.loc[pd.isna(tmp6['FN']), 'FN'] = False

#replace KNOWN_PLOIDY with blank for dubious values/calls
positives = tmp6.loc[(tmp6['KNOWN_PLOIDY'] != 'FETAL EUPLOIDY') & (tmp6['KNOWN_PLOIDY'].isnull() == False)]
outlier_sid = []

for c in chr:
    specific = positives.loc[positives['KNOWN_PLOIDY'].str.contains(f'TRISOMY {c}|MONOSOMY {c}')]
    specific.loc[:, 'ABS'] = abs(specific[f'CHR{c}_TVALUE'])
    outlier = specific.loc[(specific[f'CHR{c}_TVALUE'] < 8) & (specific['SNP_FETAL_PCT'] > 8)]
    outlier_sid += outlier['SAMPLE_ID'].tolist()

tmp6['KNOWN_PLOIDY'].loc[tmp6['SAMPLE_ID'].isin(outlier_sid)] = ''

# finally some clean up
tmp6.drop(labels=['join_helper2', 'join_helper','SAMPLEID', 'SampleId', 'SID', 'RESULT'], axis=1, inplace=True)
tmp6.rename(columns={'OrderId': 'ORDER_ID', 'State': 'STATE', 'PostalCode': 'POSTAL_CODE'}, inplace=True)
tmp6.loc[tmp6['EXTRACTIONINSTRUMENTNAME'] == 'L000461', 'EXTRACTIONINSTRUMENTNAME'] = 'L00461'
tmp6['BMIATTIMEOFDRAW'].replace(0,np.nan, inplace=True)
tmp6['BMIATTIMEOFDRAW'].replace(-1,np.nan, inplace=True)
tmp6.loc[(tmp6['PROPS_ID'].str[0] != 'A') & (tmp6['COMPANY'].isnull()), 'COMPANY'] = 'Progenity'
tmp6.loc[(tmp6['PROPS_ID'].str[0] == 'A') & (tmp6['COMPANY'].isnull()), 'COMPANY'] = 'Avero'


bad_controls = pd.read_csv('/mnt/bfx_projects/nipt_lifecycle/analysis/wip/fetal_fraction_model/02_One_off_requests/high_gof_lowt_ctrls.tsv', sep='\t', header=0)

tmp6['BAD_CONTROL'] = tmp6['PROPS_ID'].isin(bad_controls['PROPS_ID'])
tmp6.rename(columns={'PROPS_ID': 'INDIVIDUAL_ID'}, inplace=True)

#fill in KNOWN PLOIDY VALUES
calls = calls.to_frame()
calls.reset_index(inplace=True)
calls.columns = ['SAMPLE_ID', 'CALL']
calls['INDIVIDUAL_ID'] = calls['SAMPLE_ID'].str.split('_').str[-1]
calls.drop_duplicates(subset='INDIVIDUAL_ID', inplace=True)
tmp7 = tmp6.merge(calls[['INDIVIDUAL_ID', 'CALL']], how='left', left_on='INDIVIDUAL_ID', right_on='INDIVIDUAL_ID')


tmp7.to_csv('manifest_new.tsv', sep='\t', index=None)
