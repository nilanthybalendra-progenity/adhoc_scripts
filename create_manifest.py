import numpy as np
import os
import pandas as pd

from pathlib import Path

# paths where data files are

main_data_dir =   Path('/mnt/bfx_projects/nipt_lifecycle/data/wip/metadata')

meta_data =       main_data_dir / 'from_BI'
analytical_data = main_data_dir / 'analytical_data'
other_data =      main_data_dir / 'other_data'

# read in files
raw_model31_calls = pd.read_csv(analytical_data / 'calls_model3.1.tsv', sep='\t', header=0)
columns_mod_file  = ['SAMPLE_ID', 'SNP_FETAL_PCT', 'GOF', 'READ_COUNT', 'CHR13_CALL', 'CHR18_CALL', 'CHR21_CALL',
                    'CHRXY_CALL', 'CHR13_TVALUE', 'CHR18_TVALUE', 'CHR21_TVALUE', 'CHRX_TVALUE','CHRY_TVALUE',
                    'CHRY_FETAL_PCT']

model31_calls     = raw_model31_calls[columns_mod_file]
model31_calls.drop_duplicates(subset='SAMPLE_ID', keep='first', inplace=True)

columns_run_file = ['SAMPLE_ID', 'PROPS_ID', 'FLOWCELL', 'PLATE', 'WELL', 'CONTROL_SAMPLE', 'ANALYSIS_DATETIME', 'DUPLICATION_RATE']
early_prod_run = pd.read_csv(analytical_data / 'run_files_early_prod.tsv', sep='\t', header=0)
late_prod_run  = pd.read_csv(analytical_data /'run_files_late_prod.tsv', sep='\t', header=0)
prod_run_data  = pd.concat([early_prod_run[columns_run_file], late_prod_run[columns_run_file]], axis=0)
prod_run_data.drop_duplicates(subset='SAMPLE_ID', keep='first', inplace=True)

p_run_data   = pd.read_csv(meta_data / 'progenity_run_data_v02.tsv.tsv', sep='\t', header=0)
a_run_data   = pd.read_csv(meta_data / 'avero_run_data_v02.tsv.tsv', sep='\t', header=0)
run_metadata = pd.concat([p_run_data, a_run_data], axis=0)

p_sample_data   = pd.read_csv(meta_data / 'progenity_sample_data_v02.tsv.tsv', sep='\t', header=0)
a_sample_data   = pd.read_csv(meta_data / 'avero_sample_data_v02.tsv.tsv', sep='\t', header=0)
sample_metadata = pd.concat([p_sample_data, a_sample_data], axis=0)
sample_metadata.drop(labels=['COMPANY'], axis=1, inplace=True) #this is already in run_metadata
sample_metadata.drop_duplicates(subset='SAMPLEID', keep='first', inplace=True)

p_reported = pd.read_csv(meta_data / 'progenity_reported_data_v02.tsv.tsv', sep='\t', header=0)
a_reported = pd.read_csv(meta_data / 'avero_reported_data_v02.tsv.tsv', sep='\t', header=0)
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

# dropping all controls except males, some plates, and samples identified by Natalie (see exclude.txt)
male_controls  = ['C00259', 'C00260', 'C00261', 'C00262', 'C00263']
exclude_plates = ['PPU70017-1A', 'PPU70018-1B', 'PPU70020-1D']

prod_run_data = prod_run_data.loc[~prod_run_data['PLATE'].isin(exclude_plates)]
prod_run_data = prod_run_data.loc[(prod_run_data['PROPS_ID'].isin(male_controls)) | (prod_run_data['CONTROL_SAMPLE'] == 'Test')]
prod_run_data = prod_run_data.loc[~prod_run_data['SAMPLE_ID'].isin(exclude_samples['exclude'].to_list())]

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

#first drop euploid entries
euploid_values = ['FETAL EUPLOIDY', 'FETAL EUPLOIDY, FEMALE', 'FETAL EUPLOIDY, MALE', 'CHRY ABSENT', 'CHRY PRESENT']
reported['ID'] = reported['FLOWCELL'] + '_' + reported['SAMPLEID'].apply(str)
unique_sample_ids = reported.drop_duplicates(subset='ID', keep='first')

aneuploid = reported.loc[~reported['RESULT'].isin(euploid_values)]
euploid  = [sid for sid in unique_sample_ids['ID'] if sid not in aneuploid['ID']]

# remap aneuploid RESULT values
chr = ['13', '18', '21']

for c in chr:
    aneuploid.loc[(aneuploid['OUTCOMENAME'] == f'Chromosome {c}') & (aneuploid['FETALPLOIDY'] == 'Trisomy'), 'RESULT'] = f'TRISOMY {c}'
    aneuploid.loc[(aneuploid['OUTCOMENAME'] == f'Chromosome {c}') & (aneuploid['FETALPLOIDY'] == 'Monosomy'), 'RESULT'] = f'MONOSOMY {c}'

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
tmp4.loc[tmp4['PROPS_ID'].isin(male_controls), 'HOST_SEX'] = 'MALE'

#RUN NUMBER
tmp4['RERUN'] = tmp4.duplicated(subset='PROPS_ID', keep=False)

rerun_samples_sorted = tmp4.loc[(tmp4['RERUN'] == True) & (tmp4['CONTROL_SAMPLE'] == 'Test')].sort_values(by=['PROPS_ID', 'ANALYSIS_DATETIME'])
rerun_samples_sorted.loc[:, 'RUN_NUM'] = rerun_samples_sorted.groupby('PROPS_ID').cumcount()

# map count to actual values
mapping = {0: '1', 1: '2', 2: '>2'}
rerun_samples_sorted.replace({'RUN_NUM': mapping}, inplace=True)

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
tmp6.loc[pd.isna(tmp6['FN']), 'FN'] = False

#replace KNOWN_PLOIDY with blank for dubious values/calls
positives = tmp6.loc[(tmp6['KNOWN_PLOIDY'] != 'FETAL EUPLOIDY') & (tmp6['KNOWN_PLOIDY'].isnull() == False)]
outlier_sid = []

for c in chr:
    specific = positives.loc[positives['KNOWN_PLOIDY'].str.contains(f'TRISOMY {c}|MONOSOMY {c}')]
    specific['ABS'] = abs(specific[f'CHR{c}_TVALUE'])
    outlier = specific.loc[(specific[f'CHR{c}_TVALUE'] < 8) & (specific['SNP_FETAL_PCT'] > 8)]
    outlier_sid += outlier['SAMPLE_ID'].tolist()

tmp6['KNOWN_PLOIDY'].loc[tmp6['SAMPLE_ID'].isin(outlier_sid)] = ''

# finally some clean up
tmp6.drop(labels=['join_helper2', 'join_helper','SAMPLEID', 'SID', 'RESULT'], axis=1, inplace=True)
tmp6.loc[tmp6['EXTRACTIONINSTRUMENTNAME'] == 'L000461', 'EXTRACTIONINSTRUMENTNAME'] = 'L00461'
tmp6['BMIATTIMEOFDRAW'].replace(0,np.nan, inplace=True)
tmp6['BMIATTIMEOFDRAW'].replace(-1,np.nan, inplace=True)
tmp6.loc[(tmp6['PROPS_ID'].str[0] != 'A') & (tmp6['COMPANY'].isnull()), 'COMPANY'] = 'Progenity'
tmp6.loc[(tmp6['PROPS_ID'].str[0] == 'A') & (tmp6['COMPANY'].isnull()), 'COMPANY'] = 'Avero'

tmp6.to_csv('manifest.tsv', sep='\t', index=None)
