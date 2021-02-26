import numpy as np
import os
import json
import pandas as pd

from pathlib import Path


def control_list(data):
    '''Set '''
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
main_data_dir =   Path('/mnt/bfx_projects/nipt_lifecycle/data/metadata')

meta_data =       main_data_dir / 'from_BI'
analytical_data = main_data_dir / 'analytical_data'
clinical_data   = main_data_dir / 'from_clinical'
other_data =      main_data_dir / 'other_data'

# read in analytical output files
sample_cols = ['SAMPLE_ID', 'PROPS_ID', 'FLOWCELL', 'PLATE', 'WELL', 'CONTROL_SAMPLE', 'ANALYSIS_DATETIME',
                    'DUPLICATION_RATE', 'READS_TOTAL_ALIGNED_PCT', 'READS_TOTAL_COUNT','SNP_FETAL_PCT']

#samples.tsv
early_prod_sample  = pd.read_csv(analytical_data / 'poly_sample_early_prod.tsv', sep='\t', header=0)
late_prod_sample   = pd.read_csv(analytical_data /'poly_sample_02262021.tsv', sep='\t', header=0)
other_samples      = pd.read_csv(analytical_data /'HLV5WDMXX_const_T21_sample_fixed.tsv', sep='\t', header=0)
avero_validation   = pd.read_csv(analytical_data / 'avero_validation_samples.tsv', sep='\t', header=0)

avero_validation_fc = ['HJVYFDMXX', 'HKJMNDMXX', 'HKJYCDMXX', 'HKJYMDMXX']
other_fc = ['HLV5WDMXX'] #consitutional T21 flowcell

recent_prod_sample = pd.concat([late_prod_sample, avero_validation, other_samples], axis=0)

# change columns for the latest data
recent_prod_sample['PROPS_ID'] = recent_prod_sample['SAMPLE_ID']
recent_prod_sample.drop(labels=['SAMPLE_ID'], axis=1, inplace=True)
recent_prod_sample['SAMPLE_ID'] = recent_prod_sample['BFX_RESULT_ID']

# concat together all sample tsv
prod_sample_data  = pd.concat([early_prod_sample[sample_cols], recent_prod_sample[sample_cols]], axis=0)
prod_sample_data['INDEX_PLATE'] = prod_sample_data['PLATE'].str.split('-').str[-1]
prod_sample_data.drop_duplicates(subset='SAMPLE_ID', keep='first', inplace=True)

# run.tsv
early_prod_run  = pd.read_csv(analytical_data / 'poly_run_early_prod.tsv', sep='\t', header=0)
recent_prod_run = pd.read_csv(analytical_data / 'poly_run_02262021.tsv', sep='\t', header=0)
other_run       = pd.read_csv(analytical_data / 'HLV5WDMXX_const_T21_run.tsv', sep='\t', header=0)

prod_run_data = pd.concat([early_prod_run, recent_prod_run, other_run], axis=0)
prod_run_data.drop_duplicates(subset='FCID', keep='first', inplace=True)

# read in model info (need to use a specific model file for pre 2.0, as a different model was used, for the rest, the tsvs can be used)
model31_calls_for_pre_2 = pd.read_csv(analytical_data / 'calls_model3.1.tsv', sep='\t', header=0)

mod_cols  = ['SAMPLE_ID', 'GOF', 'READ_COUNT', 'CHR13_CALL', 'CHR18_CALL', 'CHR21_CALL',
                     'CHRXY_CALL', 'CHR13_PLOIDY', 'CHR18_PLOIDY', 'CHR21_PLOIDY', 'CHRX_PLOIDY', 'CHRY_PLOIDY',
                     'CHR13_TVALUE', 'CHR18_TVALUE', 'CHR21_TVALUE', 'CHRX_TVALUE', 'CHRY_TVALUE', 'CHR13_FETAL_PCT', 'CHR18_FETAL_PCT','CHR21_FETAL_PCT', 'CHRY_FETAL_PCT']

model31_calls = pd.concat([model31_calls_for_pre_2[mod_cols], recent_prod_sample[mod_cols]], axis=0).drop_duplicates(subset='SAMPLE_ID', keep='first')

# also include t-values at time of production
raw_model_cols = ['SAMPLE_ID', 'CHR13_TVALUE', 'CHR18_TVALUE', 'CHR21_TVALUE', 'CHRX_TVALUE', 'CHRY_TVALUE']
raw_model_calls = pd.concat([early_prod_sample[raw_model_cols], recent_prod_sample[raw_model_cols]], axis=0).drop_duplicates(subset='SAMPLE_ID', keep='first')
raw_model_calls.rename(columns={'CHR13_TVALUE': 'CHR13_TVALUE_PROD',
                                'CHR18_TVALUE': 'CHR18_TVALUE_PROD',
                                'CHR21_TVALUE': 'CHR21_TVALUE_PROD',
                                'CHRX_TVALUE': 'CHRX_TVALUE_PROD',
                                'CHRY_TVALUE': 'CHRY_TVALUE_PROD'}, inplace=True)

# metadata files
version = 'v18' #BI extract version number
p_run_data   = pd.read_csv(meta_data / f'progenity_run_data_{version}.tsv', sep='\t', header=0)
a_run_data   = pd.read_csv(meta_data / f'avero_run_data_{version}.tsv', sep='\t', header=0)
run_metadata = pd.concat([p_run_data, a_run_data], axis=0)

p_sample_data   = pd.read_csv(meta_data / f'progenity_sample_data_{version}.tsv', sep='\t', header=0)
a_sample_data   = pd.read_csv(meta_data / f'avero_sample_data_{version}.tsv', sep='\t', header=0)
sample_metadata = pd.concat([p_sample_data, a_sample_data], axis=0)

sample_metadata.drop(labels=['COMPANY', 'SAMPLEID.1'], axis=1, inplace=True) #compnay is already in run_metadata, there are two SAMPLEID columns
sample_metadata.drop_duplicates(subset='SAMPLEID', keep='first', inplace=True)

#fix weird STATE values
with open("states.json") as f:
    state_abbv = json.load(f)

sample_metadata['STATE'] = sample_metadata.apply(
    lambda row: row['STATE'] if row['STATE'] in state_abbv.values()
    else state_abbv[row['STATE']] if row['STATE'] in state_abbv.keys()
    else row['POSTALCODE'] if str(row['POSTALCODE']) in state_abbv.values()
    else 'CA' if str(row['STATE'])[:2] == 'CA'
    else np.nan, axis=1)

# the reported file for progenity was getting too large to output, so now partial extracts are being used
p_reported_partial_5 = pd.read_csv(meta_data / f'progenity_reported_data_PARTIAL_{version}.tsv', sep='\t', header=0)
p_reported_partial_4 = pd.read_csv(meta_data / f'progenity_reported_data_PARTIAL_v17.tsv', sep='\t', header=0)
p_reported_partial_3 = pd.read_csv(meta_data / f'progenity_reported_data_PARTIAL_v16.tsv', sep='\t', header=0)
p_reported_partial_2 = pd.read_csv(meta_data / f'progenity_reported_data_PARTIAL_v15.tsv', sep='\t', header=0)
p_reported_partial_1 = pd.read_csv(meta_data / f'progenity_reported_data_PARTIAL_v14.tsv', sep='\t', header=0)
p_reported_rest = pd.read_csv(meta_data / f'progenity_reported_data_v12.tsv', sep='\t', header=0)
p_reported = pd.concat([p_reported_partial_1, p_reported_partial_2, p_reported_partial_3, p_reported_partial_4, p_reported_partial_5, p_reported_rest], axis=0).drop_duplicates(keep='last')

a_reported = pd.read_csv(meta_data / f'avero_reported_data_{version}.tsv', sep='\t', header=0)

# a_reported_partial = pd.read_csv(meta_data / f'avero_reported_data_PARTIAL_{version}.tsv', sep='\t', header=0)
# a_reported_rest = pd.read_csv(meta_data / f'avero_reported_data_v13.tsv', sep='\t', header=0)
# a_reported = pd.read_csv(meta_data / f'avero_reported_data_v13.tsv', sep='\t', header=0)

reported   = pd.concat([p_reported, a_reported], axis=0)
reported   = reported.loc[reported['OUTCOMENAME'] != 'Chromosome XY - Twin'] #because we don't care about twin calls

exclude_samples  = pd.read_csv(other_data / 'exclude_samples.txt', sep='\t', header=0)
outcome_data     = pd.read_csv(other_data / 'outcomes.txt', sep='\t', header=0, dtype={'SID': object})
ex_well_map      = pd.read_csv(other_data / 'NIPT_well_mapping.tsv', sep='\t', header=0)
validation_truth = pd.read_csv(clinical_data / 'avero_validation_truth.tsv', sep='\t', header=0)

#aggregate plate and sample counts per flowcell and add to prod data
plate_count = prod_sample_data.groupby(['FLOWCELL'])['PLATE'].nunique()
plate_count.rename('PLATE_COUNT',inplace=True)
prod_sample_data = prod_sample_data.merge(plate_count, how='left', left_on='FLOWCELL', right_on='FLOWCELL')

sample_count = prod_sample_data.groupby(['FLOWCELL'])['SAMPLE_ID'].nunique()
sample_count.rename('SAMPLE_COUNT',inplace=True)
prod_sample_data = prod_sample_data.merge(sample_count, how='left', left_on='FLOWCELL', right_on='FLOWCELL')

# dropping some plates and samples identified by Natalie (see exclude.txt)
exclude_plates = ['PPU70017-1A', 'PPU70018-1B', 'PPU70020-1D']
prod_sample_data = prod_sample_data.loc[~prod_sample_data['PLATE'].isin(exclude_plates)]
prod_sample_data = prod_sample_data.loc[~prod_sample_data['SAMPLE_ID'].isin(exclude_samples['exclude'].to_list())]

# dropping some failed flowcells
failed_fc = ['HMLGHDMXX', 'HMK27DMXX', 'HM3F2DMXX', 'HLV2MDMXX', 'HLV7LDMXX', 'HLVW7DMXX', 'HN3KGDMXX']
prod_sample_data = prod_sample_data.loc[~prod_sample_data['FLOWCELL'].isin(failed_fc)]

# this is a flowcell that was run twice using the same fcid (before and after nipt 2.0 launch). Dropping the first run.
prod_sample_data = prod_sample_data.loc[~(prod_sample_data['FLOWCELL'].isin(['HKF7TDMXX']) & prod_sample_data['ANALYSIS_DATETIME'].str.contains('2019-09-08'))]

# start combining datasets

# join aggregated samples.tsv with the model calls
tmp = prod_sample_data.join(model31_calls.set_index('SAMPLE_ID'), sort=False, how='left', on='SAMPLE_ID')
tmp = tmp.join(raw_model_calls.set_index('SAMPLE_ID'), sort=False, how='left', on='SAMPLE_ID')

# add in run.tsv data
tmp = tmp.merge(prod_run_data[['FCID', 'RUN_PHIX_ALIGN_PCT', 'YIELD', 'WORKFLOW_VERSION']], how='left', left_on='FLOWCELL', right_on='FCID')

#join in run metadata for progenity and avero
tmp['join_helper'] = tmp['FLOWCELL'] + '_' + tmp['PLATE'] + '_' + tmp['WELL'] + tmp['PROPS_ID']
run_metadata['join_helper'] = run_metadata['FLOWCELL'] + '_' + run_metadata['PLATE'] + '_' + run_metadata['WELL'] + run_metadata['SAMPLEID']
run_metadata.drop(labels=['SAMPLEID', 'PLATE', 'WELL', 'FLOWCELL'], axis=1, inplace=True)
tmp = tmp.join(run_metadata.set_index('join_helper'), sort=False, how='left', on='join_helper')

# merge in sample metadata (demographic data) for progenity and avero
tmp2 = tmp.merge(sample_metadata, how='left', left_on='PROPS_ID', right_on='SAMPLEID')

#add in extraction well data
tmp2 = tmp2.merge(ex_well_map[['SAMPLESHEET_WELL', 'EXTRACTION_PLATE_WELL']], how='left', left_on='WELL', right_on='SAMPLESHEET_WELL')

# putting reported data in a friendly format
euploid_values = ['FETAL EUPLOIDY', 'FETAL EUPLOIDY, FEMALE', 'FETAL EUPLOIDY, MALE', 'CHRY ABSENT', 'CHRY PRESENT']
reported['ID'] = reported['FLOWCELL'] + '_' + reported['SAMPLEID'].apply(str)
unique_sample_ids = reported.drop_duplicates(subset='ID', keep='first')

aneuploid = reported.loc[~reported['RESULT'].isin(euploid_values)]

#this will speed things up (still working on this)
if os.path.isfile(other_data / f'aneuploid_calls_{version}.tsv'):
    print('Reformated aneuploid calls file exists already!')
    aneuploid_calls = pd.read_csv(other_data / f'aneuploid_calls_{version}.tsv', sep='\t', header=0, index_col=0)

else:
    # remap aneuploid RESULT values
    chr = ['13', '18', '21']
    for c in chr:
        aneuploid.loc[((aneuploid['OUTCOMENAME'] == f'Chromosome {c}') & (aneuploid['FETALPLOIDY'] == 'Trisomy')), 'RESULT'] = f'TRISOMY {c}'
        aneuploid.loc[((aneuploid['OUTCOMENAME'] == f'Chromosome {c}') & (aneuploid['FETALPLOIDY'] == 'Monosomy')), 'RESULT'] = f'MONOSOMY {c}'

    # create calls dataframe
    aneuploid_calls = aneuploid.groupby('ID')['RESULT'].apply(', '.join) #join calls if sample has more than one call

euploid_id      = [sid for sid in unique_sample_ids['ID'] if sid not in aneuploid['ID'].tolist()] #everything not in the above list but in the reported dataframe is euploid
euploid_calls   = pd.Series('FETAL EUPLOIDY', index=euploid_id)
calls = euploid_calls.append(aneuploid_calls)

# join in reported values
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
tmp4.loc[tmp4['CONTROL_SAMPLE'] == 'Control', 'HOST_SEX'] = '' #fill in knowns later

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

#flag outliers
positives = tmp7.loc[(tmp7['KNOWN_PLOIDY'] != 'FETAL EUPLOIDY') & (tmp7['KNOWN_PLOIDY'].isnull() == False)]
outlier_sid = []

for c in chr:
    specific = positives.loc[positives['KNOWN_PLOIDY'].str.contains(f'TRISOMY {c}|MONOSOMY {c}')]
    specific.loc[:, 'ABS'] = abs(specific[f'CHR{c}_TVALUE'])
    outlier = specific.loc[(specific[f'CHR{c}_TVALUE'] < 8) & (specific['SNP_FETAL_PCT'] > 8)]
    outlier_sid += outlier['SAMPLE_ID'].tolist()

print(f'Outliers: {len(outlier_sid)}')
tmp7['IS_OUTLIER'] = False
tmp7.loc[tmp7['SAMPLE_ID'].isin(outlier_sid), 'IS_OUTLIER'] = True

# dealing with controls, based on info from clinical
control_info = control_list(tmp7)

for i, row in control_info.iterrows():
    cond = tmp7['PROPS_ID'] == row['PROPS_ID']
    tmp7.loc[cond, 'HOST_SEX']       = row['HOST_SEX']
    tmp7.loc[cond, 'FETAL_SEX']      = row['FETAL_SEX']
    tmp7.loc[cond, 'CONTROL_SAMPLE'] = row['CONTROL_SAMPLE']
    tmp7.loc[cond, 'KNOWN_PLOIDY']   = row['KNOWN_PLOIDY']

# finally some clean up
tmp7.to_csv('tmp7.tsv', sep='\t', index=None)
tmp7.drop(labels=['FCID', 'join_helper2', 'join_helper', 'SID', 'RESULT', 'SAMPLEID_x', 'SAMPLEID_y', 'SAMPLESHEET_WELL'], axis=1, inplace=True)
tmp7.rename(columns={'CONTROL_SAMPLE': 'SAMPLE_TYPE', 'SAMPLETYPE': 'DNA_SOURCE'}, inplace=True)
tmp7.loc[tmp7['EXTRACTIONINSTRUMENTNAME'] == 'L000461', 'EXTRACTIONINSTRUMENTNAME'] = 'L00461'
tmp7['BMIATTIMEOFDRAW'].replace(0,np.nan, inplace=True)
tmp7['BMIATTIMEOFDRAW'].replace(-1,np.nan, inplace=True)
tmp7.loc[(tmp7['PROPS_ID'].str[0] != 'A') & (tmp7['COMPANY'].isnull()), 'COMPANY'] = 'Progenity'
tmp7.loc[(tmp7['PROPS_ID'].str[0] == 'A') & (tmp7['COMPANY'].isnull()), 'COMPANY'] = 'Avero'
tmp7.rename(columns={'PROPS_ID': 'INDIVIDUAL_ID'}, inplace=True)

# clean up failures
failed = tmp7.loc[tmp7['REPORTED_PLOIDY'].str.contains('Fail') | tmp7['REPORTED_PLOIDY'].str.contains('FAIL')]

not_failed_reported = tmp7.loc[~(tmp7['REPORTED_PLOIDY'].str.contains('Fail') | tmp7['REPORTED_PLOIDY'].str.contains('FAIL')) & (tmp7['REPORTED_PLOIDY'].isnull() == False)]
not_failed_or_reported = tmp7.loc[~(tmp7['REPORTED_PLOIDY'].str.contains('Fail') | tmp7['REPORTED_PLOIDY'].str.contains('FAIL')) & (tmp7['REPORTED_PLOIDY'].isnull() == True)]
failed['REPORTED_PLOIDY'] = failed['REPORTED_PLOIDY'].str.split(',').str[0] # get rid of duplicate comma sep entries for failures
tmp8 = pd.concat([failed, not_failed_reported, not_failed_or_reported], axis=0, sort=False)

failed = tmp8.loc[tmp8['KNOWN_PLOIDY'].str.contains('Fail') | tmp8['KNOWN_PLOIDY'].str.contains('FAIL')]
failed['IS_FAIL'] = 'TRUE'
not_failed = tmp8.loc[~(tmp8['KNOWN_PLOIDY'].str.contains('Fail') | tmp8['KNOWN_PLOIDY'].str.contains('FAIL'))]
not_failed['IS_FAIL'] = 'FALSE'
failed['KNOWN_PLOIDY'] = failed['KNOWN_PLOIDY'].str.split(',').str[0] # get rid of duplicate comma sep entries for failures
tmp9 = pd.concat([failed, not_failed], axis=0, sort=False)

# fixing weird user entry
tmp9.loc[tmp9['PLATESETUPUSERNAME'] == 'matthew.o&#39;hara', 'PLATESETUPUSERNAME']           = 'matthew.ohara'
tmp9.loc[tmp9['TARGETEDCAPTUREUSERNAME'] == 'matthew.o&#39;hara', 'TARGETEDCAPTUREUSERNAME'] = 'matthew.ohara'
tmp9.loc[tmp9['INDEXINGPCRUSERNAME'] == 'matthew.o&#39;hara', 'INDEXINGPCRUSERNAME']         = 'matthew.ohara'

# add truth info for validation flowcells
tmp9['SOURCE'] = 'PRODUCTION'
tmp9.loc[tmp9['FLOWCELL'].isin(avero_validation_fc), 'SOURCE'] = 'VALIDATION'
tmp9.loc[tmp9['FLOWCELL'].isin(other_fc), 'SOURCE'] = 'EXPERIMENT'

clinical   = tmp9.loc[tmp9['SOURCE'] == 'PRODUCTION']
validation = tmp9.loc[tmp9['SOURCE'] == 'VALIDATION']
other      = tmp9.loc[tmp9['SOURCE'] == 'EXPERIMENT']

#these are the constitutional T21 samples
other.loc[other['INDIVIDUAL_ID'].str.contains('p'), 'KNOWN_PLOIDY'] = 'TRISOMY 21'
# other.loc[other['INDIVIDUAL_ID'].str.contains('p'), 'SNP_FETAL_PCT'] = 1.0

validation.drop(labels=['INDIVIDUAL_ID', 'FLOWCELL', 'PREGNANCYTYPE', 'FETAL_SEX', 'KNOWN_PLOIDY'], axis=1, inplace=True)

# drop remnants
validation = validation.loc[validation['SAMPLE_ID'].isin(validation_truth['SAMPLE_ID'])]
val_truth = validation.set_index('SAMPLE_ID').join(validation_truth.set_index('SAMPLE_ID'), sort=False, how='left', on='SAMPLE_ID')
val_truth.reset_index(inplace=True)
full = pd.concat([clinical, val_truth, other], axis=0, sort=False)

full.to_csv('manifest_branch_adhoc.tsv', sep='\t', index=None)
