import pandas as pd


def reformat_data(file_path):

    time_data = pd.read_csv(file_path, sep='\t', header=0)
    time_data = time_data[['PLATEBARCODE', 'PLATESTARTOFSTEP', 'STEPNAME']]

    steps = pd.read_csv('/mnt/bfx_projects/nipt_lifecycle/data/metadata/from_BI/workflow_steps.tsv', sep='\t', header=0)
    print(steps)

    reformat = pd.DataFrame()

    reformat['PLATE'] = time_data['PLATEBARCODE'].drop_duplicates()

    for s in steps['STEP']:
        subset = time_data.loc[time_data['STEPNAME'] == s]
        col_name = s.replace(" - ", "_").replace(' ', '_')
        print(col_name)
        subset.rename(columns={'PLATEBARCODE':'PLATE', 'PLATESTARTOFSTEP': f'{col_name.upper()}_TIME'}, inplace=True)

        subset.drop(labels=['STEPNAME'], axis=1, inplace=True)
        reformat = reformat.join(subset.set_index('PLATE'), sort=False, how='left', on='PLATE')
        reformat = reformat.drop_duplicates(subset='PLATE')

    return reformat


def reformat_lot_data(file_path):

    lot_data = pd.read_csv(file_path, sep='\t', header=0)
    lot_data = lot_data[['PLATEBARCODE', 'FIELD_NAME', 'FIELD_DATA']]

    field_name = ['assay_target_plate_barcode', 'assay_target_plate_lot', 'assay_target_plate_lot_expiration_date', 
    'extension_ligation_plate_barcode', 'extension_ligation_plate_lot', 'extension_ligation_plate_lot_expiration_date', 
    'index_plate_barcode', 'index_plate_barcode_field', 'index_plate_lot', 'index_plate_lot_expiration_date']

    reformat = pd.DataFrame()
    reformat['PLATE'] = lot_data['PLATEBARCODE'].drop_duplicates()

    for s in field_name:
        subset = lot_data.loc[lot_data['FIELD_NAME'] == s]
      
        subset.rename(columns={'PLATEBARCODE':'PLATE', 'FIELD_DATA': f'{s.upper()}'}, inplace=True)
        subset.drop(labels=['FIELD_NAME'], axis=1, inplace=True)
        reformat = reformat.join(subset.set_index('PLATE'), sort=False, how='left', on='PLATE')
        reformat = reformat.drop_duplicates(subset='PLATE')

    return reformat


def main():
    manifest = pd.read_csv('manifest_branch_v13.tsv', sep='\t', header=0)
    avero_manifest = manifest.loc[manifest['COMPANY'] == 'Avero']
    progenity_manifest = manifest.loc[manifest['COMPANY'] == 'Progenity']

    #reformat_avero.to_csv('reformat.tsv', sep='\t', index=None)
    #join in stuff wooooooooo

    # print(len(avero_manifest))
    # print(len(progenity_manifest))
    # print(len(manifest))

    print('reformat time data')

    reformat_progenity_time = reformat_data('/mnt/bfx_projects/nipt_lifecycle/data/metadata/from_BI/progenity_workflow_data_v13.tsv')
    reformat_avero_time = reformat_data('/mnt/bfx_projects/nipt_lifecycle/data/metadata/from_BI/avero_workflow_data_v13.tsv')

    print('reformat lot data')
    reformat_progenity_lot = reformat_lot_data('/mnt/bfx_projects/nipt_lifecycle/data/metadata/from_BI/progenity_lot_data_v13.tsv')
    reformat_avero_lot = reformat_lot_data('/mnt/bfx_projects/nipt_lifecycle/data/metadata/from_BI/avero_lot_data_v13.tsv')


    avero_manifest = avero_manifest.merge(reformat_avero_time, how='left', left_on='PLATE', right_on='PLATE')
    avero_manifest = avero_manifest.merge(reformat_avero_lot, how='left', left_on='PLATE', right_on='PLATE')

    progenity_manifest = progenity_manifest.merge(reformat_progenity_time, how='left', left_on='PLATE', right_on='PLATE')
    progenity_manifest = progenity_manifest.merge(reformat_progenity_lot, how='left', left_on='PLATE', right_on='PLATE')


    # do the same for 

    manifest_final = pd.concat([avero_manifest, progenity_manifest], axis=0, sort=False)

    manifest_final.to_csv('manifest_workflow_time_lot_v13.tsv', sep='\t', index=None)


    # # add in extraction well data
    # tmp2 = tmp2.merge(ex_well_map[['SAMPLESHEET_WELL', 'EXTRACTION_PLATE_WELL']], how='left', left_on='WELL',
    #                   right_on='SAMPLESHEET_WELL')

if __name__ == "__main__":
    main()



