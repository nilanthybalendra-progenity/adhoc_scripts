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


def main():
    reformat_progenity = reformat_data('/mnt/bfx_projects/nipt_lifecycle/data/metadata/from_BI/progenity_workflow_data_v12.tsv')
    reformat_avero = reformat_data('/mnt/bfx_projects/nipt_lifecycle/data/metadata/from_BI/avero_workflow_data_v12.tsv')

    manifest = pd.read_csv('manifest_branch_v12.tsv', sep='\t', header=0)

    avero_manifest = manifest.loc[manifest['COMPANY'] == 'Avero']
    progenity_manifest = manifest.loc[manifest['COMPANY'] == 'Progenity']

    reformat_avero.to_csv('reformat.tsv', sep='\t', index=None)
    #join in stuff wooooooooo

    print(len(avero_manifest))
    print(len(progenity_manifest))
    print(len(manifest))

    avero_manifest = avero_manifest.merge(reformat_avero, how='left', left_on='PLATE', right_on='PLATE')
    progenity_manifest = progenity_manifest.merge(reformat_progenity, how='left', left_on='PLATE', right_on='PLATE')

    manifest_final = pd.concat([avero_manifest, progenity_manifest], axis=0, sort=False)

    manifest_final.to_csv('manifest_workflow_time2.tsv', sep='\t', index=None)


    # # add in extraction well data
    # tmp2 = tmp2.merge(ex_well_map[['SAMPLESHEET_WELL', 'EXTRACTION_PLATE_WELL']], how='left', left_on='WELL',
    #                   right_on='SAMPLESHEET_WELL')

if __name__ == "__main__":
    main()



