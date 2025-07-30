file_Met = '../data/THSBC_Met.csv'
file_ITS = '../data/THSBC_ITS.csv'
file_16S = '../data/THSBC_16S.csv'
file_Taxonomy = '../data/THSBC_Tax.csv'
file_MetaCycPath = '../data/THSBC_MetaCycPath.csv'
file_Hepatic = '../data/THSBC_Hepatic.csv'
file_Renal = '../data/THSBC_Renal.csv'

west_Met = '../data/WestLake_Met.csv'
west_16S = '../data/WestLake_16S.csv'
west_ITS = '../data/WestLake_ITS.csv'
west_Taxonomy = '../data/WestLake_Taxonomy.csv'

HBC_Met = '../data//HBC_Met.csv'
HBC_ITS = '../data//HBC_ITS.csv'
HBC_Taxonomy = '../data//HBC_Tax.csv'
HBC_MetaCyc = '../data//HBC_MetaCycPath.csv'

HBC_list=list('Met'=HBC_Met, 'ITS'=HBC_ITS, 'Taxonomy'=HBC_Taxonomy, 'MetaCyc'=HBC_MetaCyc)
west_list=list('Met'=west_Met, '16S'=west_16S, 'ITS'=west_ITS, 'Taxonomy'=west_Taxonomy)
file_list=list('Met'=file_Met, 'ITS'=file_ITS, '16S'=file_16S, 'Taxonomy'=file_Taxonomy, 'MetaCycPath'=file_MetaCycPath, 'Hepatic'=file_Hepatic, 'Renal'=file_Renal)

EPDS_cov_names = c('age','BMI_prep','parity','smk','drk','edu')
COM_cov_names_wk = c('age','BMI_prep','parity','smk','drk','edu', 'wk')
HDP_cov_names_batch = c('age','BMI_prep','parity','smk','drk','edu', 'wk', 'batch')