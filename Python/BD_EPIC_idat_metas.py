################################################################################################################################################

#Create_BD_EPIC_idat_csv
import pandas as pd
import os
import numpy as np
import GEOparse
from sqlalchemy import create_engine
import json

data_1 = GEOparse.get_GEO(filepath='/Volumes/storage_ext/methylation/_BD_EPIC_idat/GSE112179_idat/GSE112179_family.soft.gz')
data_1df=data_1.phenotype_data
data_1df.reset_index(inplace=True)
data_1df.drop(columns={'title', 'geo_accession', 'status', 'submission_date',
       'last_update_date', 'type', 'channel_count','characteristics_ch1.9.cell.type',
       'organism_ch1', 'taxid_ch1', 'source_name_ch1',
       'characteristics_ch1.1.sampleID', 'characteristics_ch1.2.tissuebank.id',
       'characteristics_ch1.7.tissuebank', 'treatment_protocol_ch1',
       'growth_protocol_ch1', 'molecule_ch1', 'extract_protocol_ch1',
       'label_ch1', 'label_protocol_ch1', 'hyb_protocol', 'scan_protocol',
       'description', 'data_processing', 'platform_id', 'contact_name',
       'contact_email', 'contact_institute', 'contact_address', 'contact_city',
       'contact_state', 'contact_zip/postal_code', 'contact_country',
       'supplementary_file', 'series_id', 'data_row_count'},inplace=True)
data_1df.columns
data_1df.rename(columns={'index':'Sample_Name', 'characteristics_ch1.0.tissue':'Tissue', 'characteristics_ch1.3.age':'Age',
       'characteristics_ch1.4.Sex':'Sex', 'characteristics_ch1.5.race':'Race',
       'characteristics_ch1.6.dist.dx':'Sample_Group', 'characteristics_ch1.8.pmi':'PMI'},inplace=True)
data_1df['Sex'].replace({'Male':'M','Female':'F'},inplace=True)
data_1df['Tissue'].replace({'frontal cortex of the brain':'FC'},inplace=True)
data_1df['Race'].replace({'White':'W','Black':'AF','Hispanic':'H','?':'Other'},inplace=True)
data_1df['Sample_Group'].replace({'Bipolar':'BD_NeuN','Schizophrenia':'SCZ_NeuN','Control':'N_NeuN'},inplace=True)
data_1df['Sample_Well']=pd.NA
data_1df['Pool_ID']=pd.NA
data_1df['Series_ID']='GSE112179'
path = "/Volumes/storage_ext/methylation/_BD_EPIC_idat/GSE112179_idat/GSE112179_RAW"
dir_list=os.listdir(path)
list_df = pd.DataFrame(list(set(dir_list)))
list_df[0] = list_df[0].replace(r'\._', r'', regex=True)
list_df = list_df[0].str.split('_', expand=True)
list_df.drop(3, axis=1, inplace=True)
list_df.rename({0:'Sample_Name',1:'Sentrix_ID',2:'Sentrix_Position'},axis=1,inplace=True)
list_df.drop_duplicates(keep='first', inplace=True)
data_1df=pd.merge(data_1df,list_df,on='Sample_Name')
data_1df.to_csv('/Volumes/storage_ext/methylation/_BD_EPIC_idat/GSE112179_idat/GSE112179_RAW/GSE112179_meta.csv',index=False)


data_2 = GEOparse.get_GEO(filepath='/Volumes/storage_ext/methylation/_BD_EPIC_idat/GSE129428_idat/GSE129428_family.soft.gz')
data_2df=data_2.phenotype_data
data_2df.reset_index(inplace=True)
data_2df.drop(columns={'title', 'geo_accession', 'status',
       'submission_date', 'last_update_date', 'type', 'channel_count',
       'source_name_ch1', 'organism_ch1', 'taxid_ch1','characteristics_ch1.5.tissue ph',
       'molecule_ch1', 'extract_protocol_ch1',
       'label_ch1', 'label_protocol_ch1', 'hyb_protocol', 'scan_protocol',
       'data_processing', 'platform_id', 'contact_name', 'contact_email',
       'contact_phone', 'contact_department', 'contact_institute',
       'contact_address', 'contact_city', 'contact_state',
       'contact_zip/postal_code', 'contact_country', 'supplementary_file',
       'data_row_count'},axis=1,inplace=True)
data_2df.rename(columns={'index':'Sample_Name','characteristics_ch1.0.Sex':'Sex',
                         'characteristics_ch1.1.group':'Sample_Group',
                         'characteristics_ch1.2.age at death':'Age', 
                         'characteristics_ch1.3.race':'Race',
                         'characteristics_ch1.4.tissue':'PMI', 
                         'characteristics_ch1.6.tissue':'Tissue',
                         'series_id':'Series_ID'},inplace=True)
data_2df['Tissue'].replace({'post-mortem hippocampus':'HPC'},inplace=True)
data_2df['Race'].replace({'Caucasian':'W','Black':'AF','Hispanic':'H','?':'Other'},inplace=True)
data_2df['Sex'].replace({'Male':'M','Female':'F'},inplace=True)
data_2df['Sample_Group'].replace({'Bipolar':'BD','Schizophrenia':'SCZ','Control':'N'},inplace=True)
data_2df['Sample_Well']=pd.NA
data_2df['Pool_ID']=pd.NA
split=data_2df['PMI'].str.split(':',expand=True)
split.rename(columns={1:'PMI'},inplace=True)
split.drop(columns={0},inplace=True)
split.reset_index(inplace=True)
data_2df.reset_index(inplace=True)
data_2df=pd.merge(data_2df,split,on='index')
data_2df.drop(columns={'PMI_x'},axis=1,inplace=True)
data_2df.rename(columns={'PMI_y':'PMI'},inplace=True)
path2 = "/Volumes/storage_ext/methylation/_BD_EPIC_idat/GSE129428_idat/GSE129428_RAW"
dir_list2=os.listdir(path2)
list_df2 = pd.DataFrame(list(set(dir_list2)))
list_df2[0] = list_df2[0].replace(r'\._', r'', regex=True)
list_df2 = list_df2[0].str.split('_', expand=True)
list_df2.drop(3, axis=1, inplace=True)
list_df2.rename({0:'Sample_Name',1:'Sentrix_ID',2:'Sentrix_Position'},axis=1,inplace=True)
list_df2.drop_duplicates(keep='first', inplace=True)
data_2df=pd.merge(data_2df,list_df2,on='Sample_Name')
data_2df.drop(columns={'index'},inplace=True)
data_2df = data_2df[['Sample_Name','Tissue','Age','Sex','Race','Sample_Group','PMI','Sample_Well','Pool_ID','Series_ID','Sentrix_ID','Sentrix_Position']]
data_2df.drop(data_2df[data_2df['Sample_Name']=='GSM3712776'].index,inplace=True)
data_2df.to_csv('/Volumes/storage_ext/methylation/_BD_EPIC_idat/GSE129428_idat/GSE129428_RAW/GSE129428_meta.csv',index=False)


data_3 = GEOparse.get_GEO(filepath='/Volumes/storage_ext/methylation/_BD_EPIC_idat/GSE191200_idat/GSE191200_family.soft.gz')
data_3df=data_3.phenotype_data
data_3df.reset_index(inplace=True)
data_3df.drop(columns={'title','characteristics_ch1.0.donor', 'contact_laboratory',
                       'geo_accession', 'status','submission_date', 'last_update_date', 'type', 'channel_count',
                       'source_name_ch1', 'organism_ch1', 'taxid_ch1',
                       'molecule_ch1', 'extract_protocol_ch1',
                       'label_ch1', 'label_protocol_ch1', 'hyb_protocol', 'scan_protocol',
                       'data_processing', 'platform_id', 'contact_name', 'contact_email',
                       'contact_department', 'contact_institute',
                       'contact_address', 'contact_city', 'contact_state',
                       'contact_zip/postal_code', 'contact_country', 'supplementary_file',
                       'data_row_count'},axis=1,inplace=True)
data_3df.rename(columns={'index':'Sample_Name','characteristics_ch1.3.Sex':'Sex',
                         'characteristics_ch1.4.diagnosis':'Sample_Group',
                         'characteristics_ch1.2.age':'Age', 
                         'characteristics_ch1.1.region':'Tissue',
                         'series_id':'Series_ID'},inplace=True)
data_3df['Sample_Group'].replace({'Bipolar':'BD_microglia','Schizophrenis':'SCZ_microglia','Control':'N_microglia'},inplace=True)
data_3df['Sample_Well']=pd.NA
data_3df['Pool_ID']=pd.NA

path3 = "/Volumes/storage_ext/methylation/_BD_EPIC_idat/GSE191200_idat/GSE191200_RAW"
dir_list3=os.listdir(path3)
list_df3 = pd.DataFrame(list(set(dir_list3)))
list_df3[0] = list_df3[0].replace(r'\._', r'', regex=True)
list_df3 = list_df3[0].str.split('_', expand=True)
list_df3.drop(3, axis=1, inplace=True)
list_df3.rename({0:'Sample_Name',1:'Sentrix_ID',2:'Sentrix_Position'},axis=1,inplace=True)
list_df3.drop_duplicates(keep='first', inplace=True)
data_3df=pd.merge(data_3df,list_df3,on='Sample_Name')
data_3df = data_3df[['Sample_Name','Tissue','Age','Sex','Sample_Group','Sample_Well','Pool_ID','Series_ID','Sentrix_ID','Sentrix_Position']]
data_3df_1=data_3df[data_3df['Tissue']=='GFM']
data_3df_1=data_3df_1[data_3df_1['Sample_Group'].isin(['BD_microglia','N_microglia'])]
data_3df_1['Tissue'].replace({'GFM':'FC'},inplace=True)
data_3df_1.reset_index(inplace=True)
data_3df_1.drop(columns={'index'},axis=1,inplace=True)
data_3df_1.to_csv('/Volumes/storage_ext/methylation/_BD_EPIC_idat/GSE191200_idat/GSE191200_RAW/GSE191200_meta.csv',index=False)

data_1df.drop(columns={'Race','PMI'},axis=1,inplace=True)
data_2df.drop(columns={'Race','PMI'},axis=1,inplace=True)
BD_EPIC_meta = pd.concat([data_1df, data_2df,data_3df_1], axis=0) 
BD_EPIC_meta.reset_index(inplace=True)
BD_EPIC_meta.drop(columns={'index'},axis=1,inplace=True)
BD_EPIC_meta.to_csv('/Volumes/storage_ext/methylation/_BD_EPIC_idat/BD_EPIC_meta.csv',index=False)





