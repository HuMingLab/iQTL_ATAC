import numpy as np
import pandas as pd
import re

import sys
import os

CC1=str(sys.argv[1])
CC2=str(sys.argv[2])
CC1_m=str(sys.argv[3])
CC2_m=str(sys.argv[4])
CC1_name=str(sys.argv[5])
CC2_name=str(sys.argv[6])
mm10=str(sys.argv[7])
outputf=str(sys.argv[8])

cc1_uniq=pd.read_csv(CC1, header=None,sep="\t")
cc2_uniq=pd.read_csv(CC2, header=None,sep="\t")
cc1_multi=pd.read_csv(CC1_m, header=None,sep="\t")
cc2_multi=pd.read_csv(CC2_m, header=None,sep="\t")

mm10_loc=pd.read_csv(mm10, header=None,sep="\t")

col2=['chr','mm10_start','mm10_end']
mm10_loc.columns=col2
lst = pd.Series(['chr' + str(i) for i in range(1,20)])
mm10_loc=mm10_loc[(mm10_loc.chr.isin(lst))].drop_duplicates().reset_index(drop=True)


col =['chr','mm10_start','mm10_end','uniq_mat_start','uniq_mat_end','uniq_mat_count']
cola =['chr','mm10_start','mm10_end','uniq_pat_start','uniq_pat_end','uniq_pat_count']
col1 =['chr','mm10_start','mm10_end','common_mat_start','common_mat_end','common_mat_count']
col1a =['chr','mm10_start','mm10_end','common_pat_start','common_pat_end','common_pat_count']

cc1_uniq.columns=col
cc2_uniq.columns=cola

cc1_multi.columns=col1
cc2_multi.columns=col1a

cc1_uniq['chr']=cc1_uniq['chr'].str.replace(CC1_name,'').astype(str)
cc2_uniq['chr']=cc2_uniq['chr'].str.replace(CC2_name,'').astype(str)
cc1_multi['chr']=cc1_multi['chr'].str.replace(CC1_name,'').astype(str)
cc2_multi['chr']=cc2_multi['chr'].str.replace(CC2_name,'').astype(str)

merged_cc1 = cc1_uniq.merge(cc1_multi.drop_duplicates(['chr','mm10_start','mm10_end']), how='inner', on=['chr','mm10_start','mm10_end'])

merged_cc2 = cc2_uniq.merge(cc2_multi.drop_duplicates(['chr','mm10_start','mm10_end']), how='inner', on=['chr','mm10_start','mm10_end'])

peak_merged_cc1 =mm10_loc.merge(merged_cc1, how='left', on=['chr','mm10_start','mm10_end']).fillna(0)
peak_merged_cc2 =mm10_loc.merge(merged_cc2, how='left', on=['chr','mm10_start','mm10_end']).fillna(0)

merged = peak_merged_cc1.merge(peak_merged_cc2.drop_duplicates(['chr','mm10_start','mm10_end']), how='inner', on=['chr','mm10_start','mm10_end'])


merged['total_counts']=merged['uniq_mat_count']+merged['uniq_pat_count']+((merged['common_mat_count']+merged['common_pat_count'])/2)
merged[['mm10_start','mm10_end','uniq_mat_start','uniq_mat_end','uniq_pat_start',
                    'uniq_pat_end','uniq_mat_count','uniq_pat_count','common_mat_count','common_pat_count','total_counts']]=merged[['mm10_start','mm10_end','uniq_mat_start','uniq_mat_end','uniq_pat_start',
                    'uniq_pat_end','uniq_mat_count','uniq_pat_count','common_mat_count','common_pat_count','total_counts']].astype(int)


merged_final=merged[['chr','mm10_start','mm10_end','uniq_mat_start','uniq_mat_end','uniq_pat_start',
                    'uniq_pat_end','uniq_mat_count','uniq_pat_count','common_mat_count','common_pat_count','total_counts']]
#merged_final
merged_final.to_csv(outputf, sep='\t', index=False, header=True)
