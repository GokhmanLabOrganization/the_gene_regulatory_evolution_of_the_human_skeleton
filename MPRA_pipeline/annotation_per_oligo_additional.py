import pandas as pd
import sys
import pybedtools
from pybedtools import BedTool
import numpy as np
import os

os.chdir('/home/labs/davidgo/Collaboration/Lab_Tools/liftOver/no_GUI')
from liftOver import initializer as initializer
from liftOver import liftOver_df as liftOver

os.chdir('/home/labs/davidgo/Collaboration/humanMPRA')


cells = sys.argv[1]

anno_v3 = pd.read_csv(f'./top_candidates/{cells}/humanMPRA_annotations_v3.csv', header=0)
print('Number of oligos in hMPRA:', len(anno_v3))

anno_v3 = anno_v3.drop_duplicates(subset=["oligo"], keep = "first")
print('Number of oligos in hMPRA:', len(anno_v3))

anno_v3_bed = anno_v3[['chromosome', 'start', 'end', 'id']].copy()
anno_v3_bed = anno_v3_bed[['chromosome', 'start', 'end', 'id']]

anno_v3_bed = pybedtools.BedTool.from_dataframe(anno_v3_bed)
anno_v3_bed_sort = anno_v3_bed.sort()

# SCREEN v2 element type
screen_version='V2'
screen_file_hg19 = '/home/labs/davidgo/Collaboration/GenomeAnnotation/Human/SCREEN/V2/GRCh37-lifted-cCREs.bed'
screen_annotation_file='/home/labs/davidgo/Collaboration/GenomeAnnotation/Human/SCREEN/V2/GRCh38-cCREs.bed'

activity_threshold = 4.24 # This values is based on the MAD significance.

screen_hg19_df=pd.read_csv(screen_file_hg19, sep='\t',header=None)
screen_annotation_df=pd.read_csv(screen_annotation_file, sep='\t',header=None)
screen_annotation_df.drop(columns=[0,1,2,4], inplace=True)
screen_df = pd.merge(screen_hg19_df,screen_annotation_df,left_on=3,right_on=3)
screen_df.drop(columns=['5_x'], inplace=True)
screen_df.rename(columns={"5_y": "regulatoryClass"},inplace=True)
screen_df["regulatoryClass"] = screen_df["regulatoryClass"].str.replace(',CTCF-bound','')
screen_df.drop(screen_df[screen_df.regulatoryClass =='CTCF-only'].index, inplace=True)

screen_df = screen_df.sort_values([0,1], ascending = [True, True])
screen_bed=pybedtools.BedTool.from_dataframe(screen_df)

#intersect SCREEN v2 with the oligos

classes_df = anno_v3_bed_sort.intersect(screen_bed, wao=True).to_dataframe()
classes_df = classes_df[['name', 'blockCount']]
classes_df = classes_df.groupby('name', as_index=False).agg({'blockCount': list})
classes_df.rename(columns={'blockCount': 'CRE_class'}, inplace=True)

# Add CRE class per locus

anno_v3merged = pd.merge(anno_v3, classes_df[['CRE_class', 'name']], left_on='id', right_on = "name", how='left')
anno_v3merged.drop(columns = ["name"],inplace=True)

anno_v3merged['SCREEN2_class_strict'] = (
    anno_v3merged['CRE_class']
    .apply(lambda x: ', '.join(map(str, x)) if isinstance(x, list) else x)  # Convert list to string if needed
    .apply(lambda x: set(x.split(', ')) if isinstance(x, str) else set())  # Convert to set for uniqueness
    .apply(lambda x: next(iter(x)) if len(x) == 1 else 'more than one')    # Extract single value or mark as "more than one"
)

print('DF len after adding classes:', len(anno_v3merged))

# SCREEN v4 all
SCREENv4_file = '/home/labs/davidgo/Collaboration/GenomeAnnotation/Human/SCREEN/Nadav_temp/all_cell_types/GRCh38-cCREs.bed'
print('importing chondrocyte specfic SCREEN')
SCREENv4_df = pd.read_csv(SCREENv4_file,sep='\t',header = None)
SCREENv4_df['coord'] = SCREENv4_df.apply(lambda row: f"{row[0]}:{str(row[1])}-{str(row[2])}", axis=1)

lifted_coords = liftOver(SCREENv4_df,chr_col = 0,start_col= 1, end_col= 2,source_assembly= 'hg38',target_assembly='hg19')

SCREENv4_lift = SCREENv4_df.merge(lifted_coords, left_on='coord',right_on=3,how='right')
SCREENv4_lift.drop(columns = ['0_x','1_x','2_x',4,'coord','lift_name','3_y'],inplace=True)
SCREENv4_lift.rename(columns = {5:'class','0_y':'chr','1_y':'start','2_y':'end','3_x':'ID'},inplace=True)
SCREENv4_lift = SCREENv4_lift[['chr','start','end','ID','class']]


SCREENv4_lift = SCREENv4_lift.sort_values(['chr','start'], ascending = [True, True])
SCREENv4_lift_bed = pybedtools.BedTool.from_dataframe(SCREENv4_lift)

intersect_SCREENv4 = anno_v3_bed_sort.intersect(SCREENv4_lift_bed, wao=True)
intersect_SCREENv4 = intersect_SCREENv4.to_dataframe()
intersect_SCREENv4 = intersect_SCREENv4[['name', 'itemRgb']].rename(columns = {'itemRgb':'SCREEN4_class'})
intersect_SCREENv4 = intersect_SCREENv4.groupby('name', as_index=False).agg({'SCREEN4_class': list})

intersect_SCREENv4['SCREEN4_class_strict'] = (
    intersect_SCREENv4['SCREEN4_class']
    .apply(lambda x: ', '.join(map(str, x)) if isinstance(x, list) else x)  # Convert list to string if needed
    .apply(lambda x: set(x.split(', ')) if isinstance(x, str) else set())  # Convert to set for uniqueness
    .apply(lambda x: next(iter(x)) if len(x) == 1 else 'more than one')    # Extract single value or mark as "more than one"
)


anno_v3merged = pd.merge(anno_v3merged, intersect_SCREENv4, left_on='id', right_on = "name", how='left')
anno_v3merged.drop(columns='name',inplace=True)





#SCREEN v4 chondrocytes

# SCREEN v4 all
SCREENv4chond_file = '/home/labs/davidgo/Collaboration/GenomeAnnotation/Human/SCREEN/Nadav_temp/in_vitro_differentiated_chondrocytes/ENCFF807AUZ_ENCFF466YVQ_ENCFF317LGP_ENCFF044ORH.bed'
print('importing',cell_type,'specfic SCREEN')
SCREENv4chond_df = pd.read_csv(SCREENv4chond_file,sep='\t',header = None)
SCREENv4chond_df['coord'] = SCREENv4chond_df.apply(lambda row: f"{row[0]}:{str(row[1])}-{str(row[2])}", axis=1)

lifted_coords = liftOver(SCREENv4chond_df,chr_col = 0,start_col= 1, end_col= 2,source_assembly= 'hg38',target_assembly='hg19')

SCREENv4chond_lift = SCREENv4chond_df.merge(lifted_coords, left_on='coord',right_on=3,how='right')
SCREENv4chond_lift.rename(columns = {9:'class','0_y':'chr','1_y':'start','2_y':'end','3_x':'ID'},inplace=True)
SCREENv4chond_lift = SCREENv4chond_lift[['chr','start','end','ID','class']]



SCREENv4chond_lift = SCREENv4chond_lift.sort_values(['chr','start'], ascending = [True, True])
SCREENv4chond_lift_bed = pybedtools.BedTool.from_dataframe(SCREENv4chond_lift)

intersect_SCREENv4chond = anno_v3_bed_sort.intersect(SCREENv4chond_lift_bed, wao=True)
intersect_SCREENv4chond = intersect_SCREENv4chond.to_dataframe()
intersect_SCREENv4chond = intersect_SCREENv4chond[['name', 'itemRgb']].rename(columns = {'itemRgb':'chond_SCREEN4_class'})
intersect_SCREENv4chond = intersect_SCREENv4chond.groupby('name', as_index=False).agg({'chond_SCREEN4_class': list})

intersect_SCREENv4chond['chond_SCREEN4_class_strict'] = (
    intersect_SCREENv4chond['chond_SCREEN4_class']
    .apply(lambda x: ', '.join(map(str, x)) if isinstance(x, list) else x)      # Convert list to string if needed
    .apply(lambda x: set(x.split(', ')) if isinstance(x, str) else set())       # Convert to set for uniqueness
    .apply(lambda x: x - {'Low-DNase'} if len(x) > 1 else x)                    # Remove 'Low-DNase' if more than one value exists    
    .apply(lambda x: next(iter(x)) if len(x) == 1 else 'more than one')         # Extract single value or mark as "more than one"
)


anno_v3merged = pd.merge(anno_v3merged, intersect_SCREENv4chond, left_on='id', right_on = "name", how='left')
anno_v3merged.drop(columns='name',inplace=True)



anno_v3merged.to_csv('/home/labs/davidgo/Collaboration/humanMPRA/top_candidates/chondrocytes/humanMPRA_annotations_v4.csv', header=True, index = False)




# # HAqers

# haqer_df = pd.read_csv(f'/home/labs/davidgo/Collaboration/humanMPRA/files_from_dropbox/lifted_haqers.bed', sep='\t',header=None)
# haqer_df = haqer_df.iloc[:, :4]
# haqer_bed = pybedtools.BedTool.from_dataframe(haqer_df)
# haqer_bed_sort = haqer_bed.sort()

# closest_haqer_bed = anno_v3_bed_sort.closest(haqer_bed_sort,D='ref', t="first")
# closest_haqer_df = closest_haqer_bed.to_dataframe()

# print(closest_haqer_df.head().to_string())

# closest_haqer_df = closest_haqer_df[['itemRgb', 'name']]

# merged = pd.merge(anno_v3, closest_haqer_df, left_on='id', right_on="name", how="left").drop("name", 1)
# merged = merged.rename(columns={'itemRgb': 'HAQER'})

# # ATAC_Seq
# atac = pd.read_excel(f'/home/labs/davidgo/Collaboration/humanMPRA/files_from_dropbox/Swain-Lenz_ATACseq_adipose.xlsx', sheet_name = "Derived", skiprows=1, header=0, usecols=["chr","start","end","padj"])
# atac_filter = atac[atac["padj"] < 0.05]
# atac_filter['start'] = atac_filter['start']-1
# atac_filter_bed = pybedtools.BedTool.from_dataframe(atac_filter)
# atac_filter_bed_sort = atac_filter_bed.sort()

# closest_atac_bed = anno_v3_bed_sort.closest(atac_filter_bed_sort,D='ref', t="first")
# closest_atac_df = closest_atac_bed.to_dataframe()

# print(closest_atac_df.head().to_string())

# closest_atac_df = closest_atac_df[['itemRgb', 'name']]

# merged = pd.merge(merged, closest_atac_df, left_on='id', right_on="name", how="left").drop("name", 1)
# merged = merged.rename(columns={'itemRgb': 'ATAC human/chimp gain'})

# # CNCCs - 5000 most active enhancers
# cnccs_active_enhancers = pd.read_excel(f'/home/labs/davidgo/Collaboration/humanMPRA/files_from_dropbox/Enhancers.xls', sheet_name = "top 5000 active enhancers", header=0, usecols=["chr (hg19)","start (hg19)","end (hg19)"])
# cnccs_active_enhancers['start (hg19)'] = cnccs_active_enhancers['start (hg19)']-1
# cnccs_active_enhancers_bed = pybedtools.BedTool.from_dataframe(cnccs_active_enhancers)
# cnccs_active_enhancers_bed_sort = cnccs_active_enhancers_bed.sort()

# closest_cnccs_active_enhancers_bed = anno_v3_bed_sort.closest(cnccs_active_enhancers_bed_sort, D='ref', t="first")
# closest_cnccs_active_enhancers_df = closest_cnccs_active_enhancers_bed.to_dataframe()

# print(closest_cnccs_active_enhancers_df.head().to_string())

# closest_cnccs_active_enhancers_df = closest_cnccs_active_enhancers_df.iloc[:,[3, -1]]
# closest_cnccs_active_enhancers_df = closest_cnccs_active_enhancers_df.rename(columns={closest_cnccs_active_enhancers_df.columns[-1]: 'CNCCs_5000_active_enhancers'})
# merged = pd.merge(merged, closest_cnccs_active_enhancers_df, left_on='id', right_on="name", how="left").drop("name", 1)

# # CNCCs - human biased enhancers
# cnccs_human = pd.read_excel(f'/home/labs/davidgo/Collaboration/humanMPRA/files_from_dropbox/Enhancers.xls', sheet_name = "top 1000 human-biased candidate", header=0, usecols=["chr (hg19)","start (hg19)","end (hg19)"])
# cnccs_human['start (hg19)'] = cnccs_human['start (hg19)']-1
# cnccs_human_bed = pybedtools.BedTool.from_dataframe(cnccs_human)
# cnccs_human_bed_sort = cnccs_human_bed.sort()

# closest_cnccs_human_bed = anno_v3_bed_sort.closest(cnccs_human_bed_sort, D='ref', t="first")
# closest_cnccs_human_df = closest_cnccs_human_bed.to_dataframe()

# print(closest_cnccs_human_df.head().to_string())

# closest_cnccs_human_df = closest_cnccs_human_df.iloc[:,[3, -1]]
# print(closest_cnccs_human_df.head().to_string())
# closest_cnccs_human_df = closest_cnccs_human_df.rename(columns={closest_cnccs_human_df.columns[-1]: 'CNCCs_human_biased_enhancers'})
# merged = pd.merge(merged, closest_cnccs_human_df, left_on='id', right_on="name", how="left").drop("name", 1)

# # CNCCs - chimp biased enhancers
# cnccs_chimp = pd.read_excel(f'/home/labs/davidgo/Collaboration/humanMPRA/files_from_dropbox/Enhancers.xls', sheet_name = "top 1000 chimp-biased candidate", header=0, usecols=["chr (hg19)","start (hg19)","end (hg19)"])
# cnccs_chimp['start (hg19)'] = cnccs_chimp['start (hg19)']-1
# cnccs_chimp_bed = pybedtools.BedTool.from_dataframe(cnccs_chimp)
# cnccs_chimp_bed_sort = cnccs_chimp_bed.sort()

# closest_cnccs_chimp_bed = anno_v3_bed_sort.closest(cnccs_chimp_bed_sort, D='ref', t="first")
# closest_cnccs_chimp_df = closest_cnccs_chimp_bed.to_dataframe()

# print(closest_cnccs_chimp_df.head().to_string())

# closest_cnccs_chimp_df = closest_cnccs_chimp_df.iloc[:,[3, -1]]
# closest_cnccs_chimp_df = closest_cnccs_chimp_df.rename(columns={closest_cnccs_chimp_df.columns[-1]: 'CNCCs_chimp_biased_enhancers'})
# merged = pd.merge(merged, closest_cnccs_chimp_df, left_on='id', right_on="name", how="left").drop("name", 1)

# # human specifc TFBS
# human_tfbs = pd.read_excel(f'/home/labs/davidgo/Collaboration/humanMPRA/files_from_dropbox/Glinsky2015_HumanSpecTFBS.xlsx', sheet_name = "Sheet1", header=0, usecols=["chr","start","end"])
# human_tfbs['start'] = pd.to_numeric(human_tfbs['start'], errors='raise')
# human_tfbs['end'] = pd.to_numeric(human_tfbs['end'], errors='raise')
# print(human_tfbs.shape[0])
# human_tfbs = human_tfbs[np.isfinite(human_tfbs['start'])].copy()
# print(human_tfbs.shape[0])
# human_tfbs['start'] = human_tfbs['start'].astype(int)
# human_tfbs['end'] = human_tfbs['end'].astype(int)
# human_tfbs['start'] = human_tfbs['start']-1
# human_tfbs_bed = pybedtools.BedTool.from_dataframe(human_tfbs)
# print(human_tfbs_bed.head())
# human_tfbs_bed_sort = human_tfbs_bed.sort()

# closest_human_tfbs_bed = anno_v3_bed_sort.closest(human_tfbs_bed_sort, D='ref', t="first")
# closest_human_tfbs_df = closest_human_tfbs_bed.to_dataframe()

# print(closest_human_tfbs_df.head().to_string())

# closest_human_tfbs_df = closest_human_tfbs_df.iloc[:,[3, -1]]
# print(closest_human_tfbs_df.head().to_string())
# closest_human_tfbs_df = closest_human_tfbs_df.rename(columns={closest_human_tfbs_df.columns[-1]: 'humanSpecTFBS_Glinsky'})
# merged = pd.merge(merged, closest_human_tfbs_df, left_on='id', right_on="name", how="left").drop("name", 1)

# # save df
# merged.to_csv(f'./top_candidates/chondrocytes/humanMPRA_annotations_v3.csv', header=True, index = False)