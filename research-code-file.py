#import libraries, then read in the annotation files

import pandas as pd
import numpy as np
from statistics import mean
import scipy.stats
import matplotlib.pyplot as plt

OPN1LW_df = pd.read_csv('OPN1LW_annotations.txt', sep='\t', header = 0, names = ['chr','cons','region','start','stop','junk1','junk2','junk3','id'])
OPN1MW_df = pd.read_csv('OPN1MW_annotations.txt', sep='\t', header = 0, names = ['chr','cons','region','start','stop','junk1','junk2','junk3','id'])
OPN1SW_df = pd.read_csv('OPN1SW_annotations.txt', sep='\t', header = 0, names = ['chr','cons','region','start','stop','junk1','junk2','junk3','id'])

chr7_var = pd.read_csv('chr7_variants.txt',sep='\t',header=0,names=['chr','start','stop','freqs','junk','junk2'])
chrX_var = pd.read_csv('chrX_variants.txt',sep='\t',header=0,names=['chr','start','stop','freqs','junk','junk2'])

#take a look at the data contained in the annotation files and use this to go back and add column names

print(OPN1LW_df.head())

#retrieve protein-coding regions (CDS) only

OPN1LW_CDS = OPN1LW_df[OPN1LW_df['region']=='CDS'].copy()
OPN1MW_CDS = OPN1MW_df[OPN1MW_df['region']=='CDS'].copy()
OPN1SW_CDS = OPN1SW_df[OPN1SW_df['region']=='CDS'].copy()

#retrieve the start and stop coordinates for each gene and load them to a list in order of position

OPN1LW_starts = OPN1LW_CDS['start'].tolist()
OPN1LW_stops = OPN1LW_CDS['stop'].tolist()
OPN1LW_coords = OPN1LW_starts + OPN1LW_stops
OPN1LW_coords.sort()

OPN1MW_starts = OPN1MW_CDS['start'].tolist()
OPN1MW_stops = OPN1MW_CDS['stop'].tolist()
OPN1MW_coords = OPN1MW_starts + OPN1MW_stops
OPN1MW_coords.sort()

OPN1SW_starts = OPN1SW_CDS['start'].tolist()
OPN1SW_stops = OPN1SW_CDS['stop'].tolist()
OPN1SW_coords = OPN1SW_starts + OPN1SW_stops
OPN1SW_coords.sort()

#read in the DNA sequence files

with open("chr7.fa", "r") as file1:
    f1 = file1.read()
chr7 = f1.replace("\n","")

with open("chrX.fa", "r") as file2:
    f2 = file2.read()
chrX = f2.replace("\n","")

#store the coordinates as start and stop values

SW_seq = ''
list_len = 0
while list_len < len(OPN1SW_coords):
	start = OPN1SW_coords[list_len]
	stop = OPN1SW_coords[list_len + 1]
	SW_seq += chr7[start:stop]
	list_len += 2
SW_dna_len = len(SW_seq)
print('OPN1SW is {} nucleotides long'.format(SW_dna_len))

LW_seq = ''
list_len = 0
while list_len < len(OPN1LW_coords):
	start = OPN1LW_coords[list_len]
	stop = OPN1LW_coords[list_len + 1]
	LW_seq += chrX[start:stop]
	list_len += 2
LW_dna_len = len(LW_seq)
print('OPN1LW is {} nucleotides long'.format(LW_dna_len))

MW_seq = ''
list_len = 0
while list_len < len(OPN1MW_coords):
	start = OPN1MW_coords[list_len]
	stop = OPN1MW_coords[list_len + 1]
	MW_seq += chrX[start:stop]
	list_len += 2
MW_dna_len = len(MW_seq)
print('OPN1MW is {} nucleotides long'.format(MW_dna_len))

#transcribe to DNA to RNA

LW_rna = ""
for i in LW_seq:
    if i == "A":
        LW_rna += "U"
    elif i == "C":
        LW_rna += "G"
    elif i == "G":
        LW_rna += "C"
    elif i == "T":
        LW_rna += "A"

MW_rna = ""
for i in MW_seq:
    if i == "A":
        MW_rna += "U"
    elif i == "C":
        MW_rna += "G"
    elif i == "G":
        MW_rna += "C"
    elif i == "T":
        MW_rna += "A"

SW_rna = ""
for i in SW_seq:
    if i == "A":
        SW_rna += "U"
    elif i == "C":
        SW_rna += "G"
    elif i == "G":
        SW_rna += "C"
    elif i == "T":
        SW_rna += "A"

#translate RNA to protein (dictionary was created by ChatGPT)

def translate(rna):
    table = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    protein = []
    for i in range(0, len(rna), 3):
        codon = rna [i:i+3]
        amino = table.get(codon, 'Unknown')
        protein.append(amino)
    return protein

LW_proteins = translate(LW_rna)
MW_proteins = translate(MW_rna)
SW_proteins = translate(SW_rna)

LW_protein_len = len(LW_proteins)
MW_protein_len = len(MW_proteins)
SW_protein_len = len(SW_proteins)

print('OPN1LW is {} proteins long'.format(LW_protein_len))
print('OPN1MW is {} proteins long'.format(MW_protein_len))
print('OPN1SW is {} proteins long'.format(SW_protein_len))

#extract the rows from the variant files that have the same start coordinate as the ones identified as pathogenic in dbSNP

LW_var_start = [154154734,154158844,154154602,154156440,154150812,154150897]
LW_var = chrX_var[chrX_var['start'].isin(LW_var_start)]

MW_var_start = [154191716,154187939,154195934,154190173,154182566,154195934]
MW_var = chrX_var[chrX_var['start'].isin(MW_var_start)]

SW_var_start = [128775556,128774545,128773786]
SW_var = chr7_var[chr7_var['start'].isin(SW_var_start)]

#there were not any pathogenic SNPs, so no further separation is needed

cds7var = pd.read_csv('cds7_overlapping_variants.txt', sep='\t',header=0,names=['chrom','start','stop','stuff','else','else2'])
cdsXvar = pd.read_csv('cdsX_overlapping_variants.txt', sep='\t',header=0,names=['chrom','start','stop','stuff','else','else2'])

#now lets separate these based on the specific genes

cdsLWvar = cdsXvar[(cdsXvar['start'] >= 154144243) & (cdsXvar['stop'] <= 154159032) & (cdsXvar['stuff'].str.contains('VT=SNP'))]
cdsMWvar = cdsXvar[(cdsXvar['start'] >= 154182678) & (cdsXvar['stop'] <= 154196040)& (cdsXvar['stuff'].str.contains('VT=SNP'))]
cdsSWvar = cds7var[(cds7var['start'] >= 128772540) & (cds7var['stop'] <= 128775781)& (cds7var['stuff'].str.contains('VT=SNP'))]

#parse for the population allelic frequencies (AF=)

goLW = cdsLWvar['stuff'].tolist()
goMW = cdsMWvar['stuff'].tolist()
goSW = cdsSWvar['stuff'].tolist()

newLW = [i.split(';') for i in goLW]
newMW = [i.split(';') for i in goMW]
newSW = [i.split(';') for i in goSW]

freqLW = []
for line in newLW:
    for item in line:
        if item.startswith('AF='):
            freqLW.append(item)

freqMW = []
for line in newMW:
    for item in line:
        if item.startswith('AF='):
            freqMW.append(item)

freqSW = []
for line in newSW:
    for item in line:
        if item.startswith('AF='):
            freqSW.append(item)

freqLWnum = []
for i in freqLW:
    numsplitLW = i.split('=')
    freqLWnum.append(numsplitLW[1])

freqMWnum = []
for i in freqMW:
    numsplitMW = i.split('=')
    freqMWnum.append(numsplitMW[1])

freqSWnum = []
for i in freqSW:
    numsplitSW = i.split('=')
    freqSWnum.append(numsplitSW[1])

freqLWnum = [float(i) for i in freqLWnum]
freqMWnum = [float(i) for i in freqMWnum]
freqSWnum = [float(i) for i in freqSWnum]

#average the frequencies for each gene

LWavg = mean(freqLWnum)
MWavg = mean(freqMWnum)
SWavg = mean(freqSWnum)

print('The average allelic frequency for OPN1LW is {}'.format(LWavg))
print('The average allelic frequency for OPN1MW is {}'.format(MWavg))
print('The average allelic frequency for OPN1SW is {}'.format(SWavg))

#visualize the average allelic frequencies
freq_data = [freqLWnum,freqMWnum,freqSWnum]

plt.boxplot(freq_data,labels=('OPN1LW','OPN1MW','OPN1SW'))
plt.title('Non-Pathogenic SNPs Average Allelic Frequencies')
plt.ylabel('Average Allelic Frequency')
plt.xlabel('Gene')
#plt.show()

#perform anova test to see if there is significant difference between average allelic frequencies of each gene (assume alpha = 0.05)

test = scipy.stats.f_oneway(freqLWnum,freqMWnum,freqSWnum)
print('Average allelic frequency t-test LW, MW, and SW: {}'.format(test))

#since p-val>0.05, the 3 means are not significantly different

#read in files for gene expression analysis

df1 = pd.read_csv("E-GEUV-1-query-results.fpkms.tsv",header=4,sep='\t')
df2 = pd.read_excel('just_sample_info.xlsx',header=0)

#separate populations into different dataframes

def col_move(df, keyword):
    df_new = pd.DataFrame()
    move = [col for col in df.columns if keyword in col]
    move.insert(0,df.columns[1])
    df_new = df[move].copy()

    return df_new

brit = col_move(df1,'British')
finl = col_move(df1,'Finland')
tusc = col_move(df1,'Tuscan')
utah = col_move(df1,'Utah')
yoru = col_move(df1,'Yoruba')

#separate the patients ids so they be joined on the population sample info file

def cleanup(df, substring):
    df.columns = df.columns.str.replace(substring, '')
    return df

brit_clean = cleanup(brit,'British, ')
finl_clean = cleanup(finl, 'Finland, ')
tusc_clean = cleanup(tusc, 'Tuscan, ')
utah_clean = cleanup(utah, 'Utah, ')
yoru_clean = cleanup(yoru, 'Yoruba, ')

brit_clean1 = brit_clean[(brit_clean['Gene Name']=='OPN1LW') | (brit_clean['Gene Name']=='OPN1MW') | (brit_clean['Gene Name']=='OPN1SW')]
finl_clean1 = finl_clean[(finl_clean['Gene Name']=='OPN1LW') | (finl_clean['Gene Name']=='OPN1MW') | (finl_clean['Gene Name']=='OPN1SW')]
tusc_clean1 = tusc_clean[(tusc_clean['Gene Name']=='OPN1LW') | (tusc_clean['Gene Name']=='OPN1MW') | (tusc_clean['Gene Name']=='OPN1SW')]
utah_clean1 = utah_clean[(utah_clean['Gene Name']=='OPN1LW') | (utah_clean['Gene Name']=='OPN1MW') | (utah_clean['Gene Name']=='OPN1SW')]
yoru_clean1 = yoru_clean[(yoru_clean['Gene Name']=='OPN1LW') | (yoru_clean['Gene Name']=='OPN1MW') | (yoru_clean['Gene Name']=='OPN1SW')]

#calculate the averages for each gene (swap in each population)

avg_B = brit_clean1.copy()
avg_B['Avg'] = avg_B.mean(axis=1)
checkB = avg_B.iloc[:,[0,95]]
sortedB = checkB.sort_values(by='Avg',ascending=False)

avg_T = tusc_clean1.copy()
avg_T['Avg'] = avg_T.mean(axis=1)
checkT = avg_T.iloc[:,[0,94]]
sortedT = checkT.sort_values(by='Avg',ascending=False)

avg_F = finl_clean1.copy()
avg_F['Avg'] = avg_F.mean(axis=1)
checkF = avg_F.iloc[:,[0,96]]
sortedF = checkF.sort_values(by='Avg',ascending=False)

avg_U = utah_clean1.copy()
avg_U['Avg'] = avg_U.mean(axis=1)
checkU = avg_U.iloc[:,[0,92]]
sortedU = checkU.sort_values(by='Avg',ascending=False)

#visualize gene expression for each population

SW_arr_B= np.array(brit_clean1.iloc[1,1:])
SW_arr_F= np.array(finl_clean1.iloc[1,1:])
SW_arr_T= np.array(tusc_clean1.iloc[1,1:])
SW_arr_U= np.array(utah_clean1.iloc[1,1:])
SW_arr_Y= np.array(yoru_clean1.iloc[1,1:])
all_data = [SW_arr_B,SW_arr_F,SW_arr_T,SW_arr_U,SW_arr_Y]

plt.boxplot(all_data,labels=('Brit','Finl','Tusc','Utah','Yoru'))
plt.title('Average OPN1SW Gene Expression by Population')
plt.ylabel('Average Gene Expression')
plt.xlabel('Population Group')
#plt.show()

#record mean values

avg_B = SW_arr_B.mean()
avg_T = SW_arr_T.mean()
avg_F = SW_arr_F.mean()
avg_U = SW_arr_U.mean()
avg_Y = SW_arr_Y.mean()

print('Average OPN1SW gene expression for Brit: {}'.format(avg_B))
print('Average OPN1SW gene expression for Tusc: {}'.format(avg_T))
print('Average OPN1SW gene expression for Finl: {}'.format(avg_F))
print('Average OPN1SW gene expression for Utah: {}'.format(avg_U))
print('Average OPN1SW gene expression for Yoru: {}'.format(avg_Y))

#compare avg gene expression among populations

exp_test = scipy.stats.f_oneway(SW_arr_B,SW_arr_F,SW_arr_T,SW_arr_U,SW_arr_Y)
print('One-way ANOVA test for SW gene expression of each population: {}'.format(exp_test))

#after visually reviewing the means, Yoruba seems like it may be the least similar population, so exclude and run again

exp_test2 = scipy.stats.f_oneway(SW_arr_B,SW_arr_F,SW_arr_T,SW_arr_U)
print('One-way ANOVA test for SW gene expression, Yoruba removed: {}'.format(exp_test2))

#still low p-value, exclude Utah and run again

exp_test3 = scipy.stats.f_oneway(SW_arr_B,SW_arr_F,SW_arr_T)
print('One-way ANOVA test for SW gene expression, Yoruba and Utah removed: {}'.format(exp_test3))

#now re-run analysis for male and female instead of population; cleanup sample info to have just gender and id

sample_info = df2.iloc[:,[0,4]]

#cleanup df1 to have just gene ids

exp_mat = df1.drop(df1.columns[0],axis=1)

exp_mat.columns = exp_mat.columns.str.replace('British, ', '')
exp_mat.columns = exp_mat.columns.str.replace('Tuscan, ','')
exp_mat.columns = exp_mat.columns.str.replace('Yoruba, ','')
exp_mat.columns = exp_mat.columns.str.replace('Utah, ','')
exp_mat.columns = exp_mat.columns.str.replace('Finland, ','')

notswapperoo = exp_mat.set_index(exp_mat.columns[0])
swapperoo = notswapperoo.T

#join the two dfs

joined_df = swapperoo.merge(sample_info, left_on=swapperoo.index,right_on=sample_info['Sample'],how='inner')

#put each gender in its own dataframe

female1 = joined_df[joined_df['Gender'] == 'female'].copy()
male1 = joined_df[joined_df['Gender'] == 'male'].copy()

female1 = female1.drop(['Sample','Gender'],axis=1)
male1 = male1.drop(['Sample','Gender'],axis=1)

female2 = female1.set_index(female1.columns[0])
male2 = male1.set_index(male1.columns[0])

female = female2.T
male = male2.T

#visualize gene expression for each sex

lw_fem = np.array(female[female.index=='OPN1LW'])
sw_fem = np.array(female[female.index=='OPN1SW'])
lw_fem = lw_fem.flatten()
sw_fem = sw_fem.flatten()

lw_mal = np.array(male[male.index=='OPN1LW'])
sw_mal = np.array(male[male.index=='OPN1SW'])
lw_mal = lw_mal.flatten()
sw_mal = sw_mal.flatten()

sw_comb = [sw_fem,sw_mal]

plt.boxplot(sw_comb,labels=('Female','Male'))
plt.title('Average OPN1SW Gene Expression by Sex')
plt.ylabel('Average Gene Expression')
plt.xlabel('Sex')
#plt.show()

#record average gene expression for each sex

avg1 = female
avg1['Avg'] = avg1.mean(axis=1)
check1 = avg1.iloc[:,[0,246]]
mygenes_f = check1[(check1.index=='OPN1LW') | (check1.index=='OPN1MW') | (check1.index=='OPN1SW')]
print('Average OPN1SW gene expression for females: {}'.format(mygenes_f))

avg2 = male
avg2['Avg'] = avg2.mean(axis=1)
check2 = avg2.iloc[:,[0,216]]
mygenes_m = check2[(check2.index=='OPN1LW') | (check2.index=='OPN1MW') | (check2.index=='OPN1SW')]
print('Average OPN1SW gene expression for females: {}'.format(mygenes_m))

#perform t test to compare gene expression for each sex

exp_test_sex = scipy.stats.ttest_ind(sw_mal,sw_fem,equal_var=False)
print('Two-sample t-test for OPN1SW gene expression for males and females: {}'.format(exp_test_sex))
