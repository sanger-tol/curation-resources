#!/bin/python

#! pip install python-circos
import sys
import pycircos
import collections
import matplotlib.pyplot as plt
import random
from datetime import datetime
import pandas as pd
import seaborn as sns
from Bio.SeqIO.FastaIO import SimpleFastaParser
import argparse
import os
import shutil

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--ref", required = True)
parser.add_argument("-q", "--query", required = True)
parser.add_argument("-l", "--lineage", required = True)
parser.add_argument("-ri", "--ref_id")
parser.add_argument("-qi", "--query_id")
parser.add_argument("-p", "--path")
parser.add_argument("-k", "--keep", action='store_true')
args = parser.parse_args()


if args.ref_id:
    sp1 = f'{args.ref_id}'
else:
    sp1 = 'Rf'
if args.query_id:
    sp2 = f'{args.query_id}'
else:
    sp2 = 'Qu'
lineage = f'{args.lineage}'

if args.path:
    base_path = f'{args.path}/'
else:
    base_path = os.getcwd()

names= f'{sp1}_{sp2}'

################## FORMAT BUSCO OUTPUT ##################

#Busco full table headers

#table_headers=["Busco id", "Status", "Sequence", "Gene Start", "Gene End", 'Strand', 'Score', 'Length', 'OrthoDB url', 'Description']
# Data types in busco full table
#my_types = {"Busco id": "string","Status": "string","Sequence": "string","Gene Start": "Int64","Gene End": "Int64",}


table_headers=["Busco id", "Status", "Sequence", "Gene Start", "Gene End"]
# Data types in busco full table
my_types = {"Busco id": "string","Status": "string","Sequence": "string","Gene Start": "Int64","Gene End": "Int64",}


def readfulltables(ID):
    ID_table = pd.read_csv(f'{base_path}/busco5_mini{ID}/run_{lineage}/full_table.tsv', sep="\t", skiprows=3, names=table_headers, dtype=my_types, usecols=["Busco id", "Status", "Sequence", "Gene Start", "Gene End"])
    ID_table['Sequence']=ID_table.Sequence.replace({'SUPER_':ID}, regex=True)
    ID_table = ID_table[(ID_table["Status"]=='Complete') & (ID_table.Sequence.str.startswith(ID))]
    return ID_table

def get_chroms_data(fasta, id):
    ff = open(fasta, 'r')
    ffdata = pd.DataFrame()
    chroms=[]
    end=[]
    start=[]
    for name, seq in SimpleFastaParser(ff):
        end.append(len(seq))
        chroms.append(name)
        start.append('1')
        
    ffdata['chr'] = chroms
    ffdata['start']= start
    ffdata['end'] = end

    ffdata = ffdata[ffdata['chr'].str.contains('SUPER_[a-zA-Z0-9]')]        
    ffdata['chr'] = ffdata.chr.replace({'SUPER_':id}, regex=True)
    filter = ffdata['chr'].str.contains('unloc')
    fffdata_filter = ffdata[~filter]
    return fffdata_filter

sp1_data = readfulltables(sp1)
sp2_data = readfulltables(sp2)


links = pd.merge(sp1_data, sp2_data, how ='inner', on =['Busco id']) 
links = links.rename(columns={'Sequence_x': 'chr1', 'Gene Start_x': 'start1', 'Gene End_x': 'end1','Sequence_y': 'chr2', 'Gene Start_y': 'start2', 'Gene End_y': 'end2' })
links_tidy = links.loc[:,['chr2', 'start2', 'end2','chr1','start1','end1']]
links_tidy['color'] = links_tidy['chr1']
links_tidy['end2'] = links_tidy['end2']+1000000
links_tidy['end1'] = links_tidy['end1']+1000000
links_tidy = links_tidy.rename(columns={'chr2':'chr1', 'start2':'start1', 'end2':'end1','chr1':'chr2','start1':'start2','end1':'end2'})

links_tidy.to_csv(f'{base_path}/{sp1}_{sp2}.links', sep = ',', index = False)

sp1_chroms = get_chroms_data(args.ref, sp1)
sp2_chroms = get_chroms_data(args.query, sp2)
chroms = pd.concat([sp2_chroms, sp1_chroms], axis=0)
chroms.to_csv(f'{base_path}/{sp1}_{sp2}_chrom.csv', sep = ',', index = False)

def bestSexMatch(df):
    RfZs = df[df.chr2.str.endswith(("Z", "X"))]
    zmdf = []
    if(RfZs.empty == True):
        pass
    else:
        for i in RfZs.chr2.unique():
            zmatchpcnt = (RfZs['chr2'].value_counts()[i]/len(RfZs))*100
            matcher = i, zmatchpcnt.round(2)
            zmdf.append(matcher)
        zmdf = pd.DataFrame(zmdf)
        zmdf = zmdf.rename(columns={0:'chr', 1:'pcnt'})
        maxZmatch = zmdf['pcnt'].idxmax()
        c = zmdf.iloc[maxZmatch].chr
        p = zmdf.iloc[maxZmatch].pcnt
        sexChrom = RfZs.chr1.unique()[0]
        zmdf.to_csv(f'Best_match_to_{sexChrom}_is_{c}_at_{p}%', sep = ',', index = False)

bestSexMatch(links_tidy)

################ PLOTTING ##############


now = datetime.now()
current_time = now.strftime('%Y-%m-%d_%H.%M.%S')
Garc    = pycircos.Garc
Gcircle = pycircos.Gcircle
l = links_tidy
c = chroms

chr1s = pd.Series(c['chr']).str.count(sp1)
chr1_count= chr1s.sum()
chr2s = pd.Series(c['chr']).str.count(sp2)
chr2_count= chr2s.sum()

if chr1_count > 12:
    palette_scheme = "Spectral"
else:
    palette_scheme = "Paired"

palette = sns.color_palette(palette_scheme, chr1_count)    # Testing colour palettes: Paired, Spectral, magma
chroms1_colours = []
for i in sns.color_palette(palette) :
    chroms1_colours.append(i)

chroms1_random = random.shuffle(chroms1_colours)

chroms2_colours = ['#CDC9C9'] * chr2_count
colour_list = chroms2_colours + chroms1_colours
#colour_list

# Add the colour list as a column of the croms df

c['colour'] = colour_list

colour_map = c.set_index('chr')['colour'].to_dict()

l['color'] = l['chr2'].apply(colour_map.get)

#Set chromosomes
circle = Gcircle(figsize=(8,8)) 

for i, row in c.iterrows():
    
    arc = Garc(
        arc_id=row['chr'], 
        size=row['end'], 
        interspace=2, 
        raxis_range=(935,985), 
        labelposition=80, 
        label_visible=True,
        facecolor=row['colour']
    )

    circle.add_garc(arc) 
    

circle.set_garcs()
circle.ax.set_title(f'{args.ref}_vs_{args.query}', pad = 40)

for i, row in l.iterrows():
    source = (row['chr1'], row['start1'], row['end1'], 920)
    destination = (row['chr2'], row['start2'], row['end2'], 920)
    circle.chord_plot(
        source, 
        destination, 
        facecolor=row['color']
    )


circle.figure.savefig(f"{base_path}/{args.ref}_vs_{args.query}.{current_time}.png", dpi=300, bbox_inches='tight' , pad_inches=0.4)


####### TIDY DIRECTORY  #########

def move_busco(ID):
    os.rename(f'{base_path}/busco5_mini{ID}/run_{lineage}/full_table.tsv', f'{base_path}/{ID}_BUSCO_full_table.tsv')

def tidy_busco_run(ID):
    shutil.rmtree(f'{base_path}/busco5_mini{ID}/') 

move_busco(sp1)
move_busco(sp2)

if args.keep:
    pass
else:
    tidy_busco_run(sp1)
    tidy_busco_run(sp2)
    shutil.rmtree(f'{base_path}/busco_downloads/') 
