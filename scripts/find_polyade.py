#!/bin/env/python3

# script used in Alternative polyadenylation analysis (see iso-seq_analyses.sh)

import pandas as pd
import numpy as np
from scipy import signal
import os
filename = 'intersect.unique.reads.genes.stop.txt'
df1 = pd.read_csv(filename, sep='\t')
grouped = df1.groupby(['ID=FUN_000001;'])
my_list = list()
with open('peaks.pos.tsv', mode='wt', encoding='utf-8') as myfile1, open('peaks.heigth.tsv', mode='wt', encoding='utf-8') as myfile2, open('peaks.geneID.tsv', mode='wt', encoding='utf-8') as myfile3, open('peaks.scaffold.tsv', mode='wt', encoding='utf-8') as myfile4, open('peaks.strandedness.tsv', mode='wt', encoding='utf-8') as myfile5, open('peaks.pos.highest.tsv', mode='wt', encoding='utf-8') as myfile6, open('peaks.reads.highest.tsv', mode='wt', encoding='utf-8') as myfile7, open('peaks.reads.total.tsv', mode='wt', encoding='utf-8') as myfile8, open('peaks.pos.second.tsv', mode='wt', encoding='utf-8') as myfile9, open('peaks.reads.second.tsv', mode='wt', encoding='utf-8') as myfile10, open('peaks.pos.third.tsv', mode='wt', encoding='utf-8') as myfile11, open('peaks.reads.third.tsv', mode='wt', encoding='utf-8') as myfile12:
    for key, item in grouped:
        conv_arr= grouped.get_group(key).values
        arr1 = conv_arr[:,20]
        arr1 = arr1.ravel()
        arr1, bin_edges = np.histogram(arr1, bins=3000000, range=[0,2999999])
        my_list.append(signal.find_peaks(arr1, height=4, threshold=None, distance=18))
        arr2, dict1 = signal.find_peaks(arr1, height=4, threshold=None, distance=18)
        arr2.tofile(myfile1, sep='\t')
        myfile1.write('\n')
        for key in dict1.keys():
            dict1[key].tofile(myfile2, sep='\t')
            max_value = max(dict1[key], default=None)
            index_max_value = np.where(dict1[key] == max_value)[0]
            arr2[index_max_value].tofile(myfile6, sep='\t')
            myfile6.write('\n')
            if len(sorted(list(set(dict1[key].flatten().tolist())))) < 2:
                second = None
            else:
                second = sorted(list(set(dict1[key].flatten().tolist())))[-2]
            index_second_value = np.where(dict1[key] == second)[0]
            arr2[index_second_value].tofile(myfile9, sep='\t')
            myfile9.write('\n')
            if len(sorted(list(set(dict1[key].flatten().tolist())))) < 3:
                third = None
            else:
                third = sorted(list(set(dict1[key].flatten().tolist())))[-3]
            index_third_value = np.where(dict1[key] == third)[0]
            arr2[index_third_value].tofile(myfile11, sep='\t')
            myfile11.write('\n')
        myfile2.write('\n')
        myfile3.write(conv_arr[0,15] + '\n')
        myfile4.write(conv_arr[0,1] + '\n')
        myfile5.write(conv_arr[0,5] + '\n')
        conv_arr_max = conv_arr[conv_arr[:,20] == arr2[index_max_value]]
        conv_arr_max[:,0].tofile(myfile7, sep='\t')
        myfile7.write('\n')
        conv_arr_second = conv_arr[conv_arr[:,20] == arr2[index_second_value]]
        conv_arr_second[:,0].tofile(myfile10, sep='\t')
        myfile10.write('\n')
        conv_arr_third = conv_arr[conv_arr[:,20] == arr2[index_third_value]]
        conv_arr_third[:,0].tofile(myfile12, sep='\t')
        myfile12.write('\n')
        myfile8.write(str(conv_arr.shape[0]) + '\n')

with open('peaks.tsv', mode='wt', encoding='utf-8') as myfile:
    myfile.write('\n'.join('\t'.join(map(str, tup)) for tup in my_list))
    
os.system("awk '{print NF}' peaks.pos.tsv > peaks.pergene.txt")
os.system("paste peaks.geneID.tsv peaks.strandedness.tsv peaks.reads.total.tsv peaks.scaffold.tsv peaks.pos.highest.tsv peaks.reads.highest.tsv | sort -nr +2 | head -n 500 > peaks.highest.allinfo.500.tsv")
os.system("paste peaks.geneID.tsv peaks.strandedness.tsv peaks.reads.total.tsv peaks.scaffold.tsv peaks.pos.second.tsv peaks.reads.second.tsv | sort -nr +2 | awk 'NF >= 5' | head -n 500 > peaks.second.allinfo.500.tsv")
os.system("paste peaks.geneID.tsv peaks.strandedness.tsv peaks.reads.total.tsv peaks.scaffold.tsv peaks.pos.third.tsv peaks.reads.third.tsv | sort -nr +2 | awk 'NF >= 5' | head -n 500 > peaks.third.allinfo.500.tsv")
