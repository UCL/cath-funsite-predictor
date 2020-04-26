#! /usr/bin/env python

import sys
from math import log

# run psiblast to get ascii_pssm

'''
psiblast -query test.faa -db ./blast_db/nr.db -out test_psiblast_faadb.txt -num_iterations=3 -evalue=0.001 -out_ascii_pssm ascii_pssm -inclusion_ethresh=0.001 -num_threads 4 -use_sw_tback
'''

pssmfile=sys.argv[1]

with open(pssmfile) as f:
    lines = f.readlines()
    
# calculate the shannon entropy based on wop
def calculate_shannon_entropy(list_wop):
    s=0
    for i in list_wop:
        j=int(i)
        #print(i,j)
        s = s + j
    #print(s)
  
    s2=0
    for i in list_wop:
        j=int(i)
        if j:
            a= log(j/s)/log(2)
            #print(log(2))
            b= (j*a)/s
            s2 = s2 + b
    
    s2 = s2 * -1
    #print(s,s2)
    score=float(s2)
    score=round(score,3)
    return score

# read each line for each residue and strip newline, match with the AA and pdb residue from rsa or scorecons map file, then generate the pssm and wop features
count=0
aa_list=[]

for line in lines:
    
    count=count+1
    list=line.rstrip('\n').split()
    
    if count == 3:
        aa_list=list[0:20]
        #print(aa_list)
        #print('\n')
        
    if count > 3 and list and list[0].isdigit() and count < len(lines):
        res_features=[]
        residue_count=list[0]
        aminoacid=list[1]
        list_pssm=list[2:22]
        list_wop=list[22:42]
        info_per_pos=list[42]
        gapless_match_to_pseudocounts=list[43]
            
        for (aa,pssm_i,wop_i) in zip(aa_list,list_pssm,list_wop):
            pssm_temp=aa+'_pssm'
            wop_temp=aa+'_wop'
            res_features.append(pssm_i)
            res_features.append(wop_i)
                
        entwop_score=calculate_shannon_entropy(list_wop)
            
        print(residue_count, aminoacid, entwop_score, sep='\t', end='\t')
        for entry in res_features:
            print(entry, sep='\t', end='\t')
        print(info_per_pos, gapless_match_to_pseudocounts, sep='\t')
