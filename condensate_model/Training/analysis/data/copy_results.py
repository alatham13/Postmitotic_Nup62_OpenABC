# Function to calculate radius of gyration and distance between 2nd and 2nd-last residue (to compare to FRET)
# Written by Andrew Latham
# Note the inputs at top.

import sys
import os
import math

PDBs=['N49','N98','NSP','NUL','NUS']
weights=['0.9','0.95','1.0','1.05','1.1']
runs=['run1','run2','run3','run4','run5']

main_folder='/wynton/group/sali/aplatham/ip_assembly/Nup_phase_sep/train_MOFF/train_GPU4/'

for i in range(0,len(PDBs)):
    for j in range(0,len(weights)):
        for k in range(0,len(runs)):
            sim_folder=main_folder+PDBs[i]+'/eps_'+weights[j]+'/'+runs[k]+'/'
            oldfile1='FRET_v2.txt'
            newfile1=PDBs[i]+'_'+weights[j]+'_'+runs[k]+'_FRET.txt'
            oldfile2='Rg_v2.txt'
            newfile2 = PDBs[i] + '_' + weights[j] + '_' + runs[k] + '_Rg.txt'
            print('scp -i ~/.ssh/laptop_to_wynton aplatham@dt1.wynton.ucsf.edu:'+sim_folder+oldfile1+' '+newfile1)
            os.system('scp aplatham@dt1.wynton.ucsf.edu:'+sim_folder+oldfile1+' '+newfile1)
            os.system('scp aplatham@dt1.wynton.ucsf.edu:'+sim_folder+oldfile2+' '+newfile2)
