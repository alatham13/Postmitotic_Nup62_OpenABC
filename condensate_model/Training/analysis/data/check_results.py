# Function to calculate radius of gyration and distance between 2nd and 2nd-last residue (to compare to FRET)
# Written by Andrew Latham
# Note the inputs at top.

import sys
import os
import math
import numpy as np

PDBs=['N49','N98','NSP','NUL','NUS']
#PDBs=['NUL']
weights=['0.9','0.95','1.0','1.05','1.1']
runs=['run1','run2','run3','run4','run5']
#runs=['run1','run2','run3','run4','run5','run6','run7']

main_folder='/wynton/group/sali/aplatham/ip_assembly/Nup_phase_sep/train_MOFF/train_GPU/'

count=0
count2=0
for i in range(0,len(PDBs)):
    for j in range(0,len(weights)):
        for k in range(0,len(runs)):
            newfile1=PDBs[i]+'_'+weights[j]+'_'+runs[k]+'_FRET.txt'
            newfile2 = PDBs[i] + '_' + weights[j] + '_' + runs[k] + '_Rg.txt'
            FRET=np.loadtxt(newfile1)
            Rg=np.loadtxt(newfile2)
            count=count+1
            #if np.max(FRET)>20:
            #    count2=count2+1
            #    print(np.max(FRET))
            #    print(newfile1)
            if np.max(Rg)>100:
                count2=count2+1
                print(np.max(Rg))
                print(newfile2)
print(count2)
print(count)
