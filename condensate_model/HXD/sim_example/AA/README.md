### This directory contains the inputs for running all atom simulations of an FG-Nup with HXD.

### To setup / run simulations:
#### Setup
gmx_mpi pdb2gmx -f N49_AF.pdb -o N49_processed.gro -ignh -ter
9 (CHARMM36 all-atom force field)
1 (TIP3P_CHARMM CHARMM-modified TIP3P water model)
0 (-NH3+ to N-terminous)
0 (-COO- to C-terminous)
#### Manually modify topol.top. Remove extra information / rename “PROA.itp”, move to 
gmx_mpi editconf -f N49_processed.gro -o N49_newbox.gro -c -d 2.0 -bt cubic
#### if using HXD (replace XXX with number of HXD molecules: 0% - 0, 1.5% - 69, 5.0% - 230):
gmx_mpi insert-molecules -f N49_newbox.gro -ci hexd.pdb -nmol  XXX -o N49_newbox.gro
#### Solvate (replace "SOL", 13 without HXD, 15 with HXD)
gmx_mpi solvate -cp N49_newbox.gro -cs spc216.gro -o N49_solvated.gro -p topol.top
gmx_mpi grompp -f ions.mdp -c N49_solvated.gro -p topol.top -o ions.tpr
gmx_mpi genion -s ions.tpr -p topol.top -o ions.gro -pname SOD -nname CLA -neutral -conc 0.162
13
15
#### Run sims
gmx_mpi grompp -f step4.0_minimization.mdp -c ions.gro -p topol.top -o min.tpr -r ions.gro
$GMXBIN/gmx_mpi mdrun -deffnm min -dlb yes -ntomp $NSLOTS -pin on
gmx_mpi grompp -f step4.1_equilibration.mdp -c min.gro -r ions.gro -p topol.top -o nvt.tpr
$GMXBIN/gmx_mpi mdrun -deffnm nvt -dlb yes -ntomp $NSLOTS -pin on -nb gpu
gmx_mpi grompp -f step4.2_equilibration.mdp -c nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
$GMXBIN/gmx_mpi mdrun -deffnm npt -dlb yes -ntomp $NSLOTS -pin on -nb gpu
gmx_mpi grompp -f step5_production.mdp -c npt.gro -t npt.cpt -p topol.top -o prod.tpr
$GMXBIN/gmx_mpi mdrun -deffnm prod -dlb yes -ntomp $NSLOTS -pin on -nb gpu

### To analyze simulations
- echo 0 | gmx_mpi trjconv -f prod.xtc -s prod.tpr -pbc mol -o prod_nojump1.xtc -e 2000000.000 -b 200000.000
- python ../../analysis/code/process_AA.py prod_nojump1.xtc