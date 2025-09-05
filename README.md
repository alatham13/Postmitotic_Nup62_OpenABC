\brief Code for modeling the NPC with Nup62 depletion experiments upon postmitotic exit

# Info
Code for modeling of protein condensates in OpenABC, https://github.com/ZhangGroup-MITChemistry/OpenABC.

_Author(s)_: Andrew Latham

_License_: MIT License

_Publications_:
- W Zhang, et al., Phase-separating nucleoporins are required for dilating nuclear membrane pores into selective transport channels after mitosis
- S Liu, et al., OpenABC enables flexible, simplified, and efficient GPU accelerated simulations of biomolecular condensates, PLoS Comput. Biol. 19(9): e1011442. 

_Requirements_:
- openmm (version 7.5.1 was used, requires < 7.6)
- MDanalysis (version 2.4.3 was used)
- MDtraj (version 1.9.7 was used)
- openmm-plumed (version 1.0 was used)
- Hardware: simulations were run on A40 GPUs and took up to 11 days. OpenMM allows for simulations in CPU, CUDA, or OpenCL platforms, although performance may vary on other platforms.

_Installation_:
Installation of prerequisites can take up to 1 hour. Note that the root directory for this repository (Postmitotic_Nup62_OpenABC) should be added to the path in order to use the python scripts in the openabc2 folder.

Set up a conda environment for your installation
```
conda create --name py39_test python=3.9
```

Install OpenMM
```
conda install -c conda-forge openmm=7.5.1
```

Install MDTraj
```
conda install -c conda-forge mdtraj=1.9.7
```

Install MDanalysis
```
conda install conda-forge::mdanalysis=2.4.3
```

Install Plumed
```
conda install openmm-plumed=1.0
```

_Contents_:
A model of Nup62 phase separation used to rationalize the force exerted by the Nup62 complex during assembly. This code is broken down into 2 folders:
1. openabc2 - a modified version of openabc, enabling scaled versions of the MOFF potential
2. condensate_model - the simulations used to analyze Nup phase separation
