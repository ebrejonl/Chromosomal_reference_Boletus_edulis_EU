# Boletus edulis chromosomal assembly <img src='https://github.com/user-attachments/assets/5d2d4735-d930-4a44-9abc-af66fcadc332' align="right" height="300" /></a>




### This is the repository that accompanies the publication: 
#### A haplotype-resolved chromosomal reference genome for the porcini mushroom Boletus edulis (will add doi of preprint soon)

This is a snakemake pipeline for all analysis, with software packaged in containers for reproducibility.

#### to run the pipeline, you need the following:
- [Snakemake](https://snakemake.readthedocs.io)
- a cluster system like [slurm](https://slurm.schedmd.com/documentation.html)
- [singularity/apptainer/Docker](https://apptainer.org/documentation/)
- 

The reference genome assembly and newly generated short reads can be found at PRJNA1187522 (available at publication)

#### The entire pipeline looks like this:
![workflow_dag](https://github.com/user-attachments/assets/11553806-7d36-4596-860a-3deb889e6ae5)




#### Individual steps can be found in the rules folder in their .smk files











