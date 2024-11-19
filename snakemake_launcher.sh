## script to launch the pipeline on the cluster
#snakemake --jobs 50 \
#  -p \
#  --default-resources mem_mb=50000 threads=1 \
#  --use-singularity \
#  --singularity-args "--bind /prj/porcini/CONTAINERS_LOCAL/genotyping_wgs_data_v1.sif" \ ## local path of containers
#  --singularity-args "--bind /prj/porcini/CONTAINERS_LOCAL/pop_genomics_v0.1.sif" \
#  --singularity-args "--bind /prj/porcini/CONTAINERS_LOCAL/container_comp_genomics_v0.3.sif" \
#  --use-conda \
#  --rerun-triggers mtime \
#  --latency-wait 1000 \
#  --rerun-incomplete \
#  --keep-going \
#  --cluster '
#    sbatch \
#      --export ALL \
#      -n {threads} \
#      -e JOBS_err/{name}.{jobid}.err \
#      -o JOBS_err/{name}.{jobid}.out \
#      --mem={resources.mem_mb}' \
#      --jn Ref.{name}.{jobid}.sh #--dry-run


## script for laptop
snakemake --cores 17 --use-conda --rerun-incomplete --use-singularity \
    --singularity-args "--bind /home/etiennebrejon/Desktop/work/Dropbox/Reference_genome_21_02_2024/14_11/Chromosomal_reference_Boletus_edulis_EU"   #--dry-run #--rulegraph
