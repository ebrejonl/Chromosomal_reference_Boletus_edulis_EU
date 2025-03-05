
rule AT_content:
    input:"aligned/inter.hic"
    output:"Supplementary_Figure2.pdf"
    container: c_R 
    threads: 12
    shell:
     """
     Rscript Supp_figure.R
     """
