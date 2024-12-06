rule Telomeres:
    input:
        fasta = "../Haplotype2_renamed_reordered_Chr_only.fasta"
    output:
        "Data/Annotation/TELO_telomeric_repeat_windows.tsv"
    container: c_synteny 
    shell:
        """
        tidk explore --minimum 5 --distance 0.1 --maximum 60 {input.fasta} 
        tidk search --string CCTAA --output {output}
        tidk plot -t TELO_telomeric_repeat_windows.tsv -o Telo # for plotting
        """

