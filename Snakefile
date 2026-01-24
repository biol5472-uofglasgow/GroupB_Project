#Snakemake file for Gene Model Summariser

configfile: "biol5472_groupB_config.yaml"
from snakemake.io import directory

OUTDIR = config.get("outdir", "results/Gene_Model_Summariser")

rule all:
    input:
        f"{OUTDIR}/summary_report.txt",
        f"{OUTDIR}/run.json",
        f"{OUTDIR}/plots"

rule run_Gene_Model_Summariser:
    input:
        script=config["script"],
        gff=config["gff"],
        fasta=config["fasta"],
    output:
        summary_report=f"{OUTDIR}/summary_report.txt",
        run_json=f"{OUTDIR}/run.json",
        plots_dir=f"{OUTDIR}/plots",
    conda:
        "environment.yml"
    shell:
        r"""
        mkdir -p {OUTDIR}
        python {input.script} \
          --gff {input.gff} \
          --fasta {input.fasta} \
          --outdir {OUTDIR}
        """
