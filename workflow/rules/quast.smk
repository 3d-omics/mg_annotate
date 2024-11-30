rule quast:
    """Run quast over one the dereplicated mags"""
    input:
        [RESULTS / f"drep.{secondary_ani}.fa.gz" for secondary_ani in SECONDARY_ANIS],
    output:
        directory(RESULTS / "quast"),
    log:
        RESULTS / "quast.log",
    conda:
        "../environments/quast.yml"
    threads: 4
    shell:
        """
        quast \
            --output-dir {output} \
            --threads {threads} \
            {input} \
        2> {log} 1>&2
        """


rule quast__all:
    input:
        rules.quast.output,
