rule checkm2:
    """Run CheckM2 over the dereplicated mags"""
    input:
        mags=MAGS,
        db=features["databases"]["checkm2"],
    output:
        report=RESULTS / "checkm2.quality_report.tsv",
        tmp_dir=temp(directory(RESULTS / "checkm2.quality_report")),
    log:
        RESULTS / "checkm2.quality_report.log",
    conda:
        "../environments/checkm2.yml"
    threads: 24
    resources:
        mem_mb=16 * 1024,
        runtime=24 * 60,
    shell:
        """
        checkm2 predict \
            --threads {threads} \
            --input {input.mags} \
            --extension .fa \
            --output-directory {output.tmp_dir} \
            --database_path {input.db}/uniref100.KO.1.dmnd \
            --remove_intermediates \
        2>> {log} 1>&2

        cp \
            --verbose \
            {output.tmp_dir}/quality_report.tsv \
            {output.report} \
        2>> {log} 1>&2
        """


rule checkm2__all:
    input:
        rules.checkm2.output,
