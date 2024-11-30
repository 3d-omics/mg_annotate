rule checkm2__predict__:
    """Run CheckM2 over the dereplicated mags"""
    input:
        mags=MAGS,
        db=features["databases"]["checkm2"],
    output:
        work_dir=temp(directory(CHECKM2)),
    log:
        RESULTS / "checkm2.predict.log",
    conda:
        "../environments/checkm2.yml"
    params:
        out_dir=CHECKM2,
    shell:
        """
        checkm2 predict \
            --database_path {input.db}/uniref100.KO.1.dmnd \
            --extension .fa \
            --input {input.mags} \
            --output-directory {output.work_dir} \
            --remove_intermediates \
            --threads {threads} \
        2>> {log} 1>&2
        """


rule checkm2__quality_report__:
    input:
        work_dir=CHECKM2,
    output:
        RESULTS / "checkm2.quality_report.tsv",
    log:
        RESULTS / "checkm2.quality_report.log",
    conda:
        "base"
    shell:
        """
        cp \
            --verbose \
            {input.work_dir}/quality_report.tsv \
            {output} \
        2> {log} 1>&2
        """


rule checkm2:
    input:
        rules.checkm2__quality_report__.output,


localrules:
    checkm2__quality_report__,
