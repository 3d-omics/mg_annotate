rule checkm2:
    """Run CheckM2 over the dereplicated mags"""
    input:
        mags=MAGS,
        db=features["databases"]["checkm2"],
    output:
        RESULTS / "checkm2.quality_report.tsv",
    log:
        RESULTS / "checkm2.quality_report.log",
    conda:
        "checkm2.yml"
    params:
        out_dir=CHECKM,
    shell:
        """
        rm \
            --force \
            --recursive \
            --verbose \
            {params.out_dir} \
        2> {log} 1>&2

        checkm2 predict \
            --database_path {input.db}/uniref100.KO.1.dmnd \
            --extension .fa \
            --input {input.mags} \
            --output-directory {params.out_dir} \
            --remove_intermediates \
            --threads {threads} \
        2>> {log} 1>&2

        mv \
            --verbose \
            {params.out_dir}/quality_report.tsv \
            {output} \
        2>> {log} 1>&2

        rm \
            --force \
            --recursive \
            --verbose \
            {params.out_dir} \
        2>> {log} 1>&2
        """
