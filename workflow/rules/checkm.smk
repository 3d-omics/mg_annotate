rule checkm__lineage_wf__:
    input:
        mags=MAGS,
    output:
        work_dir=temp(directory(RESULTS / "checkm.lineage_wf")),
    log:
        RESULTS / "checkm.lineage_wf.log",
    conda:
        "checkm.yml"
    shell:
        """
        checkm lineage_wf \
            --threads {threads} \
            --extension fa \
            {input.mags} \
            {output.work_dir} \
        2> {log} 1>&2
        """


rule checkm__qa__:
    input:
        work_dir=RESULTS / "checkm.lineage_wf",
    output:
        tsv=RESULTS / "checkm.qa.tsv",
    log:
        RESULTS / "checkm.qa.log",
    conda:
        "checkm.yml"
    shell:
        """
        checkm qa \
            --out_format 2 \
            --file {output.tsv} \
            --tab_table \
            --threads {threads} \
            {input.work_dir}/lineage.ms \
            {input.work_dir} \
        2> {log} 1>&2
        """


rule checkm__archive__:
    input:
        work_dir=RESULTS / "checkm.lineage_wf",
    output:
        archive=RESULTS / "checkm.lineage_wf.tar.gz",
    log:
        RESULTS / "checkm.archive.log",
    shell:
        """
        tar \
            --create \
            --file {output.archive} \
            --use-compress-program="pigz --processes {threads}" \
            --verbose \
            {input.work_dir} \
        2>> {log} 1>&2
        """


rule checkm:
    input:
        rules.checkm__qa__.output,
        rules.checkm__archive__.output,
