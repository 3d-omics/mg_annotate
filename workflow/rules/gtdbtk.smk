rule gtdbtk__classify_wf__:
    """Run GTDB-Tk over the dereplicated genomes."""
    input:
        mags=MAGS,
        database=features["databases"]["gtdbtk"],
    output:
        work_dir=directory(GTDBTK),
    log:
        RESULTS / "gtdbtk.log",
    conda:
        "../environments/gtdbtk.yml"
    shell:
        """
        export GTDBTK_DATA_PATH="{input.database}"

        gtdbtk classify_wf \
            --genome_dir {input.mags} \
            --extension fa \
            --out_dir {output.work_dir} \
            --cpus {threads} \
            --skip_ani_screen \
        2> {log} 1>&2
        """


rule gtdbtk__join_bac_and_ar__:
    input:
        work_dir=GTDBTK,
    output:
        summary=RESULTS / "gtdbtk.summary.tsv",
        bac_tree=RESULTS / "gtdbtk.backbone.bac120.classify.tree",
        ar_tree=touch(RESULTS / "gtdbtk.ar53.classify.tree"),
    log:
        GTDBTK / "gtdbtk.join.log",
    conda:
        "../environments/csvkit.yml"
    shell:
        """
        if [[ -f {input.work_dir}/gtdbtk.ar122.summary.tsv ]] ; then

            csvstack \
                --tabs \
                {input.work_dir}/gtdbtk.bac120.summary.tsv \
                {input.work_dir}/gtdbtk.ar53.summary.tsv \
            | csvformat \
                --out-tabs \
            > {output.summary} \
            2> {log}

            cp \
                --verbose \
                {input.work_dir}/classify/gtdbtk.ar53.classify.tree \
                {output.ar_tree} \
            2>> {log}

        else

            cp \
                --verbose \
                {input.work_dir}/gtdbtk.bac120.summary.tsv \
                {output.summary} \

        fi

        cp \
            --verbose \
            {input.work_dir}/classify/gtdbtk.backbone.bac120.classify.tree \
            {output.bac_tree} \
        2>> {log} 1>&2
        """


rule gtdbtk:
    input:
        rules.gtdbtk__join_bac_and_ar__.output,
