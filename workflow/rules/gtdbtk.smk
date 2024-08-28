rule gtdbtk__classify_wf__:
    """Run GTDB-Tk over the dereplicated genomes."""
    input:
        mags=MAGS,
        database=features["databases"]["gtdbtk"],
    output:
        identify=temp(directory(GTDBTK / "identify")),
        classify=temp(directory(GTDBTK / "classify")),
        align=temp(directory(GTDBTK / "align")),
        json=GTDBTK / "gtdbtk.json",
        bac120=GTDBTK / "bac120.summary.tsv",
        bac_tree=GTDBTK / "classify" / "gtdbtk.backbone.bac120.classify.tree",
        ar53=touch(GTDBTK / "ar53.summary.tsv"),
        ar_tree=touch(GTDBTK / "classify" / "gtdbtk.ar53.classify.tree"),
    log:
        RESULTS / "gtdbtk.log",
    conda:
        "__environment__.yml"
    params:
        out_dir=GTDBTK,
    shell:
        """
        export GTDBTK_DATA_PATH="{input.database}"

        gtdbtk classify_wf \
            --genome_dir {input.mags} \
            --extension fa \
            --out_dir {params.out_dir} \
            --cpus {threads} \
            --skip_ani_screen \
        2>> {log}
        """


rule gtdbtk__join_bac_and_ar__:
    input:
        bac_summary=GTDBTK / "bac120.summary.tsv",
        ar_summary=GTDBTK / "ar53.summary.tsv",
        bac_tree=GTDBTK / "classify" / "gtdbtk.backbone.bac120.classify.tree",
        ar_tree=GTDBTK / "classify" / "gtdbtk.ar53.classify.tree",
    output:
        summary=RESULTS / "gtdbtk.summary.tsv",
        bac_tree=RESULTS / "gtdbtk.backbone.bac120.classify.tree",
        ar_tree=RESULTS / "gtdbtk.ar53.classify.tree",
    log:
        GTDBTK / "gtdbtk.join.log",
    conda:
        "__environment__.yml"
    shell:
        """
        if [[ -s {input.ar_summary} ]] ; then

            ( csvstack \
                --tabs \
                {input.bac_summary} \
                {input.ar_summary} \
            | csvformat \
                --out-tabs \
            > {output.summary} \
            ) 2> {log}

            cp \
                --verbose \
                {input.ar_tree} \
                {output.ar_tree} \
            2>> {log} 1>&2

        else

            cp \
                --verbose \
                {input.bac_summary} \
                {output.summary} \
            2> {log} 1>&2

            touch {output.ar_tree} 2>> {log} 1>&2

        fi

        cp \
            --verbose \
            {input.bac_tree} \
            {output.bac_tree} \
        2>> {log}
    """


rule gtdbtk:
    input:
        rules.gtdbtk__join_bac_and_ar__.output,
