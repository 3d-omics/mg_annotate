rule gtdbtk:
    """Run GTDB-Tk over the dereplicated genomes."""
    input:
        mags=MAGS,
        database=features["databases"]["gtdbtk"],
    output:
        summary=RESULTS / "gtdbtk.summary.tsv",
        align=RESULTS / "gtdbtk.align.tar.gz",
        classify=RESULTS / "gtdbtk.classify.tar.gz",
        identify=RESULTS / "gtdbtk.identify.tar.gz",
        bac_tree=RESULTS / "gtdbtk.backbone.bac120.classify.tree",
        ar_tree=touch(RESULTS / "gtdbtk.ar53.tree"),
    log:
        RESULTS / "gtdbtk.log",
    conda:
        "__environment__.yml"
    params:
        out_dir=RESULTS,
        ar53=RESULTS / "gtdbtk.ar53.summary.tsv",
        bac120=RESULTS / "gtdbtk.bac120.summary.tsv",
        ar_tree=RESULTS / "classify" / "gtdbtk.ar53.tree",
        bac_tree=RESULTS / "classify" / "gtdbtk.backbone.bac120.classify.tree",
    resources:
        attempt=get_attempt,
    retries: 5
    shell:
        """
        rm \
            --recursive \
            --force \
            {params.out_dir}/gtdbtk.align \
            {params.out_dir}/gtdbtk.classify \
            {params.out_dir}/gtdbtk.identify \
            {params.out_dir}/gtdbtk.log \
            {params.out_dir}/gtdbtk.json \
            {params.out_dir}/gtdbtk.summary.tsv \
            {params.out_dir}/gtdbtk.warnings.log \
            {params.ar53} \
            {params.bac120} \
        2>> {log}.{resources.attempt} 1>&2

        export GTDBTK_DATA_PATH="{input.database}"

        gtdbtk classify_wf \
            --genome_dir {input.mags} \
            --extension fa \
            --out_dir {params.out_dir} \
            --cpus {threads} \
            --skip_ani_screen \
        2>> {log}.{resources.attempt} 1>&2

        if [[ -f {params.ar53} ]] ; then

            ( csvstack \
                --tabs \
                {params.bac120} \
                {params.ar53} \
            | csvformat \
                --out-tabs \
            > {output.summary} \
            ) 2>> {log}.{resources.attempt}

            cp \
                --verbose \
                {params.ar_tree} \
                {output.ar_tree} \
            2>> {log} 1>&2

        else

            cp \
                --verbose \
                {params.bac120} \
                {output.summary} \
            2>> {log}.{resources.attempt}

            touch {output.ar_tree} 2>> {log} 1>&2

        fi

        cp \
            --verbose \
            {params.bac_tree} \
            {output.bac_tree} \
        2>> {log}

        for folder in align classify identify ; do
            tar \
                --create \
                --directory {params.out_dir} \
                --file {params.out_dir}/gtdbtk.${{folder}}.tar.gz \
                --remove-files \
                --use-compress-program="pigz --processes {threads}" \
                --verbose \
                ${{folder}} \
            2>> {log}.{resources.attempt} 1>&2
        done

        mv \
            {log}.{resources.attempt} \
            {log}
        """
