rule drep__:
    """Dereplicate all the bins using dRep."""
    input:
        genomes=MAGS,
    output:
        fasta=RESULTS / "drep.{secondary_ani}.fa.gz",
        tarball=RESULTS / "drep.{secondary_ani}.tar.gz",
        # dereplicated_genomes=directory(DREP / secondary_ani / "dereplicated_genomes"),
        # data=DREP / secondary_ani / "data.tar.gz",
        # data_tables=DREP / secondary_ani / "data_tables.tar.gz",
    log:
        RESULTS / "drep.{secondary_ani}.log",
    conda:
        "__environment__.yml"
    params:
        out_dir=lambda w: f"{DREP}/drep.{w.secondary_ani}",
        secondary_ani=lambda w: w.secondary_ani,
        minimum_completeness=params["drep"]["minimum_completeness"],
        maximum_contamination=params["drep"]["maximum_contamination"],
    resources:
        attempt=get_attempt,
    # retries: 5
    shell:
        """
        rm \
            --recursive \
            --force \
            {params.out_dir}/data_tables \
            {params.out_dir}/data \
            {params.out_dir}/dereplicated_genomes \
            {params.out_dir}/figures \
            {params.out_dir}/log \
        2>> {log}.{resources.attempt} 1>&2

        dRep dereplicate \
            {params.out_dir} \
            --S_ani {params.secondary_ani} \
            --completeness {params.minimum_completeness} \
            --contamination {params.maximum_contamination} \
            --genomes {input.genomes}/*.fa \
            --processors {threads} \
        2>> {log}.{resources.attempt} 1>&2

        ( cat \
            {params.out_dir}/dereplicated_genomes/*.fa \
        | bgzip \
            --compress-level 9 \
            --threads {threads} \
        > {output.fasta} \
        ) 2>> {log}.{resources.attempt} 1>&2

        cat \
            {params.out_dir}/dereplicated_genomes/*.fa \
        | bgzip \
            --compress-level 9 \
            --threads {threads} \
        > {output.fasta} \
        2>> {log}.{resources.attempt} 1>&2


        tar \
            --create \
            --directory {params.out_dir} \
            --file {output.tarball} \
            --remove-files \
            --use-compress-program="pigz --processes {threads}" \
            --verbose \
            ${{folder}} \
        2>> {log}.{resources.attempt} 1>&2

        cp --force {log}.{resources.attempt} {log}
        """


rule drep:
    input:
        [RESULTS / f"drep.{secondary_ani}.fa.gz" for secondary_ani in SECONDARY_ANIS],
