rule drep__dereplicate__:
    """Dereplicate all the bins using dRep."""
    input:
        genomes=MAGS,
    output:
        out_dir=temp(directory(RESULTS / "drep.{secondary_ani}.dir")),
    log:
        RESULTS / "drep.{secondary_ani}.log",
    conda:
        "__environment__.yml"
    params:
        secondary_ani=lambda w: w.secondary_ani,
        minimum_completeness=params["drep"]["minimum_completeness"],
        maximum_contamination=params["drep"]["maximum_contamination"],
    shell:
        """
        rm \
            --recursive \
            --force \
            --verbose \
            {output.out_dir}/data_tables \
            {output.out_dir}/data \
            {output.out_dir}/dereplicated_genomes \
            {output.out_dir}/figures \
            {output.out_dir}/log \
        2> {log} 1>&2

        dRep dereplicate \
            {output.out_dir} \
            --S_ani         {params.secondary_ani} \
            --completeness  {params.minimum_completeness} \
            --contamination {params.maximum_contamination} \
            --genomes       {input.genomes}/*.fa \
            --processors    {threads} \
        2>> {log} 1>&2
        """


rule drep__get_fasta__:
    input:
        out_dir=RESULTS / "drep.{secondary_ani}.dir",
    output:
        fasta=RESULTS / "drep.{secondary_ani}.fa.gz",
    log:
        RESULTS / "drep.{secondary_ani}.fa.log",
    conda:
        "__environment__.yml"
    shell:
        """
        ( cat \
            {input.out_dir}/dereplicated_genomes/*.fa \
        | bgzip \
            --compress-level 9 \
            --threads {threads} \
        > {output.fasta} \
        ) 2> {log}
        """


rule drep__tarball__:
    input:
        out_dir=RESULTS / "drep.{secondary_ani}.dir",
    output:
        tarball=RESULTS / "drep.{secondary_ani}.tar.gz",
    log:
        RESULTS / "drep.{secondary_ani}.tar.log",
    conda:
        "__environment__.yml"
    shell:
        """
        tar \
            --create \
            --directory {input.out_dir} \
            --file {output.tarball} \
            --remove-files \
            --use-compress-program="pigz --processes {threads}" \
            --verbose \
            ${{folder}} \
        2>> {log} 1>&2
        """


rule drep:
    input:
        [RESULTS / f"drep.{secondary_ani}.tar.gz" for secondary_ani in SECONDARY_ANIS],
        [RESULTS / f"drep.{secondary_ani}.fa.gz" for secondary_ani in SECONDARY_ANIS],
