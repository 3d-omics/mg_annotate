rule drep__quality_report__:
    input:
        RESULTS / "checkm2.quality_report.tsv",
    output:
        temp(RESULTS / "drep.quality_report.tsv"),
    log:
        RESULTS / "drep.quality_report.log",
    conda:
        "__environment__.yml"
    shell:
        """
        echo \
            -E "genome,completeness,contamination" \
        > {output} \
        2> {log}

        ( tail -n+1 {input} \
        | cut -f 1-3 \
        | awk \
            '{{print $1 ".fa," $2 "," $3}}' \
        ) \
        >> {output} \
        2>> {log}
        """


rule drep__dereplicate__:
    """Dereplicate all the bins using dRep."""
    input:
        genomes=MAGS,
        quality_report=RESULTS / "drep.quality_report.tsv",
    output:
        work_dir=temp(directory(RESULTS / "drep.{secondary_ani}.dir")),
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
            {output.work_dir}/data_tables \
            {output.work_dir}/data \
            {output.work_dir}/dereplicated_genomes \
            {output.work_dir}/figures \
            {output.work_dir}/log \
        2> {log} 1>&2

        dRep dereplicate \
            {output.work_dir} \
            --S_ani         {params.secondary_ani} \
            --completeness  {params.minimum_completeness} \
            --contamination {params.maximum_contamination} \
            --genomes       {input.genomes}/*.fa \
            --processors    {threads} \
            --genomeInfo    {input.quality_report} \
        2>> {log} 1>&2
        """


rule drep__get_fasta__:
    input:
        work_dir=RESULTS / "drep.{secondary_ani}.dir",
    output:
        fasta=RESULTS / "drep.{secondary_ani}.fa.gz",
    log:
        RESULTS / "drep.{secondary_ani}.fa.log",
    conda:
        "__environment__.yml"
    shell:
        """
        ( cat \
            {input.work_dir}/dereplicated_genomes/*.fa \
        | bgzip \
            --compress-level 9 \
            --threads {threads} \
        > {output.fasta} \
        ) 2> {log}
        """


rule drep__tarball__:
    input:
        work_dir=RESULTS / "drep.{secondary_ani}.dir",
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
            --directory {input.work_dir} \
            --file {output.tarball} \
            --remove-files \
            --use-compress-program="pigz --processes {threads}" \
            --verbose \
            {input.work_dir} \
        2>> {log} 1>&2
        """


rule drep:
    input:
        [RESULTS / f"drep.{secondary_ani}.tar.gz" for secondary_ani in SECONDARY_ANIS],
        [RESULTS / f"drep.{secondary_ani}.fa.gz" for secondary_ani in SECONDARY_ANIS],
