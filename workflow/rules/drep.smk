rule drep__quality_report:
    input:
        RESULTS / "checkm2.quality_report.tsv",
    output:
        temp(RESULTS / "drep.quality_report.tsv"),
    log:
        RESULTS / "drep.quality_report.log",
    conda:
        "../environments/drep.yml"
    localrule: True
    shell:
        """
        echo \
            -E "genome,completeness,contamination" \
        > {output} \
        2> {log}

        ( tail -n+2 {input} \
        | cut -f 1-3 \
        | awk \
            '{{print $1 ".fa," $2 "," $3}}' \
        ) \
        >> {output} \
        2>> {log}
        """


rule drep__dereplicate:
    """Dereplicate all the bins using dRep."""
    input:
        genomes=MAGS,
        quality_report=RESULTS / "drep.quality_report.tsv",
    output:
        work_dir=temp(directory(RESULTS / "drep.{secondary_ani}.dir")),
    log:
        RESULTS / "drep.{secondary_ani}.log",
    conda:
        "../environments/drep.yml"
    params:
        secondary_ani=lambda w: w.secondary_ani,
        minimum_completeness=params["drep"]["minimum_completeness"],
        maximum_contamination=params["drep"]["maximum_contamination"],
    threads: 24
    resources:
        mem_mb=double_ram(4 * 1024),
        runtime=6 * 60,
    shell:
        """
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


rule drep__get_fasta:
    input:
        work_dir=RESULTS / "drep.{secondary_ani}.dir",
    output:
        fasta=RESULTS / "drep.{secondary_ani}.fa.gz",
    log:
        RESULTS / "drep.{secondary_ani}.fa.log",
    conda:
        "../environments/drep.yml"
    threads: 24
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


rule drep__tarball:
    input:
        work_dir=RESULTS / "drep.{secondary_ani}.dir",
    output:
        tarball=RESULTS / "drep.{secondary_ani}.tar.gz",
    log:
        RESULTS / "drep.{secondary_ani}.tar.log",
    conda:
        "../environments/drep.yml"
    threads: 24
    shell:
        """
        tar \
            --create \
            --file {output.tarball} \
            --use-compress-program="pigz --processes {threads}" \
            --verbose \
            {input.work_dir} \
        2>> {log} 1>&2
        """


rule drep__all:
    input:
        [RESULTS / f"drep.{secondary_ani}.tar.gz" for secondary_ani in SECONDARY_ANIS],
        [RESULTS / f"drep.{secondary_ani}.fa.gz" for secondary_ani in SECONDARY_ANIS],
