rule dereplicate_dram_setup_db:
    """Set up the DRAM database."""
    input:
        features["dram_database"],
    output:
        touch(DREP_DRAM / "dram_db_setup.done"),
    log:
        DREP_DRAM / "dram_db_setup.log",
    conda:
        "dram.yml"
    shell:
        "python workflow/scripts/dram_setup.py {input} 2> {log}"


# rule dereplicate_dram_annotate:
#     """Annotate dereplicated genomes with DRAM."""
#     input:
#         drep_folder=DREP / "dereplicated_genomes",
#         mock_db=DREP_DRAM / "dram_db_setup.done",
#         gtdbtk_tree=GTDDBTK / "gtdbtk.bac120.summary.tree",
#     output:
#         outdir=DREP_DRAM / "annotate" / "dereplicated_genomes",
#         annotations=DREP_DRAM / "annotate" / "dereplicated_genomes" / "annotations.tsv",
#         trnas=touch(DREP_DRAM / "annotate" / "dereplicated_genomes" / "trnas.tsv"),
#         rrnas=touch(DREP_DRAM / "annotate" / "dereplicated_genomes" / "rrnas.tsv"),
#     log:
#         DREP_DRAM / "annotate" / "dereplicated_genomes.log",
#     conda:
#         "dram.yml"
#     params:
#         min_contig_size=1500,
#     threads: 24
#     resources:
#         mem_mb=double_ram(params["dereplicate"]["dram"]["memory_gb"]),
#         runtime=24 * 60,
#     retries: 5
#     shell:
#         """
#         rm -rfv {output.outdir} 2> {log} 1>&2

#         DRAM.py annotate \
#             --input_fasta "{input.drep_folder}/*.fa" \
#             --output_dir {output.outdir} \
#             --threads {threads} \
#             --min_contig_size {params.min_contig_size} \
#         2>> {log} 1>&2
#         """


rule dereplicate_dram_annotate:
    """Annotate dereplicate genomes with DRAM"""
    input:
        dereplicated_genomes=DREP / "dereplicated_genomes",
        mock_db=DREP_DRAM / "dram_db_setup.done",
        gtdbtk_summary=DREP_GTDBTK / "gtdbtk.summary.tsv",
    output:
        annotation=DREP_DRAM / "annotations.tsv",
        trnas=DREP_DRAM / "trnas.tsv",
        rrnas=DREP_DRAM / "rrnas.tsv",
    log:
        DREP_DRAM / "annotate.log",
    conda:
        "dram.yml"
    threads: 24
    params:
        min_contig_size=1500,
        tmp_dir=DREP_DRAM / "annotate",
    resources:
        mem_mb=double_ram(params["dereplicate"]["dram"]["memory_gb"]),
        runtime=48 * 60,
    retries: 5
    shell:
        """
        rm -rfv {params.tmp_dir} 2> {log} 1>&2
        mkdir --parents {params.tmp_dir} 2>>{log} 1>&2

        parallel \
            DRAM.py annotate \
                --input_fasta {{}} \
                --output_dir {params.tmp_dir}/{{/.}} \
                --threads 1 \
                --gtdb_taxonomy {input.gtdb_summary} \
        ::: {input.dereplicated_genomes}/*.fa

        for i in {annotations,trnas,rrnas} ; do
            ( csvstack \
                --tabs \
                {params.tmp_dir}*/annotations.tsv \
            | csvformat \
                --out-tabs \
            > $i.tsv \
            ) 2>> {log}
        done

        tar cvf - {params.tmp_dir} | pigz > {params.tmp_dir}.tar.gz
        """


rule dereplicate_dram_distill:
    """Distill DRAM annotations."""
    input:
        annotations=DREP_DRAM / "annotations.tsv",
        trnas=DREP_DRAM / "trnas.tsv",
        rrnas=DREP_DRAM / "rrnas.tsv",
        mock_db=DREP_DRAM / "dram_db_setup.done",
    output:
        genome=DREP_DRAM / "genome_stats.tsv",
        metabolism=DREP_DRAM / "metabolism_summary.xlsx",
        product_html=DREP_DRAM / "product.html",
        product_tsv=DREP_DRAM / "product.tsv",
    log:
        DREP_DRAM / "distill.log2",
    conda:
        "dram.yml"
    resources:
        mem_mb=double_ram(params["dereplicate"]["dram"]["memory_gb"]),
        runtime=24 * 60,
    params:
        outdir=DREP_DRAM,
    retries: 5
    shell:
        """
        DRAM.py distill \
            --input_file {input.annotations} \
            --rrna_path {input.rrnas} \
            --trna_path {input.trnas} \
            --output_dir {params.outdir} \
        2> {log} 1>&2
        """


rule dereplicate_dram:
    """Run DRAM on dereplicated genomes."""
    input:
        rules.dereplicate_dram_distill.output,


localrules:
    dereplicate_dram_setup_db,


# dram in parallel
# join anntations
# run distill
