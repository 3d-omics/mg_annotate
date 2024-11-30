include: "dram_functions.smk"


rule dram__setup:
    """
    Set up the databases from DRAM, no matter what the config file says.
    """
    input:
        dram_db=features["databases"]["dram"],
    output:
        touch(RESULTS / "dram.setup.txt"),
    log:
        RESULTS / "dram.setup.log",
    conda:
        "../environments/dram.yml"
    shell:
        """
        DRAM-setup.py set_database_locations \
            --amg_database_loc          {input.dram_db}/amg_database.*.tsv \
            --dbcan_fam_activities_loc  {input.dram_db}/CAZyDB.*.fam-activities.txt \
            --dbcan_loc                 {input.dram_db}/dbCAN-HMMdb-V*.txt \
            --dbcan_subfam_ec_loc       {input.dram_db}/CAZyDB.*.fam.subfam.ec.txt \
            --description_db_loc        {input.dram_db}/description_db.sqlite \
            --etc_module_database_loc   {input.dram_db}/etc_mdoule_database.*.tsv \
            --function_heatmap_form_loc {input.dram_db}/function_heatmap_form.*.tsv \
            --genome_summary_form_loc   {input.dram_db}/genome_summary_form.*.tsv \
            --kofam_hmm_loc             {input.dram_db}/kofam_profiles.hmm \
            --kofam_ko_list_loc         {input.dram_db}/kofam_ko_list.tsv \
            --module_step_form_loc      {input.dram_db}/module_step_form.*.tsv \
            --peptidase_loc             {input.dram_db}/peptidases.*.mmsdb \
            --pfam_hmm_loc              {input.dram_db}/Pfam-A.hmm.dat.gz \
            --pfam_loc                  {input.dram_db}/pfam.mmspro \
            --viral_loc                 {input.dram_db}/refseq_viral.*.mmsdb \
            --vog_annotations_loc       {input.dram_db}/vog_annotations_latest.tsv.gz \
            --vogdb_loc                 {input.dram_db}/vog_latest_hmms.txt \
        2>> {log} 1>&2
        """


rule dram__annotate:
    """Annotate dereplicate genomes with DRAM"""
    input:
        fasta=RESULTS / "dram.mags" / "{mag_id}.fa",
        dram_db=features["databases"]["dram"],
        setup=RESULTS / "dram.setup.txt",
    output:
        work_dir=temp(directory(RESULTS / "dram.annotate" / "{mag_id}")),
    log:
        RESULTS / "dram.annotate" / "{mag_id}.log",
    conda:
        "../environments/dram.yml"
    params:
        min_contig_size=params["dram"]["annotate"]["min_contig_size"],
    resources:
        mem_mb=8 * 1024,
        runtime=24 * 60,
    shell:
        """
        rm \
            --recursive \
            --force \
            --verbose \
            {output.work_dir} \
        2> {log} 1>&2

        DRAM.py annotate \
            --input_fasta {input.fasta} \
            --output_dir {output.work_dir} \
            --threads 1 \
        2>> {log} 1>&2
        """


for file in ["annotations", "trnas", "rrnas"]:

    rule:
        name:
            f"dram__annotate__aggregate_{file}"
        input:
            collect_dram_annotate,
        output:
            RESULTS / f"dram.{file}.tsv.gz",
        log:
            RESULTS / f"dram.{file}.log",
        conda:
            "../environments/dram.yml"
        params:
            work_dir=RESULTS / "dram.annotate",
        shell:
            f"( csvtk concat {{params.work_dir}}/*/{file} "
            f"| sed -r 's/[[:alnum:]]+:bin_[0-9]+_([[:alnum:]]+:bin_[0-9]+@contig_[0-9]+)/\1/g' "
            f"| bgzip --compress-level 9 "
            f"> {{output}} "
            f") 2> {{log}}"


for file in ["genes.gff", "genes.fna", "genes.faa", "scaffolds.fna"]:

    rule:
        name:
            f"dram__annotate__concatenate_{file}"
        input:
            collect_dram_annotate,
        output:
            RESULTS / f"dram.{file}.gz",
        log:
            RESULTS / f"dram.{file}.log",
        conda:
            "../environments/dram.yml"
        params:
            work_dir=RESULTS / "dram.annotate",
        shell:
            f"( cat {{params.work_dir}}/*/{file} "
            f"| sed -r 's/[[:alnum:]]+:bin_[0-9]+_([[:alnum:]]+:bin_[0-9]+@contig_[0-9]+)/\1/g' "
            f"| bgzip --compress-level 9 "
            f"> {{output}}"
            f") 2> {{log}}"


rule dram__annotate__aggregate_genbank:
    """Aggregate all DRAM genbank files"""
    input:
        collect_dram_annotate,
    output:
        RESULTS / "dram.genbank.gbk.gz",
    log:
        RESULTS / "dram.genbank.log",
    conda:
        "../environments/dram.yml"
    params:
        work_dir=RESULTS / "dram.annotate",
    resources:
        runtime=6 * 60,
    shell:
        """
        ( cat \
            --verbose \
            {params.work_dir}/*/genbank/*.gbk \
        | sed \
            -r 's/[[:alnum:]]+:bin_[0-9]+_([[:alnum:]]+:bin_[0-9]+@contig_[0-9]+)/\1/g' \
        | bgzip \
            --compress-level 9 \
        > {output} \
        ) 2> {log}
        """


rule dram__annotate__archive:
    """
    Create tarball once annotations are merged done
    """
    input:
        annotations=RESULTS / "dram.annotations.tsv.gz",
        trnas=RESULTS / "dram.trnas.tsv.gz",
        rrnas=RESULTS / "dram.rrnas.tsv.gz",
        gtf=RESULTS / "dram.genes.gff.gz",
        fna=RESULTS / "dram.genes.fna.gz",
        faa=RESULTS / "dram.genes.faa.gz",
        scaffolds=RESULTS / "dram.scaffolds.fna.gz",
        genbank=RESULTS / "dram.genbank.gbk.gz",
    output:
        tarball=RESULTS / "dram.annotate.tar.gz",
    log:
        RESULTS / "dram.archive.log",
    conda:
        "../environments/dram.yml"
    threads: 24
    params:
        out_dir=RESULTS,
        work_dir=RESULTS / "dram.annotate",
    shell:
        """
        tar \
            --create \
            --file {output.tarball} \
            --use-compress-program="pigz --processes {threads}" \
            --verbose \
            {params.work_dir} \
        2>> {log} 1>&2

        rm -rfv {params.work_dir}
        """


rule dram__distill:
    """Distill DRAM annotations."""
    input:
        annotations=RESULTS / "dram.annotations.tsv.gz",
        trnas=RESULTS / "dram.trnas.tsv.gz",
        rrnas=RESULTS / "dram.rrnas.tsv.gz",
        dram_db=features["databases"]["dram"],
        setup=RESULTS / "dram.setup.txt",
    output:
        work_dir=temp(directory(RESULTS / "dram.distill")),
    log:
        RESULTS / "dram.distill.log",
    conda:
        "../environments/dram.yml"
    resources:
        mem_mb=16 * 1024,
        runtime=24 * 60,
    shell:
        """
        DRAM.py distill \
            --input_file {input.annotations} \
            --output_dir {output.work_dir} \
            --rrna_path  {input.rrnas} \
            --trna_path  {input.trnas} \
        2>> {log} 1>&2
        """


rule dram__distill__archive:
    input:
        work_dir=RESULTS / "dram.distill",
    output:
        genome=RESULTS / "dram.genome_stats.tsv",
        metabolism=RESULTS / "dram.metabolism_summary.xlsx",
        product_tsv=RESULTS / "dram.product.tsv",
    log:
        RESULTS / "dram.distill_archive.log",
    conda:
        "../environments/dram.yml"
    params:
        out_dir=RESULTS,
    threads: 24
    localrule: True
    shell:
        """
        for file in genome_stats.tsv metabolism_summary.xlsx product.tsv ; do

            cp \
                --verbose \
                {input.work_dir}/$file \
                {params.out_dir}/dram.$file \

        done 2> {log} 1>&2
        """


rule dram__all:
    """Run DRAM on dereplicated genomes."""
    input:
        rules.dram__annotate__archive.output,
        rules.dram__distill__archive.output,
