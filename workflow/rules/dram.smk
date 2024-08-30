rule dram__setup__:
    """
    Set up the databases from DRAM, no matter what the config file says.
    """
    input:
        dram_db=features["databases"]["dram"],
    output:
        touch(RESULTS / "dram.setup"),
    log:
        RESULTS / "dram.setup.log",
    conda:
        "__environment__.yml"
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


rule dram__annotate__:
    """Annotate dereplicate genomes with DRAM"""
    input:
        fasta=MAGS / "{mag_id}.fa",
        dram_db=features["databases"]["dram"],
        setup=RESULTS / "dram.setup",
    output:
        work_dir=temp(directory(RESULTS / "dram.annotate" / "{mag_id}")),
    log:
        RESULTS / "dram.annotate" / "{mag_id}.log",
    conda:
        "__environment__.yml"
    params:
        min_contig_size=params["dram"]["annotate"]["min_contig_size"],
        work_dir=lambda w: RESULTS / "dram.annotate" / f"{w.mag_id}",
    shell:
        """
        rm -rf {params.work_dir}

        DRAM.py annotate \
            --input_fasta {input.fasta} \
            --output_dir {output.work_dir} \
            --threads 1 \
        2>> {log} 1>&2
        """


def collect_dram_annotate(wildcards):
    checkpoint_output = checkpoints.mags.get().output[0]
    checkpoint_output = checkpoints.mags.get().output[0]
    mag_ids = glob_wildcards(MAGS / "{mag_id}.fa").mag_id
    return [RESULTS / "dram.annotate" / mag_id for mag_id in mag_ids]


rule dram__annotate__aggregate_annotations__:
    """Aggregate DRAM annotations"""
    input:
        collect_dram_annotate,
    output:
        RESULTS / "dram.annotations.tsv.gz",
    log:
        RESULTS / "dram.annotate.aggregate.log",
    conda:
        "__environment__.yml"
    params:
        work_dir=RESULTS / "dram.annotate",
    shell:
        """
        ( csvstack \
            --tabs \
            {params.work_dir}/*/annotations.tsv \
        | csvformat \
            --out-tabs \
        | bgzip \
            --compress-level 9 \
        > {output} ) \
        2> {log}
        """


rule dram__annotate__aggregate_trnas__:
    """Aggregate DRAM tRNAs"""
    input:
        collect_dram_annotate,
    output:
        RESULTS / "dram.trnas.tsv",
    log:
        RESULTS / "dram.trnas.log",
    conda:
        "__environment__.yml"
    params:
        work_dir=RESULTS / "dram.annotate",
    shell:
        """
        ( csvstack \
            --tabs \
            {params.work_dir}/*/trnas.tsv \
        | csvformat \
            --out-tabs \
        > {output} ) \
        2> {log}
        """


rule dram_annotate__aggregate_rrnas__:
    """Aggregate DRAM rRNAs"""
    input:
        collect_dram_annotate,
    output:
        RESULTS / "dram.rrnas.tsv",
    log:
        RESULTS / "dram.rrnas.log",
    conda:
        "__environment__.yml"
    params:
        work_dir=RESULTS / "dram.annotate",
    shell:
        """
        ( csvstack \
            --tabs \
            {params.work_dir}/*/rrnas.tsv \
        | csvformat \
            --out-tabs \
        > {output} ) \
        2> {log}
        """


rule dram__annotate_archive__:
    """
    Create tarball once annotations are merged done
    """
    input:
        work_dirs=collect_dram_annotate,
        annotations=RESULTS / "dram.annotations.tsv.gz",
        trnas=RESULTS / "dram.trnas.tsv",
        rrnas=RESULTS / "dram.rrnas.tsv",
    output:
        tarball=RESULTS / "dram.annotate.tar.gz",
    log:
        RESULTS / "dram.archive.log",
    conda:
        "__environment__.yml"
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


rule dram__distill__:
    """Distill DRAM annotations."""
    input:
        annotations=RESULTS / "dram.annotations.tsv.gz",
        trnas=RESULTS / "dram.trnas.tsv",
        rrnas=RESULTS / "dram.rrnas.tsv",
        dram_db=features["databases"]["dram"],
        setup=RESULTS / "dram.setup",
    output:
        work_dir=temp(directory(RESULTS / "dram.distill")),
    log:
        RESULTS / "distill.log",
    conda:
        "__environment__.yml"
    shell:
        """
        DRAM.py distill \
            --input_file {input.annotations} \
            --output_dir {output.work_dir} \
            --rrna_path  {input.rrnas} \
            --trna_path  {input.trnas} \
        2>> {log} 1>&2
        """


rule dram__distill_archive__:
    input:
        work_dir=RESULTS / "dram.distill",
    output:
        genome=RESULTS / "dram.genome_stats.tsv",
        metabolism=RESULTS / "dram.metabolism_summary.xlsx",
        product_html=RESULTS / "dram.product.html",
        product_tsv=RESULTS / "dram.product.tsv",
    log:
        RESULTS / "dram.distill_archive.log",
    conda:
        "__environment__.yml"
    params:
        out_dir=RESULTS,
        dram_dir=DRAM,
    shell:
        """
        for file in genome_stats.tsv metabolism_summary.xlsx product.html product.tsv ; do

            cp \
                --verbose \
                {input.work_dir}/$file \
                {params.out_dir}/dram.$file \

        done 2> {log} 1>&2
        """


rule dram:
    """Run DRAM on dereplicated genomes."""
    input:
        rules.dram__annotate_archive__.output,
        rules.dram__distill_archive__.output,


localrules:
    dram__distill_archive__,
