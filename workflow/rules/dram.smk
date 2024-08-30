rule dram__annotate__:
    """Annotate dereplicate genomes with DRAM in parallel"""
    input:
        fasta=MAGS / "{mag_id}.fa",
        dram_db=features["databases"]["dram"],
    output:
        out_dir=temp(directory(RESULTS / "dram.annotate" / "{mag_id}")),
    log:
        RESULTS / "dram.annotate" / "{mag_id}.log",
    conda:
        "__environment__.yml"
    params:
        min_contig_size=params["dram"]["annotate"]["min_contig_size"],
        out_dir=RESULTS,
        tmp_dir=RESULTS / "dram.annotate",
        parallel_retries=5,
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

        DRAM.py annotate \
            --input_fasta {input.fasta} \
            --output_dir {output.out_dir} \
            --threads 1 \
        2>> {log} 1>&2
        """


def collect_dram_annotate_annotations(wildcards):
    return [
        RESULTS / "dram.annotate" / mag_id / "annotations.tsv"
        for mag_id in glob_wildcards(RESULTS / "dram.annotate" / "{mag_id}.fa").mag_id
    ]


def collect_dram_annotate_trnas(wildcards):
    return [
        RESULTS / "dram.annotate" / mag_id / "trnas.tsv"
        for mag_id in glob_wildcards(RESULTS / "dram.annotate" / "{mag_id}.fa").mag_id
    ]


def collect_dram_annotate_rrnas(wildcards):
    return [
        RESULTS / "dram.annotate" / mag_id / "rrnas.tsv"
        for mag_id in glob_wildcards(RESULTS / "dram.annotate" / "{mag_id}.fa").mag_id
    ]


rule dram__annotate_aggregate_annotations:
    input:
        collect_dram_annotate_annotations,
    output:
        RESULTS / "dram.annotations.tsv.gz",
    log:
        RESULTS / "dram.aggregate_annotations.log",
    conda:
        "__environment__.yml"
    shell:
        """
        ( csvstack \
            --tabs \
            {input} \
        | csvformat \
            --out-tabs \
        | bgzip \
            --compress-level 9 \
        > {output} ) \
        2> {log}
        """


rule dram__annotate_aggregate_trnas:
    input:
        collect_dram_annotate_trnas,
    output:
        RESULTS / "dram.trnas.tsv",
    log:
        RESULTS / "dram.aggregate_trnas.log",
    conda:
        "__environment__.yml"
    shell:
        """
        ( csvstack \
            --tabs \
            {input} \
        | csvformat \
            --out-tabs \
        > {output} ) \
        2> {log}
        """


rule dram__annotate_aggregate_rrnas:
    input:
        collect_dram_annotate_rrnas,
    output:
        RESULTS / "dram.rrnas.tsv",
    log:
        RESULTS / "dram.aggregate_rrnas.log",
    conda:
        "__environment__.yml"
    shell:
        """
        ( csvstack \
            --tabs \
            {input} \
        | csvformat \
            --out-tabs \
        > {output} ) \
        2> {log}
        """


rule dram__annotate_archive__:
    input:
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
        """


rule dram__distill__:
    """Distill DRAM annotations."""
    input:
        annotations=RESULTS / "dram.annotations.tsv.gz",
        trnas=RESULTS / "dram.trnas.tsv",
        rrnas=RESULTS / "dram.rrnas.tsv",
        dram_db=features["databases"]["dram"],
    output:
        work_dir=temp(directory(RESULTS / "dram.distill")),
    log:
        RESULTS / "distill.log",
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
        rules.dram__distill_archive__.output,


localrules:
    dram__distill_archive__,
