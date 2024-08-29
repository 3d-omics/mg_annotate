rule dram__annotate__:
    """Annotate dereplicate genomes with DRAM in parallel"""
    input:
        mags=MAGS,
        dram_db=features["databases"]["dram"],
    output:
        out_dir=temp(directory(RESULTS / "dram.annotate")),
    log:
        RESULTS / "dram.annotate.log",
    conda:
        "__environment__.yml"
    params:
        min_contig_size=params["dram"]["annotate"]["min_contig_size"],
        out_dir=RESULTS,
        tmp_dir=RESULTS / "dram.annotate",
        parallel_retries=5,
    shell:
        """
        rm \
            --recursive \
            --force \
            --verbose \
            {params.tmp_dir} \
        2>> {log} 1>&2

        mkdir \
            --parents \
            {params.tmp_dir} \
        2>>{log} 1>&2

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

        ( find \
            results/mags \
            -name "*.fa" \
            -exec ls -al {{}} \\; \
        | sort \
            --key 5 \
            --numeric-sort \
            --reverse \
        | awk \
            '{{print $9}}' \
        | parallel \
            --jobs {threads} \
            DRAM.py annotate \
                --input_fasta {{}} \
                --output_dir {params.tmp_dir}/{{/.}} \
                --threads 1 \
        ) 2>> {log} 1>&2
        """


rule dram__annotate_archive__:
    input:
        out_dir=RESULTS / "dram.annotate",
    output:
        annotation=RESULTS / "dram.annotations.tsv.gz",
        trnas=RESULTS / "dram.trnas.tsv",
        rrnas=RESULTS / "dram.rrnas.tsv",
        tarball=RESULTS / "dram.annotate.tar.gz",
    log:
        RESULTS / "dram.archive.log",
    conda:
        "__environment__.yml"
    params:
        out_dir=RESULTS,
    shell:
        """
        for file in annotations trnas rrnas ; do

            ( csvstack \
                --tabs \
                {input.out_dir}/*/$file.tsv \
            | csvformat \
                --out-tabs \
            > {params.out_dir}/dram.$file.tsv \
            ) 2>> {log}

        done

        bgzip \
            --compress-level 9 \
            --threads {threads} \
            {params.out_dir}/dram.annotations.tsv \
        2>> {log} 1>&2

        tar \
            --create \
            --file {output.tarball} \
            --use-compress-program="pigz --processes {threads}" \
            --verbose \
            {input.out_dir} \
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
        rm \
            --recursive \
            --force \
            --verbose \
            {output.work_dir} \
        2> {log} 1>&2

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

            mv \
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
