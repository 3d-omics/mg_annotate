rule dram__annotate__:
    """Annotate dereplicate genomes with DRAM"""
    input:
        mags=MAGS,
        gtdbtk_summary=RESULTS/ "gtdbtk.summary.tsv",
        dram_db=features["databases"]["dram"],
    output:
        annotation=RESULTS / "dram.annotations.tsv.gz",
        trnas=RESULTS / "dram.trnas.tsv",
        rrnas=RESULTS / "dram.rrnas.tsv",
        tarball=RESULTS / "annotate.tar.gz",
    log:
        RESULTS / "dram.annotate.log",
    conda:
        "__environment__.yml"
    params:
        min_contig_size=params["dram"]["annotate"]["min_contig_size"],
        out_dir=RESULTS,
        tmp_dir=RESULTS / "annotate",
        parallel_retries=5,
    shell:
        """
        rm \
            --recursive \
            --force \
            --verbose {params.tmp_dir} \
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

        ( ls -1S {input.mags} \
        | parallel \
            --jobs $(( {threads} )) \
            DRAM.py annotate \
                --input_fasta {{}} \
                --output_dir {params.tmp_dir}/{{/.}} \
                --threads 1 \
                --gtdb_taxonomy {input.gtdbtk_summary} \
        ) 2>> {log} 1>&2

        for file in annotations trnas rrnas ; do
            
            ( csvstack \
                --tabs \
                {params.tmp_dir}/*/$file.tsv \
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
            --directory {params.out_dir} \
            --file {output.tarball} \
            --remove-files \
            --use-compress-program="pigz --processes {threads}" \
            --verbose \
            annotate \
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
        genome=RESULTS / "dram.genome_stats.tsv",
        metabolism=RESULTS / "dram.metabolism_summary.xlsx",
        product_html=RESULTS / "dram.product.html",
        product_tsv=RESULTS / "dram.product.tsv",
    log:
        DRAM / "distill.log",
    conda:
        "__environment__.yml"
    params:
        outdir_tmp=DRAM / "distill",
        outdir=DRAM,
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
            --rrna_path {input.rrnas} \
            --trna_path {input.trnas} \
            --output_dir {params.outdir_tmp} \
        2> {log} 1>&2

        mv \
            --verbose \
            {params.outdir_tmp}/genome_stats.tsv \
            {params.outdir}/dram.genome_stats.tsv \
        2>> {log} 1>&2

        mv \
            --verbose \
            {params.outdir_tmp}/metabolism_summary.xlsx \
            {params.outdir}/dram.metabolism_summary.xlsx \
        2>> {log} 1>&2
        
        mv \
            --verbose \
            {params.outdir_tmp}/product.html \
            {params.outdir}/dram.product.html \
        2>> {log} 1>&2
        
        mv \
            {params.outdir_tmp}/product.tsv \
            {params.outdir}/dram.product.tsv \
        2>> {log} 1>&2

        rmdir \
            --verbose \
            {params.outdir_tmp} \
        2>> {log} 1>&2
        """


rule dram:
    """Run DRAM on dereplicated genomes."""
    input:
        rules.dram__distill__.output,