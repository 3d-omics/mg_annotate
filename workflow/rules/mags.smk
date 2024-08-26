rule mags:
    input:
        features["mag_catalogue_dir"],
    output:
        directory(MAGS),
    conda:
        "__environment__.yml"
    log:
        "results/mags.log",
    shell:
        """
        mkdir --parents {output} 2> {log}

        for file in {input}/*.fa.gz ; do
            gzip \
                --decompress \
                --stdout \
                $file \
            > {output}/$(basename $file .gz)
        done 2>> {log}
        """
