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

        rsync \
            -Pravt \
            {input}/ \
            {output} \
        2> {log} 1>&2

        find \
            {output} \
            -name "*.fa.gz" \
            -exec gzip -d {{}} \\; \
        2>> {log} 1>&2
        """
