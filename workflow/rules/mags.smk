checkpoint mags:
    input:
        features["mag_catalogue_dir"],
    output:
        directory(MAGS),
    conda:
        "base"
    log:
        RESULTS / "mags.log",
    localrule: True
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


rule mags__all:
    input:
        rules.mags.output,
