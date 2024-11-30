def collect_dram_annotate(wildcards):
    checkpoint_output = checkpoints.mags.get().output[0]
    mag_ids = glob_wildcards(MAGS / "{mag_id}.fa").mag_id
    return [RESULTS / "dram.annotate" / mag_id for mag_id in mag_ids]
