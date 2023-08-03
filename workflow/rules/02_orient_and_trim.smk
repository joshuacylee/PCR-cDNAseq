# pychopper - reads orientation

rule pychopper:
    input:
        "data/{sample}.fastq.gz"
    output:
        "orienting/{sample}.fastq.gz"
    params:
        report = "orienting/{sample}.pdf",
        kit    = config["kit"],
        method = config["pychopper_method"],
        other  = config["pychopper_other_params"]
    conda:
        "../envs/cdna_mapping.yml"
    threads:
        config["pychopper_threads"]
    log:
        "orienting/{sample}.log"
    shell:
        """
        catfishq --max_n 0 {input} > {output}_in_tmp.fastq;
        cdna_classifier.py -t {threads} -r {params.report} \
            -k {params.kit} -m {params.method} {params.other} \
            {output}_in_tmp.fastq {output}_tmp 2> {log}
        rm {output}_in_tmp.fastq
        gzip -c {output}_tmp > {output}; rm {output}_tmp
        """

rule trim:
    input:
        rules.pychopper.output
    output:
        "trimming/{sample}.fastq.gz"
    params:
        short = "trimming/{sample}_too-short.fastq.gz",
        five_end_trim_seq     = lambda wc: "" + config["cutadapt_five_end_trim_seq"] + "",
        five_end_max_error    = config["cutadapt_five_end_max_error"],
        five_end_min_overlap  = config["cutadapt_five_end_min_overlap"],
        three_end_trim_seq    = lambda wc: "" + config["cutadapt_three_end_trim_seq"] + "",
        three_end_max_error   = config["cutadapt_three_end_max_error"],
        three_end_min_overlap = config["cutadapt_three_end_min_overlap"],
        min_length            = config["cutadapt_min_length"],
        preset                = lambda wc: "" + "--rename='{header} trimmed_seq_name={adapter_name}'"+ ""
    conda:
        "../envs/cdna_mapping.yml"
    threads:
        config["cutadapt_threads"]
    log:
        "trimming/{sample}.log"
    shell:
        """
        cutadapt -j {threads} --minimum-length {params.min_length} \
            -a "X{params.five_end_trim_seq};max_error_rate={params.five_end_max_error};min_overlap={params.five_end_min_overlap};optional...{params.three_end_trim_seq}X;max_error_rate={params.three_end_max_error};min_overlap={params.three_end_min_overlap};optional" \
            {params.preset} -o {output} --too-short-output {params.short} \
            {input} > {log}
        """
