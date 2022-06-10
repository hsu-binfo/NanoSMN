from pathlib import Path

configfile: "config.yml"

INPUTDIR = Path(config["inputdir"])
OUTDIR = Path(config["outdir"])
LOGDIR = Path(config["logdir"])
all_outputs = []

# SAMPLES = set(glob_wildcards(INPUTDIR/"samples"/config["input_fn_pattern"]).sample)
SAMPLES = ['SMN15_D15347_NP11001_BC15']
TYPES = ('SMN1', 'SMN2')
# print(TYPES)

# all_pre_filter_stat = expand(str(OUTDIR/"{sample}/{sample}.stat"), sample=SAMPLES)
# all_outputs.append(all_pre_filter_stat)

# all_NanoFilt = expand(str(OUTDIR/"{sample}/{sample}.filtered.fq.gz"), sample=SAMPLES)
# all_outputs.append(all_NanoFilt)

# all_post_filter_stat = expand(str(OUTDIR/"{sample}/{sample}.filtered.stat"), sample=SAMPLES)
# all_outputs.append(all_post_filter_stat)

# all_graph_alignment = expand(str(OUTDIR/"{sample}/{sample}.aln.gaf"), sample=SAMPLES)
# all_outputs.append(all_graph_alignment)

# all_cnv_plot = expand(str(OUTDIR/"{sample}/{sample}.cnv.png"), sample=SAMPLES)
# all_outputs.append(all_cnv_plot)

# all_read_list = expand(str(OUTDIR/"{sample}/SMN1_reads.lst"), sample=SAMPLES)
# all_outputs.append(all_read_list)

# all_split_read = expand(str(OUTDIR/"{sample}/{sample}.filtered.{type}.fq"), sample=SAMPLES, type=TYPES)
# all_outputs.append(all_split_read)

# all_alignment = expand(str(OUTDIR/"{sample}/{sample}.{type}.bam"), sample=SAMPLES, type=TYPES)
# all_outputs.append(all_alignment)

# all_variant_calling = expand(str(OUTDIR/"{sample}/{type}/variant_calls.final.vcf.gz"), sample=SAMPLES, type=TYPES)
# all_outputs.append(all_variant_calling)

# all_variant_annotation = expand(str(OUTDIR/"{sample}/{type}.variant_calls.final.vep"), sample=SAMPLES, type=TYPES)
# all_outputs.append(all_variant_annotation)

all_final_report = expand(str(OUTDIR/"{sample}/report.html"), sample=SAMPLES)
all_outputs.append(all_final_report)

rule all:
  input:
    all_outputs

#############################
# Pre-processing
#############################
rule pre_filter_stat:
  input:
    fastq=INPUTDIR/"samples"/config["input_fn_pattern"].format(sample="{sample}"),
  output:
    stat=OUTDIR/"{sample}/{sample}.stat",
  log:
    stdout=str(LOGDIR/"NanoFilt/{sample}.pre_filter_stat.stdout.log"),
    stderr=str(LOGDIR/"NanoFilt/{sample}.pre_filter_stat.stderr.log"),
  shadow:
    "shallow"
  conda:
    "envs/main.yml"
  threads: 10
  shell:
    """
    NanoStat \
      --fastq {input.fastq} \
      --name {output.stat}
      > {log.stdout} \
      2> {log.stderr}
    """

NanoFilt_config = config["NanoFilt"]
rule NanoFilt:
  input:
    fastq=INPUTDIR/"samples"/config["input_fn_pattern"].format(sample="{sample}"),
  output:
    filtered_fastq=OUTDIR/"{sample}/{sample}.filtered.fq" if config["keep_output"] else temp(OUTDIR/"NanoFilt/{sample}.filtered.fq"),
  log:
    # stdout=str(LOGDIR/"NanoFilt/{sample}.stdout.log"),
    stderr=str(LOGDIR/"NanoFilt/{sample}.stderr.log"),
  shadow:
    "shallow"
  conda:
    "envs/main.yml"
  threads: 10
  params:
    extra=NanoFilt_config["extra"],
  shell:
    """
    zcat {input.fastq} | \
    NanoFilt \
      -q 10 \
      {params.extra} \
      > {output.filtered_fastq} \
      2> {log.stderr}
    """

rule post_filter_stat:
  input:
    fastq=OUTDIR/"{sample}/{sample}.filtered.fq".format(sample="{sample}"),
  output:
    stat=OUTDIR/"{sample}/{sample}.filtered.stat",
  log:
    stdout=str(LOGDIR/"NanoFilt/{sample}.post_filter_stat.stdout.log"),
    stderr=str(LOGDIR/"NanoFilt/{sample}.post_filter_stat.stderr.log"),
  shadow:
    "shallow"
  conda:
    "envs/main.yml"
  threads: 10
  shell:
    """
    NanoStat \
      --fastq {input.fastq} \
      --name {output.stat}
      > {log.stdout} \
      2> {log.stderr}
    """

#############################
# Graph Alignment
#############################
graph_alignment_config = config["graph_alignment"]
rule graph_alignment:
  input:
    filtered_fastq=OUTDIR/"{sample}/{sample}.filtered.fq".format(sample="{sample}"),
    ref_gfa=INPUTDIR/"SMN_amplicon_align.gfa",
  output:
    gaf=OUTDIR/"{sample}/{sample}.aln.gaf" if config["keep_output"] else temp(OUTDIR/"{sample}/{sample}.aln.gaf"),
  log:
    stdout=str(LOGDIR/"GraphAligner/{sample}.stdout.log"),
    stderr=str(LOGDIR/"GraphAligner/{sample}.stderr.log"),
  shadow:
    "shallow"
  conda:
    "envs/main.yml"
  threads: 10
  params:
    extra=graph_alignment_config["extra"],
  shell:
    """
    GraphAligner \
      -t {threads} \
      -g {input.ref_gfa} \
      -f {input.filtered_fastq} \
      -a {output.gaf} \
      -x vg \
      {params.extra} \
      > {log.stdout} \
      2> {log.stderr}
    """

#############################
# CNV
#############################
rule cnv_plot:
  input:
    gaf=OUTDIR/"{sample}/{sample}.aln.gaf".format(sample="{sample}"),
  output:
    cnv_plot=OUTDIR/"{sample}/{sample}.cnv.png",
  log:
    stdout=str(LOGDIR/"cnv_plot/{sample}.stdout.log"),
    stderr=str(LOGDIR/"cnv_plot/{sample}.stderr.log"),
  shadow:
    "shallow"
  conda:
    "envs/main.yml"
  threads: 10
  params:
    sample="{sample}".format(sample="{sample}"),
    path=OUTDIR
  shell:
    """
    python utils/draw_cnv_plot.py \
      --path {params.path} \
      --sample {params.sample} \
      > {log.stdout} \
      2> {log.stderr}
    """

#############################
# Variant Calling
#############################
rule read_list:
  input:
    filtered_fastq=OUTDIR/"{sample}/{sample}.filtered.fq".format(sample="{sample}"),
    gaf=OUTDIR/"{sample}/{sample}.aln.gaf",
  output:
    SMN1_list=OUTDIR/"{sample}/SMN1_reads.lst",
    SMN2_list=OUTDIR/"{sample}/SMN2_reads.lst",
  log:
    stdout=str(LOGDIR/"variant_calling/{sample}.read_list.stdout.log"),
    stderr=str(LOGDIR/"variant_calling/{sample}.read_list.stderr.log"),
  shadow:
    "shallow"
  conda:
    "envs/main.yml"
  threads: 10
  params:
    sample="{sample}",
    path=OUTDIR
  shell:
    """
    python utils/generate_read_list_by_type.py \
      --path {params.path} \
      --sample {params.sample} \
      > {log.stdout} \
      2> {log.stderr}
    """

rule split_read:
  input:
    filtered_fastq=OUTDIR/"{sample}/{sample}.filtered.fq".format(sample="{sample}"),
    SMN_list=OUTDIR/"{sample}/{type}_reads.lst",
  output:
    SMN_fastq=OUTDIR/"{sample}/{sample}.filtered.{type}.fq" if config["keep_output"] else temp(OUTDIR/"{sample}/{sample}.filtered.{type}.fq"),
  log:
    # stdout=str(LOGDIR/"variant_calling/{sample}.split_read.stdout.log"),
    stderr=str(LOGDIR/"variant_calling/{sample}.split_read.{type}.stderr.log"),
  shadow:
    "shallow"
  conda:
    "envs/main.yml"
  threads: 10
  shell:
    """
    seqtk subseq \
      {input.filtered_fastq} \
      {input.SMN_list} \
      > {output.SMN_fastq} \
      2> {log.stderr}
    """

rule alignment:
  input:
    SMN_fastq=OUTDIR/"{sample}/{sample}.filtered.{type}.fq".format(sample="{sample}", type="{type}"),
  output:
    bam=OUTDIR/"{sample}/{sample}.{type}.bam" if config["keep_output"] else temp(OUTDIR/"{sample}/{sample}.{type}.bam"),
  log:
    stdout=str(LOGDIR/"variant_calling/{sample}.alignment.{type}.stdout.log"),
    stderr=str(LOGDIR/"variant_calling/{sample}.alignment.{type}.stderr.log"),
  shadow:
    "shallow"
  conda:
    "envs/main.yml"
  threads: 10
  params:
    reference=lambda w: INPUTDIR/"{type}_amplicon.fa".format(type=w.type),
  shell:
    """
    minimap2 \
      -a \
      {params.reference} \
      {input.SMN_fastq} \
      2> {log.stderr} | \
    samtools sort \
      -o {output.bam} \
    """

rule variant_calling:
  input:
    bam=OUTDIR/"{sample}/{sample}.{type}.bam".format(sample="{sample}", type="{type}"),
  output:
    vcf=OUTDIR/"{sample}/{type}/variant_calls.final.vcf.gz"
  log:
    stdout=str(LOGDIR/"variant_calling/{sample}.variant_calling.{type}.stdout.log"),
    stderr=str(LOGDIR/"variant_calling/{sample}.variant_calling.{type}.stderr.log"),
  shadow:
    "shallow"
  conda:
    "envs/main.yml"
  threads: 10
  params:
    outdir=directory(OUTDIR/"{sample}/{type}"),
    reference=lambda w: INPUTDIR/"{type}_amplicon.fa".format(type=w.type),
  shell:
    """
    samtools index {input.bam}
    NanoCaller \
      -bam {input.bam} \
      -p ont \
      -o {params.outdir} \
      -ref {params.reference} \
      -chrom {wildcards.type} \
      -cpu {threads} \
      > {log.stdout} \
      2> {log.stderr}
    """

rule variant_annotation:
  input:
    vcf=OUTDIR/"{sample}/{type}/variant_calls.final.vcf.gz".format(sample="{sample}", type="{type}"),
  output:
    vep=OUTDIR/"{sample}/{type}.variant_calls.final.vep",
  log:
    stdout=str(LOGDIR/"variant_calling/{sample}.variant_annotation.{type}.stdout.log"),
    stderr=str(LOGDIR/"variant_calling/{sample}.variant_annotation.{type}.stderr.log"),
  shadow:
    "shallow"
  conda:
    "envs/main.yml"
  threads: 10
  params:
    sample="{sample}",
    # type=lambda w: "{type}".format(type=w.type),
    type="{type}",
    path=OUTDIR,
  shell:
    """
    python utils/run_variant_annotation.py \
      --path {params.path} \
      --sample {params.sample} \
      --type {params.type} \
      > {log.stdout} \
      2> {log.stderr}
    """

#############################
# Final Report
#############################
rule final_report:
  input:
    OUTDIR/"{sample}/SMN1.variant_calls.final.vep",
    OUTDIR/"{sample}/SMN2.variant_calls.final.vep",
    OUTDIR/"{sample}/{sample}.cnv.png",
    OUTDIR/"{sample}/{sample}.stat",
    OUTDIR/"{sample}/{sample}.filtered.stat",
  output:
    report=OUTDIR/"{sample}/report.html",
  log:
    stdout=str(LOGDIR/"final_report/{sample}.stdout.log"),
    stderr=str(LOGDIR/"final_report/{sample}.stderr.log"),
  shadow:
    "shallow"
  conda:
    "envs/main.yml"
  threads: 10
  params:
    sample="{sample}",
    path=OUTDIR
  shell:
    """
    python utils/generate_report.py \
      --path {params.path} \
      --sample {params.sample} \
      > {log.stdout} \
      2> {log.stderr}
    """