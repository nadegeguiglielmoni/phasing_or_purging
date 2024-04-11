# Phasing or purging: tackling the genome assembly of a highly heterozygous animal species in the era of high-accuracy long reads

## Quality control

```sh
jellyfish count -C -m 21 -s 50000000000 hifi_reads.fastq -o hifi.jf
jellyfish histo hifi.jf > hifi.histo

genomescope.R -i hifi.histo -o genomescope_k21 -k 21
```

```sh
NanoPlot --fastq ont_reads.fastq.gz -c blue --tsv_stats -o nanoplot_ont
NanoPlot --fastq hifi_reads.fastq.gz -c red --tsv_stats -o nanoplot_hifi
```

## *De novo* assembly

```sh
nextdenovo run.cfg
```
Nanopore Q20+ reads as ONT corrected
```sh
[General]
job_type = local
job_prefix = nextDenovo
task = all
rewrite = yes
deltmp = yes
parallel_jobs = 20
input_type = corrected
read_type = ont # clr, ont, hifi
input_fofn = input.fofn
workdir = nextdenovo_v25_hifi_default_g250Mb_ont_q20_reads

[correct_option]
read_cutoff = 1k
genome_size = 250m # estimated genome size
sort_options = -m 10g -t 20
minimap2_options_raw = -t 8
pa_correction = 5
correction_options = -p 30

[assemble_option]
minimap2_options_cns = -t 8
nextgraph_options = -a 1
```

Nanopore Q20+ reads as HiFi
```sh
[General]
job_type = local
job_prefix = nextDenovo
task = all
rewrite = yes
deltmp = yes
parallel_jobs = 20
input_type = raw
read_type = hifi # clr, ont, hifi
input_fofn = input.fofn
workdir = nextdenovo_v25_hifi_default_g250Mb_ont_q20_reads

[correct_option]
read_cutoff = 1k
genome_size = 250m # estimated genome size
sort_options = -m 10g -t 20
minimap2_options_raw = -t 8
pa_correction = 5
correction_options = -p 30

[assemble_option]
minimap2_options_cns = -t 8
nextgraph_options = -a 1
```
