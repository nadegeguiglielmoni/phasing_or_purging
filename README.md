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

[Flye](https://github.com/fenderglass/Flye) v2.9
```sh
flye -o flye_v29_nano-hq_ont --nano-hq ont_reads.fastq.gz
flye -o flye_v29_nano-corr_ont --nano-corr ont_reads.fastq.gz
flye -o flye_v29_nano-corr_ont-q20 --nano-corr ont_q20_reads.fastq.gz
flye -o flye_v29_pacbio-hifi_ont-q20 --pacbio-hifi ont_q20_reads.fastq.gz
flye -o flye_v29_pacbio-hifi_hifi --pacbio-hifi hifi_reads.fastq.gz
```

[hifiasm](https://github.com/chhylp123/hifiasm) v0.19.4
```sh
hifiasm -l 0 -o hifiasm_out long_reads.fastq.gz
hifiasm -l 3 -o hifiasm_out long_reads.fastq.gz
awk '/^S/{print ">"$2;print $3}' hifiasm_out.bp.p_ctg.gfa > hifiasm_out.bp.p_ctg.fasta
```
```sh
filtlong --min_length 30000 ont_reads.fastq.gz 
hifiasm -l 0 -o hifiasm_out --ul ont_reads.min30kb.fastq.gz hifi_reads.fastq.gz
hifiasm -l 3 -o hifiasm_out --ul ont_reads.min30kb.fastq.gz hifi_reads.fastq.gz
awk '/^S/{print ">"$2;print $3}' hifiasm_out.bp.p_ctg.gfa > hifiasm_out.bp.p_ctg.fasta
```

[NextDenovo](https://github.com/Nextomics/NextDenovo) v2.5

```sh
ls long_reads.fastq.gz > input.fofn
nextdenovo run.cfg
```
Nanopore reads or Nanopore Q20+ reads as ONT corrected
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
workdir = nextdenovo_v25

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

Nanopore Q20+ reads as HiFi or PacBio HiFi reads
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
workdir = nextdenovo_v25

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

[PECAT](https://github.com/lemene/PECAT) v0.0.3

Nanopore or Nanopore Q20+ reads
```sh
project=pecat_v003_cfg_ont_g250m
reads=ont_reads.fastq.gz
genome_size= 250000000
threads=20
cleanup=0
compress=0
grid=auto
prep_min_length=3000
prep_output_coverage=80
corr_iterate_number=1
corr_block_size=4000000000
corr_filter_options=--filter0=l=5000:al=2500:alr=0.5:aal=8000:oh=3000:ohr=0.3
corr_correct_options=--score=weight:lc=10 --aligner edlib:bs=1000:mc=6  --min_coverage 4 --filter1 oh=1000:ohr=0.01 --candidate n=600:f=30 --min_identity 90 --min_local_identity 80
corr_rd2rd_options=-x ava-ont -f 0.005 -I 10G
corr_output_coverage=60
align_block_size=12000000000
align_rd2rd_options=-X -g3000 -w30 -k19 -m100 -r500 -I 10G -f 0.002
align_filter_options=--filter0=l=5000:aal=6000:aalr=0.5:oh=3000:ohr=0.3 --task=extend --filter1=oh=300:ohr=0.03 --min_identity 0.95
asm1_assemble_options=--max_trivial_length 10000

phase_method=2
phase_rd2ctg_options=-x map-ont -w 10 -k19  -c -p 0.5 -r 1000 -I 10G -K 8G
phase_use_reads=1
phase_phase_options= --coverage lc=20 --phase_options icr=0.1:icc=8:sc=10
phase_filter_options = --threshold 1000

phase_clair3_command = docker run -i -v `pwd -P`:`pwd -P` -v /tmp:/tmp hkubal/clair3:latest /opt/bin/run_clair3.sh
phase_clair3_rd2ctg_options=-x map-ont -w10 -k19 -c -p 0.5 -r 1000 -I 10G -K 8G
phase_clair3_use_reads=0
phase_clair3_phase_options= --coverage lc=20 --phase_options icr=0.1:icc=3:sc=10 --filter i=90
phase_clair3_filter_options = --threshold 2500 --rate 0.05
phase_clair3_options=--platform=ont --model_path=/opt/models/r941_prom_sup_g5014  --include_all_ctgs

asm2_assemble_options=--reducer0 "best:cmp=2,0.1,0.1|phase:sc=3" --contig_format prialt,dual

polish_map_options = -x map-ont -w10 -k19 -I 10g -K 8G -a
polish_cns_options =
polish_use_reads=0
polish_filter_options=--filter0 oh=2000:ohr=0.2:i=96

polish_medaka = 1
polish_medaka_command = singularity exec -B `pwd -P`:`pwd -P`  -B /tmp:/tmp medaka_v1.7.2.sif medaka
polish_medaka_map_options = -x map-ont -w10 -k19 -I 10g -K 8G
polish_medaka_cns_options = --model r1041_e82_260bps_sup_g632
polish_medaka_filter_options=--filter0 oh=2000:ohr=0.2:i=96
```

PacBio HiFi
```sh
project=pecat_v003_cfg_hifi_g250m
reads=hifi_reads.fastq.gz
genome_size= 250000000
threads=30
cleanup=0
grid=auto
prep_min_length=3000
prep_output_coverage=60
corr_iterate_number=1
corr_block_size=4000000000
corr_correct_options=--score=weight:lc=8 --aligner diff:s=100 --min_coverage 1 --filter1 oh=100 --min_identity 96 --min_local_identity 95
corr_filter_options=--filter0=l=5000:al=2500:alr=0.5:aal=5000:oh=3000:ohr=0.3
corr_rd2rd_options=-X -g3000 -w30 -k19 -m100 -r500 -f 0.002 -K8G -I 8G
corr_output_coverage=60
align_block_size=4000000000
align_rd2rd_options=-X -g3000 -w30 -k19 -m100 -r500 -f 0.002 -K8G -I 8G
align_filter_options=--filter0=l=5000:aal=6000:aalr=0.5:oh=3000:ohr=0.3 --task=extend --filter1=oh=50:ohr=0.01  --aligner diff:s=100 --min_identity 0.90
asm1_assemble_options= --min_identity 0.99 --min_coverage 1

phase_method=2
phase_rd2ctg_options=-x map-hifi  -c -p 0.5 -r 1000 -K8G
phase_use_reads=1
phase_phase_options= --coverage lc=8 --phase_options icr=0.02:icc=3:sc=4 --filter=i=95.00:alr=0.80:oh=100:ohr=0.01:ilid=100

phase_clair3_command=docker run -i -v `pwd -P`:`pwd -P` -v /tmp:/tmp hkubal/clair3:latest /opt/bin/run_clair3.sh
phase_clair3_use_reads=0
phase_clair3_options=--platform=hifi --model_path=/opt/models/hifi  --include_all_ctgs
phase_clair3_rd2ctg_options=-x map-hifi  -c -p 0.5 -r 1000 -K8G
phase_clair3_phase_options=--coverage lc=8 --phase_options icr=0.02:icc=2:sc=4 --filter i=95
phase_clair3_filter_options=

asm2_assemble_options= --reducer0 "best:cmp=2,0.1,0.1|phase:sc=2" --min_identity 0.99 --max_trivial_length 10000 --contig_format dual,prialt --min_coverage 1

polish_map_options=-x map-hifi -I8G -K8G -a
polish_filter_options=--filter0 oh=500:ohr=0.05:i=98
polish_cns_options=
polish_medaka_command = singularity exec -B `pwd -P`:`pwd -P`  -B /tmp:/tmp medaka_v1.7.2.sif medaka
```

[Verkko](https://github.com/marbl/verkko) v0.0.3
```sh
verkko -d verkko_v141 --hifi hifi_reads.filtlong.fastq.gz --nano ont_reads.fastq.gz --min-ont-length 30000
```

## Decontamination

[BLAST](https://www.ncbi.nlm.nih.gov/geo/query/blast.html) v2.10.0
```sh
blastn -query assembly.fasta -db nt -outfmt "6 qseqid staxids bitscore std sscinames scomnames" \
               -max_hsps 1 -evalue 1e-25 -out assembly.fasta.blast.out
```

[minimap2](https://github.com/lh3/minimap2) v2.24
```sh
# for Nanopore reads
minimap2 -ax map-ont assembly.fasta ont_reads.trimmed.fastq.gz | samtools sort > minimap2.assembly.bam

# for HiFi reads
minimap2 -ax map-hifi assembly.fasta hifi_reads.fastq.gz | samtools sort > minimap2.assembly.bam
```

[BUSCO](https://gitlab.com/ezlab/busco) v5.4.7
```sh
docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v5.4.7_cv1 busco -i assembly.fasta -o busco_metazoa_odb10_assembly -m genome -l metazoa_odb10
```

[BlobToolKit](https://blobtoolkit.genomehubs.org/) v4.3.2
```sh
blobtools add --fasta assembly.fasta --cov minimap2.assembly.bam --hits assembly.fasta.blast.out --busco busco_metazoa_odb10_assembly/run_metazoa_odb10/full_table.tsv --taxdump taxdump --create blobdir_out
blobtools view blobdir_out
```

## Haplotig purging

```sh
# for Nanopore reads
minimap2 -x map-ont assembly.fasta ont_reads.trimmed.fastq.gz | gzip -c - > minimap2.assembly.paf.gz

# for HiFi reads
minimap2 -x map-hifi assembly.fasta hifi_reads.fastq.gz | gzip -c - > minimap2.assembly.paf.gz
        
pbcstat minimap2.assembly.paf.gz 
calcuts PB.stat > cutoffs 2>calcults.log

hist_plot.py -c cutoffs PB.stat purge_dups.${pri_asm}.png

split_fa assembly.fasta > assembly.fasta.split
minimap2 -xasm5 -DP assembly.fasta.split assembly.fasta.split | gzip -c - > assembly.fasta.split.self.paf.gz

purge_dups -2 -T cutoffs -c PB.base.cov assembly.fasta.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs dups.bed assembly.fasta
```

## Assembly evaluation

[assembly-stats](https://github.com/sanger-pathogens/assembly-stats) v1.0.1
```sh
assembly-stats assembly.fasta
```

[BUSCO](https://gitlab.com/ezlab/busco) v5.4.7
```sh
docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v5.4.7_cv1 busco -i assembly.fasta -o busco_metazoa_odb10_assembly -m genome -l metazoa_odb10
docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v5.4.7_cv1 busco -i assembly.fasta -o busco_nematoda_odb10_assembly -m genome -l nematoda_odb10
```

[KAT](https://github.com/TGAC/KAT) v2.4.2
```sh
kat comp -o kat_comp hifi_reads.fastq.gz assembly.fasta
```
