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

