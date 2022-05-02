# ViralCC: leveraging metagenomic proximity-ligation to retrieve complete viral genomes and reveal active co-host systems

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [A test dataset to demo ViralCC](#a-test-dataset-to-demo-viralcc)
- [Instruction to process raw data](#instruction-to-process-raw-data)
- [Instruction to run ViralCC](#instruction-to-run-viralcc)
- [Instruction of reproducing results in ViralCC paper](https://github.com/dyxstat/Reproduce_ViralCC)
- [Contacts and bug reports](#contacts-and-bug-reports)
- [Copyright and License Information](#copyright-and-license-information)
- [Issues](https://github.com/dyxstat/ViralCC/issues)

# Overview
`ViralCC` is a new open-source metagenomic Hi-C-based binning pipeline to recover high-quality viral genomes. 
`ViralCC` not only considers the Hi-C interaction graph, but also puts forward a novel host proximity graph of viral contigs 
as a complementary source of information to the remarkably sparse Hi-C interaction map. The two graphs are then integrated together, 
followed by the Leiden graph clustering using the integrative graph to generate draft viral genomes.

**If you want to reproduce results in our ViralCC paper, please read our instructions [here](https://github.com/dyxstat/Reproduce_ViralCC).**

# System Requirements
## Hardware requirements
`ViralCC` requires only a standard computer with enough RAM to support the in-memory operations.

## Software requirements
### OS Requirements
`ViralCC` v1.0.0 is supported and tested in *MacOS* and *Linux* systems.

### Python Dependencies
`ViralCC` mainly depends on the Python scientific stack.

```
numpy
scipy
pysam
scikit-learn
pandas
Biopython
leidenalg
```


# Installation Guide
We recommend using [**conda**](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html) to install `ViralCC`. 
Typical installation time is 1-5 minutes depending on your system.

### Clone the repository with git:
```
git clone https://github.com/dyxstat/ViralCC.git
```

Once complete, enter the repository folder and then create a `ViralCC` environment using conda.


### Enter the ViralCC folder:
```
cd ViralCC
```

### Construct environment in the linux or MacOS system:
```
conda env create -f viralcc_linux_env.yaml
```
or
```
conda env create -f viralcc_osx_env.yaml
```

### Enter the environment:
```
conda activate ViralCC_env
```

# A test dataset to demo ViralCC
We provide a small simulated dataset, located under the Test directory, to demo and test the software:
```
Test/final.contigs.fa
Test/MAP_SORTED.bam
Test/viral_contigs.txt
```
Run `ViralCC` on the testing dataset:
```
python ./viralcc.py pipeline -v Test/final.contigs.fa Test/MAP_SORTED.bam Test/viral_contigs.txt Test/out_test
```

The expected run time for demo is several seconds and the expected output are in the 'Test/out_test' directory:
```
Test/out_test/cluster_viral_contig.txt
Test/out_test/prokaryotic_contig_info.csv
Test/out_test/VIRAL_BIN/VIRAL_BIN0000.fa
Test/out_test/VIRAL_BIN/VIRAL_BIN0001.fa
Test/out_test/viralcc.log
Test/out_test/viral_contig_info.csv
```



# Instruction to process raw data
Follow the instructions in this section to process the raw shotgun and Hi-C data and generate the input for `ViralCC`:

### Clean raw shotgun and Hi-C reads

Adaptor sequences are removed by `bbduk` from the `BBTools` suite with parameter `ktrim=r k=23 mink=11 hdist=1 minlen=50 tpe tbo` and reads are quality-trimmed using `bbduk` with parameters `trimq=10 qtrim=r ftm=5 minlen=50`. Additionally, the first 10 nucleotides of Hi-C reads are trimmed by `bbduk` with parameter `ftl=10`. Identical PCR optical and tile-edge duplicates for Hi-C reads were removed by the script `clumpify.sh` from `BBTools` suite.

### Assemble shotgun reads

For the shotgun library, de novo metagenome assembly is produced by an assembly software, such as MEGAHIT.
```
megahit -1 SG1.fastq.gz -2 SG2.fastq.gz -o ASSEMBLY --min-contig-len 1000 --k-min 21 --k-max 141 --k-step 12 --merge-level 20,0.95
```

### Align Hi-C paired-end reads to assembled contigs

Hi-C paired-end reads are aligned to assembled contigs using a DNA mapping software, such as BWA MEM. Then, samtools with parameters ‘view -F 0x904’ is applied to remove unmapped reads, supplementary alignments, and secondary alignments. BAM file needs to be sorted **by name** using 'samtools sort'.
```
bwa index final.contigs.fa
bwa mem -5SP final.contigs.fa hic_read1.fastq.gz hic_read2.fastq.gz > MAP.sam
samtools view -F 0x904 -bS MAP.sam > MAP_UNSORTED.bam
samtools sort -n MAP_UNSORTED.bam -o MAP_SORTED.bam
```
### Identify viral contigs from assembled contigs

Assembled contigs were screened by a viral sequence detection software, such as VirSorter to identify viral contigs.
```
wrapper_phage_contigs_sorter_iPlant.pl -f final.contigs.fa --db 1 --wdir virsorter_output --data-dir virsorter-data
```


# Instruction to run ViralCC
```
python ./viralcc.py pipeline [Parameters] FASTA_file BAM_file VIRAL_file OUTPUT_directory
```
### Parameters
```
--min-len: Minimum acceptable contig length (default 1000)
--min-mapq: Minimum acceptable alignment quality (default 30)
--min-match: Accepted alignments must be at least N matches (default 30)
--min-k: Lower bound of k for determining the host poximity graph (default 4)
--random-seed: Random seed for the Leiden clustering (default 42)
--cover (optional): Cover existing files. Otherwise, an error will be returned if the output file is detected to exist.
-v (optional): Verbose output about more specific details of the ViralCC procedure.
```
### Input File

* FASTA_file: a fasta file of the assembled contig (e.g. Test/final.contigs.fa)
* BAM_file: a bam file of the Hi-C alignment (e.g. Test/MAP_SORTED.bam)
* VIRAL_file: a txt file containing the names of identified viral contigs in one column **without header** (e.g. Test/viral_contigs.txt)


### Output File

* VIRAL_BIN: folder containing the fasta files of draft viral bins
* cluster_viral_contig.txt: clustering results with 2 columns, the first is the viral contig name, and the second is the group number.
* viral_contig_info.csv: information of viral contigs with three columns (contig name, contig length, and GC-content)
* prokaryotic_contig_info.csv: information of non-viral contigs with three columns (contig name, contig length, and GC-content)
* viralcc.log: log file of ViralCC

### Example
```
python ./viralcc.py pipeline -v final.contigs.fa MAP_SORTED.bam viral_contigs.txt out_test
```

# Contacts and bug reports
If you have any questions or suggestions, welcome to contact Yuxuan Du (yuxuandu@usc.edu).


# Copyright and License Information
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.






