# ViralCC: leveraging metagenomic proximity-ligation to retrieve complete viral genomes and reveal active co-host systems

## Introduction
ViralCC is a new open-source metagenomic Hi-C-based binning pipeline to recover complete viral genomes. 


## Install and Setup
### conda
We recommend using conda to run ViralCC.

After installing Anaconda (or miniconda), Users can clone the repository with git:
```
git clone https://github.com/dyxstat/ViralCC.git
```

Once complete, you can enter the repository folder and then create a HiCBin environment using conda:
```
# enter the HiCBin folder
cd HiCBin
# Construct environment
conda env create -f HiCBin_linux_env.yaml 
or
conda env create -f HiCBin_osx_env.yaml
# Enter the environment
conda activate HiCBin_env
```

Normalization method in HiCBin depends on R package '[glmmTMB](https://github.com/glmmTMB/glmmTMB)'. Though the R package can be installed by 'conda install -c conda-forge r-glmmtmb', you will meet one potential warning derived from the dependency version (https://github.com/glmmTMB/glmmTMB/issues/615): 'Error in .Call("FreeADFunObject", ptr, PACKAGE = DLL) : "FreeADFunObject" not available for .Call() for package "glmmTMB"' and we are not sure whether this warning would influence the noramlization results. To get rid of this warning, we strongly recommend you to install the source version of package 'glmmTMB' directly in R:

```
# Enter the R
R
# Download the R package and you may need to select a CRAN mirror for the installation
install.packages("glmmTMB", type="source")
```

Finally, you can test the pipeline, and testing result are in test/out/hicbin.log:
```
python ./hicbin.py test test/out
```

## Initial data preparation
### 1.Preprocess Raw reads
Adaptor sequences are removed by bbduk from the BBTools suite with parameter ‘ktrim=r k=23 mink=11 hdist=1 minlen=50 tpe tbo’ and reads are quality-trimmed using bbduk with parameters ‘trimq=10 qtrim=r ftm=5 minlen=50’. Then, the first 10 nucleotides of each read are trimmed by bbduk with parameter ‘ftl=10’.
### 2.Shotgun assembly
For the shotgun library, de novo metagenome assembly is produced by an assembly software, such as MEGAHIT.
```
megahit -1 SG1.fastq.gz -2 SG2.fastq.gz -o ASSEMBLY --min-contig-len 1000 --k-min 21 --k-max 141 --k-step 12 --merge-level 20,0.95
```
### 3.Align Hi-C paired-end reads to assembled contigs
Hi-C paired-end reads are mapped to assembled contigs using BWA-MEM with parameters ‘-5SP’. Then, samtools with parameters ‘view -F 0x904’ is applied to remove unmapped reads (0x4) and supplementary (0x800) and secondary (0x100) alignments and then sort BAM files by name.
```
bwa index final.contigs.fa
bwa mem -5SP final.contigs.fa hic_read1.fastq.gz hic_read2.fastq.gz > MAP.sam
samtools view -F 0x904 -bS MAP.sam > MAP_UNSORTED.bam
samtools sort MAP_UNSORTED.bam -o . > MAP_SORTED.bam
```
### 4.Assign taxonomy to contigs by TAXAassign
The taxonomic assignment of contigs was resolved using NCBI’s Taxonomy and its nt database by TAXAassign(v0.4) with parameters ‘-p -c 20 -r 10 -m 98 -q 98 -t 95 -a “60,70,80,95,95,98” -f’. 
### 5.Calculate the coverage of assembled contigs
Firstly, BBmap from the BBTools suite is applied to map the shotgun reads to the assembled contigs with parameters ‘bamscript=bs.sh; sh bs.sh’. The coverage of contigs is computed using script: ‘jgi summarize bam contig depths’ from MetaBAT2 v2.12.1.
```
bbmap.sh in1=SG1.fastq.gz in2=SG2.fastq.gz ref=final.contigs.fa out=SG_map.sam bamscript=bs.sh; sh bs.sh
jgi_summarize_bam_contig_depths --outputDepth coverage.txt --pairedContigs pair.txt SG_map_sorted.bam
```

## HiCBin analysis
### Implement the binning pipeline of HiCBin 
```
python ./hicbin.py pipeline [Parameters] FASTA_file BAM_file TAX_file COV_file OUTPUT_directory
```
#### Parameters
```
-e: Case-sensitive enzyme name. Use multiple times for multiple enzymes 
    (Optional; If no enzyme is input, HiCzin_LC mode will be automatically employed to do normalization)
--min-len: Minimum acceptable contig length (default 1000)
--min-mapq: Minimum acceptable alignment quality (default 30)
--min-match: Accepted alignments must be at least N matches (default 30)
--min-signal: Minimum acceptable signal (default 2)
--thres: Maximum acceptable fraction of incorrectly identified valid contacts in spurious contact detection (default 0.05)
--min-binsize: Minimum bin size used in output (default 150000)
--cover: Cover existing files
-v: Verbose output
```
#### Input File

* *FASTA_file*: a fasta file of the assembled contig (e.g. final.contigs.fa)
* *BAM_file*: a bam file of the Hi-C alignment (e.g. MAP_SORTED.bam)
* *TAX_file*: a csv file of contigs' taxonomy assignment by TAXAassign (e.g. contig_tax.csv)
* *COV_file*: a txt file of contigs' coverage information computed by script: ‘jgi summarize bam contig depths’ from MetaBAT2 (e.g. depth.txt)


#### Example
```
python ./hicbin.py pipeline -e HindIII -e NcoI -v final.contigs.fa MAP_SORTED.bam contig_tax.csv depth.txt out
```
If the restriction enzymes employed in the experiment are unspecified, use
```
python ./hicbin.py pipeline -v final.contigs.fa MAP_SORTED.bam contig_tax.csv depth.txt out
```
The results of the HiCBin software are all in the 'out' directory. The initial draft genomic bins are in 'out/BIN' and 'hicbin.log' file contains the specific implementation information of HiCBin.

### Implement the post-processing step of HiCBin
Initial draft genomic bins are assessed using [CheckM](https://github.com/Ecogenomics/CheckM).
Then the post-processing step of HiCBin is conducted for partially contaminated bins with completeness larger than 50% and contamination larger than 10%.
in order to purify the contaminated bins. Please make sure that the OUTPUT_directory is the same as your directory in the pipeline action.
```
python ./hicbin.py recluster --cover -v FASTA_file Contaminated_Bins_file OUTPUT_directory
```
#### Input File
* *FASTA_file*: a fasta file of the assembled contig (e.g. final.contigs.fa).
* *Contaminated_Bins_file*: a csv file of the names of the partially contaminated bins; Bin names are arranged in columns and *don't include the file formats .fa at the end of each name*

Example of a Contaminated_Bins_file:
```
BIN0000
BIN0001
BIN0005
...
```

#### Example
```
python ./hicbin.py recluster --cover -v final.contigs.fa contaminated_bins.csv out
```
The generated sub bins are in 'out/SUB_BIN' and the specific implementation details of the post-processing step will be added to the 'hicbin.log' file.

### Merge the generated bins from BIN and SUB_BIN
```
python ./merge_bins.py contaminated_bins.csv out
```
Then, all bins are merged in out/FINAL_BIN directory, which are the final bins generated by HiCBin.

## Contacts and bug reports
If you have any questions or suggestions, welcome to contact Yuxuan Du (yuxuandu@usc.edu).

## Citation
If you use our HiCBin software, please cite the following papers:

Du, Y., Laperriere, S.M., Fuhrman, J., Sun, F.: Normalizing Metagenomic Hi-C Data and Detecting Spurious Contacts Using Zero-Inflated Negative Binomial Regression. J Comput Biol. 29, 106–120 (2022). doi:10.1089/cmb.2021.0439

Du, Y., Sun, F. HiCBin: binning metagenomic contigs and recovering metagenome-assembled genomes using Hi-C contact maps. Genome Biol 23, 63 (2022). https://doi.org/10.1186/s13059-022-02626-w

## Acknowlegement
We employ some coding sections from bin3C [1] and solidbin [2] in our pipeline. If there are some interest conflicts, please contact Yuxuan Du (yuxuandu@usc.edu).

[1] DeMaere, M., Darling, A. bin3C: exploiting Hi-C sequencing data to accurately resolve metagenome-assembled genomes. Genome Biology 20, 46 (2019).

[2] Wang, Ziye and Wang, Zhengyang and Lu, Yang Young and Sun, Fengzhu and Zhu, Shanfeng. SolidBin: improving metagenome binning with semi-supervised normalized cut. Bioinformatics 35, 21 (2019).


## Copyright and License Information
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.






