import gzip
import io
import os
import sys
import subprocess


def make_dir(path, exist_ok=False):
    """
    Convenience method for making directories with a standard logic.
    An exception is raised when the specified path exists and is not a directory.
    :param path: target path to create
    :param exist_ok: if true, an existing directory is ok. Existing files will still cause an exception
    """
    if not os.path.exists(path):
        os.mkdir(path)
    elif not exist_ok:
        raise IOError('output directory already exists!')
    elif os.path.isfile(path):
        raise IOError('output path already exists and is a file!')




def count_fasta_sequences(file_name):
    """
    Estimate the number of fasta sequences in a file by counting headers. Decompression is automatically attempted
    for files ending in .gz. Counting and decompression is by why of subprocess calls to grep and gzip. Uncompressed
    files are also handled. This is about 8 times faster than parsing a file with BioPython and 6 times faster
    than reading all lines in Python.

    :param file_name: the fasta file to inspect
    :return: the estimated number of records
    """
    if file_name.endswith('.gz'):
        proc_uncomp = subprocess.Popen(['gzip', '-cd', file_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc_read = subprocess.Popen(['grep', r'^>'], stdin=proc_uncomp.stdout, stdout=subprocess.PIPE)
    else:
        proc_read = subprocess.Popen(['grep', r'^>', file_name], stdout=subprocess.PIPE)

    n = 0
    for _ in proc_read.stdout:
        n += 1
    return n


def gen_bins(fastafile,resultfile,outputdir):
    # read fasta file
    sequences={}
    if fastafile.endswith("gz"):
        with gzip.open(fastafile,'r') as f:
            for line in f:
                line=str(line,encoding="utf-8")
                if line.startswith(">"):
                    if " " in line:
                        seq,others=line.split(' ', 1)
                        sequences[seq] = ""
                    else :
                        seq=line.rstrip("\n")
                        sequences[seq] = ""
                else:
                    sequences[seq] += line.rstrip("\n")
    else:
        with open(fastafile,'r') as f:
            for line in f:
                if line.startswith(">"):
                    if " " in line:
                        seq,others=line.split(' ', 1)
                        sequences[seq] = ""
                    else :
                        seq=line.rstrip("\n")
                        sequences[seq] = ""
                else:
                    sequences[seq] += line.rstrip("\n")
    dic={}
    with open(resultfile,"r") as f:
        for line in f:
            contig_name,cluster_name=line.strip().split('\t')
            try:
                dic[cluster_name].append(contig_name)
            except:
                dic[cluster_name]=[]
                dic[cluster_name].append(contig_name)
    print("Writing bins in \t{}".format(outputdir))
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    
    bin_name=0
    for _,cluster in dic.items():
        if bin_name < 10:
            bin = 'VIRAL_BIN'+ '000' + str(bin_name) + '.fa'
        elif bin_name >= 10 and bin_name < 100:
            bin = 'VIRAL_BIN'+ '00' + str(bin_name) + '.fa'
        elif bin_name >= 100 and bin_name < 1000:
            bin = 'VIRAL_BIN'+ '0' + str(bin_name) + '.fa'
        else:
            bin = 'VIRAL_BIN'+str(bin_name) + '.fa'
        binfile=os.path.join(outputdir,"{}".format(bin))
        with open(binfile,"w") as f:
            for contig_name in cluster:
                contig_name=">"+contig_name
                try:
                    sequence=sequences[contig_name]
                except:
                    continue
                f.write(contig_name+"\n")
                f.write(sequence+"\n")
        bin_name+=1

