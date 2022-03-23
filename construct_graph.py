#!/usr/bin/env python
# coding: utf-8

from collections import OrderedDict, namedtuple
import Bio.SeqIO as SeqIO
from Bio.SeqUtils import GC
from itertools import combinations
from math import exp, log, sqrt
import pandas as pd
import numpy as np
import pysam
import scipy.sparse as scisp
import tqdm
import csv
import os
from utils import count_fasta_sequences, open_input
import logging

# package logger
logger = logging.getLogger(__name__)

#######offset is the enumeration of the length######################
#localid is the index of the list#
#refid is the index of the fasta file, which is a global index# 
SeqInfo = namedtuple('SeqInfo', ['localid', 'refid', 'name', 'length', 'GC']) #create a new class of tuple: SeqInfo

   

class Sparse2DAccumulator(object):
######create a 2D coo sparse matrix###########
    def __init__(self, N):
        self.shape = (N, N)
        self.mat = {}
        ##mat is a dictionary here
        # fixed counting type
        self.dtype = np.uint32

    def setitem(self, index, value):
        assert len(index) == 2 and index[0] >= 0 and index[1] >= 0, 'invalid index: {}'.format(index)
        assert isinstance(value, (int, np.int)), 'values must be integers'
        self.mat[index] = value

    def getitem(self, index):
        if index in self.mat:
            return self.mat[index]
        else:
            return 0

    def get_coo(self, symm=True):
        """
        Create a COO format sparse representation of the accumulated values.

        :param symm: ensure matrix is symmetric on return
        :return: a scipy.coo_matrix sparse matrix
        """
        coords = [[], []]
        data = []
        m = self.mat
        for i, j in m.keys(): ##m.keys() will return a tuple of two values
            coords[0].append(i)
            coords[1].append(j)
            data.append(m[i, j])

        m = scisp.coo_matrix((data, coords), shape=self.shape, dtype=self.dtype)

        if symm:
            m += scisp.tril(m.T, k=-1)

        return m.tocoo()



class ContactMatrix:

    def __init__(self, bam_file, seq_file, viral_file,  path , min_mapq=30, min_len=1000, min_match=30, min_k=4):

        ########input the parameter################
        ########################################################################
        ########################################################################
        #bam_file: alignment info of Hi-C library on contigs in bam#
        #seq_file: store the assembly contigs in fasta#
        #min_mapq: minimal mapping quality(default 30)#
        #min_len: minimal length of contigs(default 1000bp)

        self.bam_file = bam_file
        self.seq_file = seq_file
        self.viral_file = viral_file
        self.path = path
        self.min_mapq = min_mapq
        self.min_len = min_len
        self.min_match = min_match

        #fasta_info store the info from fasta file#
        #seq_info store the information of contigs from bam file#
        #cov_info store the information of coverage#
        #seq_map store the contact map#
        self.fasta_info = {}
        self.seq_info = []
        self.seq_map = None
        self.map_viral = None
        self.map_row = None
        self.viral_same_host = None
        self.map_combine = None
        self.viral_info = []
        self.prokaryotic_info = []
        

        #total reads: number of alignment in bam file#
        self.total_reads = None
        
        logger.info('Reading fasta file...')
        with open_input(seq_file) as multi_fasta:
            # fasta_info is a dictionary of preparation of seq_info and seq_info is the true results
            # get an estimate of sequences for progress
            fasta_count = count_fasta_sequences(seq_file)
            for seqrec in tqdm.tqdm(SeqIO.parse(multi_fasta, 'fasta'), total=fasta_count , desc='Analyzing both prokaryotic and viral contigs in reference fasta'):
                if len(seqrec) < min_len:
                    continue
                self.fasta_info[seqrec.id] = {'length': len(seqrec), 'GC': GC(seqrec.seq)}

        logger.debug('There are totally {} contigs in reference fasta'.format(fasta_count))
        
        # now parse the header information from bam file
        ###########input the bam data###############
        #########seq_info contain the global contig information and we don't change the seq_info after being created##############
        with pysam.AlignmentFile(bam_file, 'rb') as bam:
        ##test that BAM file is the correct sort order
            if 'SO' not in bam.header['HD'] or bam.header['HD']['SO'] != 'queryname':
                raise IOError('BAM file must be sorted by read name')

            # determine the set of active sequences
            # where the first filtration step is by length
            logger.info('Filtering contigs by minimal length({})...'.format(self.min_len))
            ref_count = {'seq_missing': 0, 'too_short': 0}

            offset = 0
            localid = 0
            
            for n, (rname, rlen) in enumerate(zip(bam.references, bam.lengths)):
            # minimum length threshold
                if rlen < min_len:
                    ref_count['too_short'] += 1
                    continue

                try:
                    fa = self.fasta_info[rname]
    
                except KeyError:
                    logger.info('Contig "{}" was not present in reference fasta'.format(rname))
                    ref_count['seq_missing'] += 1
                    #self.missing.append(rname)
                    continue

                assert fa['length'] == rlen, 'Sequence lengths in {} do not agree: bam {} fasta {}'.format(rname, fa['length'], rlen)

                self.seq_info.append(SeqInfo(localid , n , rname, rlen, fa['GC']))
                localid = localid + 1
                offset += rlen
            
            self.total_len = offset
            self.total_seq = localid
            
            del self.fasta_info, offset, localid

            if self.total_seq == 0:
                logger.info('No sequences in BAM found in FASTA')
                raise ImportError('No sequences in BAM found in FASTA')
            
            logger.debug('{} contigs miss and {} contigs are too short'.format(ref_count['seq_missing'] , ref_count['too_short']))
            logger.debug('Accepted {} contigs covering {} bp'.format(self.total_seq, self.total_len))
            
            
            logger.info('Counting reads in bam file...')
            self.total_reads = bam.count(until_eof=True)
            logger.debug('BAM file contains {0} alignments'.format(self.total_reads))

            logger.info('Handling the alignments...')
            self._bin_map(bam)
                

        assert self.seq_map.shape[0] == len(self.seq_info), 'Filter error'
        
        #########fiter viral contig and corresponding viral contacts#########
        viral_contig_name = pd.read_csv(viral_file, sep='\t' , header=None).values[: , 0]
        
        viral_contig_index = []
        prokaryotic_contig_index = []
        #######Distinguish the index of viral contigs and prokaryotic contigs###########
        for i in range(len(self.seq_info)):
            contig_temp = self.seq_info[i]
            if contig_temp.name in viral_contig_name:
                viral_contig_index.append(i)
            else:
                prokaryotic_contig_index.append(i)
        
        del contig_temp
        
        for i , idn in enumerate(viral_contig_index):
            viral_seq = self.seq_info[idn]
            assert viral_seq.localid == idn, 'the local index does not match the contact matrix index'
            self.viral_info.append(SeqInfo(i , viral_seq.refid , viral_seq.name , viral_seq.length , viral_seq.GC))

        for i , idn in enumerate(prokaryotic_contig_index):
            prokaryotic_seq = self.seq_info[idn]
            assert prokaryotic_seq.localid == idn, 'the local index does not match the contact matrix index'
            self.prokaryotic_info.append(SeqInfo(i , prokaryotic_seq.refid , prokaryotic_seq.name , prokaryotic_seq.length , prokaryotic_seq.GC))
            
        del self.seq_info
        
        logger.info('There are {} viral contigs'.format(len(viral_contig_index)))
        logger.info('There are {} potential host contigs'.format(len(prokaryotic_contig_index)))
        logger.info('Write information of viral contigs and potential host contigs')
        self._write_viral_info()
        self._write_prokaryotic_info()
        
        
        self.map_viral = self.seq_map.tocsr()
        self.map_viral = self.map_viral[viral_contig_index , :]
        self.map_viral = self.map_viral.tocsc()
        self.map_viral = self.map_viral[: , viral_contig_index]
        self.map_viral = self.map_viral.todense()
        
        self.map_row = self.seq_map.tocsr()
        self.map_row = self.map_row[viral_contig_index , :]
        self.map_row = self.map_row.tocsc()
        self.map_row = self.map_row[: , prokaryotic_contig_index]
        self.map_row = self.map_row.tocsr()
        
        del self.seq_map
        
        numV = len(self.viral_info)
        numE = np.count_nonzero(self.map_viral)/2
        
        for k in range(min_k , 50):
            count = 0
            self.viral_same_host = np.zeros((numV,numV))
            for i in range(numV):
                for j in range(i+1 , numV):
                    host_i = self.map_row[i , :].tocoo().col
                    if host_i.shape[0] < 1:
                        break
                    host_j = self.map_row[j , :].tocoo().col
                    if host_j.shape[0] < 1:
                        break
                    if np.intersect1d(host_i , host_j).shape[0] >= k:
                        count += 1
                        self.viral_same_host[i , j] = 1
                        self.viral_same_host[j , i] = 1
            if count <= numE:
                break
        
        logger.info('the threshold of shared host contig is {}'.format(k))
        logger.info('there are {} edges in the host proximity graph'.format(count))
        logger.info('there are {} edges in the Hi-C interaction graph'.format(numE))
        
        
        self.map_combine = np.zeros((numV,numV))
        map_temp = scisp.coo_matrix(self.map_viral+self.viral_same_host)
        sources = map_temp.row
        targets = map_temp.col
        for i,j in zip(sources , targets):
            self.map_combine[i , j] = 1
        
        logger.info('Integrate the Hi-C interaction graph and the host proximity graph')
        logger.info('Integrative graph construction finished and there are {} edges in the integrative graph'.format(np.count_nonzero(self.map_combine)/2))
        del self.map_viral, self.map_row, prokaryotic_contig_index, viral_contig_index, self.viral_same_host, map_temp
        
        
    def _bin_map(self, bam):
        """
        Accumulate read-pair observations from the supplied BAM file.
        Maps are initialized here. Logical control is achieved through initialisation of the
        ContactMap instance, rather than supplying this function arguments.

        :param bam: this instance's open bam file.
        """
        import tqdm

        def _simple_match(r):
            return r.mapping_quality >= _mapq

        def _strong_match(r):
            if r.mapping_quality < _mapq or r.cigarstring is None:
                return False
            cig = r.cigartuples[-1] if r.is_reverse else r.cigartuples[0]
            return cig[0] == 0 and cig[1] >= self.min_match

        # set-up match call
        _matcher = _strong_match if self.min_match else _simple_match

        def next_informative(_bam_iter, _pbar):
            while True:
                r = next(_bam_iter)
                _pbar.update()
                if not r.is_unmapped and not r.is_secondary and not r.is_supplementary:
                    return r
       
        _seq_map = Sparse2DAccumulator(self.total_seq)


        with tqdm.tqdm(total=self.total_reads) as pbar:

            # locals for read filtering
            _mapq = self.min_mapq

            _idx = self.make_reverse_index('refid') #from global index to local index#
            _len = bam.lengths
     
            counts = OrderedDict({
                'accepted pairs': 0,
                'map_same_contig pairs': 0,
                'ref_excluded pairs': 0,
                'poor_match pairs': 0,
                'single read':0})

            bam.reset()
            bam_iter = bam.fetch(until_eof=True)
            self.index1 = 0
            while True:
                self.index1 += 1
                try:
                    r1 = next_informative(bam_iter, pbar)
                    while True:
                        # read records until we get a pair
                        r2 = next_informative(bam_iter, pbar)
                        if r1.query_name == r2.query_name:
                            break
                        r1 = r2 ###if we don't get a pair, next _bam_iter
                        counts['single read'] += 1
                        
                except StopIteration:
                    break

                if r1.reference_id not in _idx or r2.reference_id not in _idx:
                    counts['ref_excluded pairs'] += 1
                    continue

                if r1.reference_id == r2.reference_id:
                    counts['map_same_contig pairs'] += 1
                    continue

                if not _matcher(r1) or not _matcher(r2):
                    counts['poor_match pairs'] += 1
                    continue

                # get internal indices
                ix1 = _idx[r1.reference_id]
                ix2 = _idx[r2.reference_id]

                # maintain just a half-matrix
                if ix2 < ix1:
                    ix1, ix2 = ix2, ix1

                counts['accepted pairs'] += 1
                ix = (ix1 , ix2)
                if _seq_map.getitem(ix):
                    temp_value = _seq_map.getitem(ix) + 1
                    _seq_map.setitem(ix , temp_value)
                else:
                    _seq_map.setitem(ix , 1)
                
        self.seq_map = _seq_map.get_coo()
        del _seq_map, r1, r2, _idx

        logger.debug('Pair accounting: {}'.format(counts))


    def make_reverse_index(self, field_name):
        """
        Make a reverse look-up (dict) from the chosen field in seq_info to the internal index value
        of the given sequence. Non-unique fields will raise an exception.

        :param field_name: the seq_info field to use as the reverse.
        :return: internal array index of the sequence
        """
        rev_idx = {}
        for n, seq in enumerate(self.seq_info):
            fv = getattr(seq, field_name)
            if fv in rev_idx:
                raise RuntimeError('field contains non-unique entries, a 1-1 mapping cannot be made')
            rev_idx[fv] = n
        return rev_idx



    def _write_prokaryotic_info(self):
        with open(os.path.join(self.path , 'prokaryotic_contig_info.csv'),'w') as out:
            for seq in self.prokaryotic_info:
                out.write(str(seq.name)+ ',' +str(seq.length)+ ',' +str(seq.GC))
                out.write('\n')
                
    def _write_viral_info(self):
        with open(os.path.join(self.path , 'viral_contig_info.csv'),'w') as out:
            for seq in self.viral_info:
                out.write(str(seq.name)+ ',' +str(seq.length)+ ',' +str(seq.GC))
                out.write('\n')


