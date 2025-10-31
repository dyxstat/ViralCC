#!/usr/bin/env python
# coding: utf-8

from collections import OrderedDict, namedtuple
from collections.abc import Iterable
from Bio.Restriction import Restriction
import Bio.SeqIO as SeqIO
import numpy as np
import pysam
import scipy.sparse as scisp
import tqdm
import os
from utils import count_fasta_sequences, open_input
import logging

# package logger
logger = logging.getLogger(__name__)

#######offset is the enumeration of the length######################
#localid is the index of the list#
#refid is the index of the fasta file, which is a global index# 
SeqInfo = namedtuple('SeqInfo', ['localid', 'refid', 'name', 'sites', 'length', 'covcc']) #create a new class of tuple: SeqInfo



class SiteCounter(object):

    def __init__(self, enzyme_names,  is_linear=True):
        """
        Simple class to count the total number of enzymatic cut sites for the given
        list if enzymes.
        :param enzyme_names: a list of enzyme names (proper case sensitive spelling a la NEB)
        :param tip_size: when using tip based counting, the size in bp
        :param is_linear: Treat sequence as linear.
        """
        if isinstance(enzyme_names, str):
            enzyme_names = [enzyme_names]
        assert (isinstance(enzyme_names, Iterable) and
                not isinstance(enzyme_names, str)), 'enzyme_names must of a collection of names'
        self.enzymes = [getattr(Restriction , en) for en in enzyme_names] ##get the enzyme name from the standard module
        self.is_linear = is_linear


    def count_sites(self, seq):
        """
        Count the number of sites found in the given sequence, where sites from
        all specified enzymes are combined

        :param seq: Bio.Seq object
        :return: the total number of sites
        """

        return sum(len(en.search(seq, self.is_linear)) for en in self.enzymes)




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

    def __init__(self, bam_file, enzymes, seq_file, path , min_mapq=30, min_len=1000, min_match=30, min_signal=1):

        ########input the parameter################
        ########################################################################
        ########################################################################
        #bam_file: alignment info of Hi-C library on contigs in bam#
        #enzymes: name of restriction enzymes used in Hi-C experiments#
        #seq_file: store the assembly contigs in fasta#
        #path: output path#
        #min_mapq: minimal mapping quality(default 30)#
        #min_len: minimal length of contigs(default 1000bp)

        self.bam_file = bam_file
        self.enzymes = enzymes
        self.seq_file = seq_file
        self.path = path
        self.min_mapq = min_mapq
        self.min_len = min_len
        self.min_match = min_match
        self.min_signal = min_signal
        #fasta_info store the info from fasta file#
        #seq_info store the information of contigs from bam file#
        #seq_map store the contact map#
        self.fasta_info = {}
        self.seq_info = []
        self.seq_map = None
        self.total_reads = None

        
        logger.info('Reading fasta file...')
        with open_input(seq_file) as multi_fasta:
            # prepare the site counter for the given experimental conditions
            # fasta_info is a dictionary of preparation of seq_info and seq_info is the true results
            site_counter = SiteCounter(enzymes,  is_linear=True)
            # get an estimate of sequences for progress
            fasta_count = count_fasta_sequences(seq_file)
            for seqrec in tqdm.tqdm(SeqIO.parse(multi_fasta, 'fasta'), total=fasta_count , desc='Analyzing contigs in reference fasta'):
                if len(seqrec) < min_len:
                    continue
                self.fasta_info[seqrec.id] = {'sites': site_counter.count_sites(seqrec.seq),
                                         'length': len(seqrec)}

        logger.debug('There are {} contigs in reference fasta'.format(fasta_count))

        
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
                    continue

                assert fa['length'] == rlen, 'Sequence lengths in {} do not agree: bam {} fasta {}'.format(rname, fa['length'], rlen)

                self.seq_info.append(SeqInfo(localid , n , rname, fa['sites'], rlen, 0)) ######initially set the covcc coverage to be zero
                localid = localid + 1
                offset += rlen

            ####### total length of contigs##########
            ####### total_seq is number of contigs###
            self.total_len = offset
            self.total_seq = localid

            del self.fasta_info

            if self.total_seq == 0:
                raise ImportError('No sequences in BAM file can be found in FASTA file')
            
            logger.debug('{} contigs miss and {} contigs are too short'.format(ref_count['seq_missing'] , ref_count['too_short']))
            logger.debug('Accepted {} contigs covering {} bp'.format(self.total_seq, self.total_len))
            logger.info('Counting reads in bam file...')

 
            self.total_reads = bam.count(until_eof=True)
            logger.debug('BAM file contains {0} alignments'.format(self.total_reads))

            logger.info('Handling the alignments...')
            self._bin_map(bam)
                

        logger.info('Filtering contigs according to minimal signal({})...'.format(self.min_signal))
        contig_id = self.max_offdiag()
        logger.debug('{} contigs remain'.format(len(contig_id)))
        
        self.seq_map = self.seq_map.tolil()
        seq_temp = [] ###temporately store the sequence information#######
        for i , idn in enumerate(contig_id):
            seq = self.seq_info[idn]
            assert seq.localid == idn, 'the local index does not match the contact matrix index'
            seq_temp.append(SeqInfo(i , seq.refid , seq.name , seq.sites , seq.length , self.seq_map[idn , idn]))

        self.seq_info = seq_temp
        del seq_temp

        self.seq_map = self.seq_map.tocsr()
        self.seq_map = self.seq_map[contig_id , :]
        self.seq_map = self.seq_map.tocsc()
        self.seq_map = self.seq_map[: , contig_id]
        self.seq_map = self.seq_map.tocoo()
        del contig_id
        
        assert self.seq_map.shape[0] == len(self.seq_info), 'Filter error'
        ######Ouput raw Hi-C contact matrix##########
        scisp.save_npz(os.path.join(self.path, 'Raw_contact_matrix.npz'), self.seq_map.tocsr())
        
        ########change the diaganol entries of Hi-C matrix to zero#######
        self.seq_map = self.seq_map.tolil()
        self.seq_map.setdiag(0)
        
        #######Compute the Hi-C signals for each contig########
        self.seq_map = self.seq_map.tocsr()
        self.row_sum = np.matrix.tolist(self.seq_map.sum(axis=0))[0]
        self._write_contig_info()
        
    def _bin_map(self, bam):
        """
        Accumulate read-pair observations from the supplied BAM file.
        Maps are initialized here. Logical control is achieved through initialisation of the
        ContactMap instance, rather than supplying this function arguments.

        :param bam: this instance's open bam file.
        """

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
                'accepted map_different_contig pairs': 0,
                'accepted map_same_contig pairs': 0,
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

                if not _matcher(r1) or not _matcher(r2):
                    counts['poor_match pairs'] += 1
                    continue

                if r1.reference_id == r2.reference_id:
                    counts['accepted map_same_contig pairs'] += 1
                else:
                    counts['accepted map_different_contig pairs'] += 1
                    
                # get internal indices
                ix1 = _idx[r1.reference_id]
                ix2 = _idx[r2.reference_id]

                # maintain just a half-matrix
                if ix2 < ix1:
                    ix1, ix2 = ix2, ix1

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


    def _write_contig_info(self):
        with open(os.path.join(self.path , 'tmp' , 'contig_info.csv'),'w') as out:
            for i , seq in enumerate(self.seq_info):
                out.write(str(seq.name)+ ',' +str(seq.sites)+ ',' +str(seq.length) + ',' +str(seq.covcc) + ',' + str(self.row_sum[i]))
                out.write('\n')
                
        with open(os.path.join(self.path , 'contig_info.csv'),'w') as out:
            out.write(str('Contig name')+ ','  + str('Number of restriction sites')+ ',' +str('Contig length') + '\n')
            for seq in self.seq_info:
                out.write(str(seq.name)+ ',' + str(seq.sites)+ ',' + str(seq.length))
                out.write('\n')
        
        
    def max_offdiag(self):
        """
        Determine the maximum off-diagonal values of a given symmetric matrix. As this
        is assumed to be symmetric, we consider only the rows.

        :param _m: a scipy.sparse matrix
        :return: the off-diagonal maximum values
        """
        _m = self.seq_map
        assert scisp.isspmatrix(_m), 'Input matrix is not a scipy.sparse object'
        _m = _m.tolil(True)
        _m.setdiag(0)
        _sig = np.asarray(_m.tocsr().max(axis=0).todense()).ravel()
        _contig_id = []
        for i in range(_m.shape[0]):
            if _sig[i] >= self.min_signal:
                _contig_id.append(i)
        del _m
        return _contig_id
