#########The structure of the main script is modified from bin3C########
from raw_contact import ContactMatrix as RawContact
from construct_graph import ContactMatrix
from bin import ClusterBin
from exceptions import ApplicationException
from utils import make_dir, gen_bins
import logging
import sys
import argparse
import os
import time
import warnings

##Ignore the warning information of package deprecation##
warnings.filterwarnings("ignore")

__version__ = '1.0.0, released at 03/2022'

if __name__ == '__main__':
    
    def mk_version():
        return 'ViralCC v{}'.format(__version__)

    def out_name(base, suffix):
        return '{}{}'.format(base, suffix)

    def ifelse(arg, default):
        if arg is None:
            return default
        else:
            return arg

    runtime_defaults = {
        'min_len': 1000,
        'min_signal': 1,
        'min_mapq': 30,
        'min_match': 30,
        'min_k': 4,
        'random_seed': 42
    }

    global_parser = argparse.ArgumentParser(add_help=False)
    global_parser.add_argument('-V', '--version', default=False, action='store_true', help='Show the application version')
    global_parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    global_parser.add_argument('--cover', default=False, action='store_true', help='Cover existing files')
    global_parser.add_argument('--log', help='Log file path [OUTDIR/viralcc.log]')


    parser = argparse.ArgumentParser(description='ViralCC: a metagenomic proximity-based tool to retrieve complete viral genomes')

    subparsers = parser.add_subparsers(title='commands', dest='command', description='Valid commands',
                                       help='choose an analysis stage for further options')
    cmd_raw = subparsers.add_parser('raw', parents=[global_parser],
                                      description='Raw contacts.')

    cmd_pl = subparsers.add_parser('pipeline', parents=[global_parser],
                                      description='Retrieve complete viral genomes.')
    
    cmd_link = subparsers.add_parser('link', parents=[global_parser],
                                      description='Find MGE-host linkages from contact matrix.')

    '''
    Generating raw contacts subparser input
    '''
    cmd_raw.add_argument('--min-len', type=int,
                           help='Minimum acceptable contig length [1000]')
    cmd_raw.add_argument('--min-signal', type=int,
                           help='Minimum acceptable Hi-C signal [1]')
    cmd_raw.add_argument('--min-mapq', type=int,
                           help='Minimum acceptable mapping quality [30]')
    cmd_raw.add_argument('--min-match', type=int,
                           help='Accepted alignments must being N matches [30]')
    cmd_raw.add_argument('-e', '--enzyme', metavar='NEB_NAME', action='append',
                           help='Case-sensitive enzyme name. Use multiple times for multiple enzymes')
    cmd_raw.add_argument('FASTA', help='Reference fasta sequence')
    cmd_raw.add_argument('BAM', help='Input bam file in query order')
    cmd_raw.add_argument('OUTDIR', help='Output directory')


    '''
    pipeline subparser input
    '''
    cmd_pl.add_argument('--min-len', type=int,
                               help='Minimum acceptable reference length [3000]')
    cmd_pl.add_argument('--min-mapq', type=int,
                           help='Minimum acceptable mapping quality [30]')
    cmd_pl.add_argument('--min-match', type=int,
                           help='Accepted alignments must being N matches [30]')
    cmd_pl.add_argument('--min-k', type=int,
                           help='Lower bound of k for determining host poximity graph [4]')
    cmd_pl.add_argument('--random-seed', type=int,
                           help='Random seed to run the Leiden algorithm [42]')
    cmd_pl.add_argument('FASTA', help='Reference fasta sequence')
    cmd_pl.add_argument('BAM', help='Input bam file in query order')
    cmd_pl.add_argument('VIRAL', help='Viral contig labels')
    cmd_pl.add_argument('OUTDIR', help='Output directory')
    
    '''
    link subparser input
    '''
    cmd_link.add_argument('--viral-bin', required=True,
                          help='Viral bin folder (ViralCC output)')
    cmd_link.add_argument('--host-bin', required=True,
                          help='Host bin folder')
    cmd_link.add_argument('--contact', required=True,
                          help='Contact matrix file (npz format)')
    cmd_link.add_argument('--contig-info', required=True,
                          help='Contig information file')
    cmd_link.add_argument('--viral-annotation', default=None,
                          help='Optional: Viral bin annotation (VIRGO output)')
    cmd_link.add_argument('--host-annotation', default=None,
                          help='Optional: Host bin annotation (GTDBTK output)')
    cmd_link.add_argument('OUTDIR', help='Output directory')

    args = parser.parse_args()

    if args.version:
        print(mk_version())
        sys.exit(0)

    try:
        make_dir(args.OUTDIR, args.cover)
    except IOError:
        print('Error: cannot find out directory or the directory already exists')
        sys.exit(1)

    logging.captureWarnings(True)
    logger = logging.getLogger('main')

    # root log listens to everything
    root = logging.getLogger('')
    root.setLevel(logging.DEBUG)

    # log message format
    formatter = logging.Formatter(fmt='%(levelname)-8s | %(asctime)s | %(name)7s | %(message)s')

    # Runtime console listens to INFO by default
    ch = logging.StreamHandler()
    if args.verbose:
        ch.setLevel(logging.DEBUG)
    else:
        ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    root.addHandler(ch)

    # File log listens to all levels from root
    if args.log is not None:
        log_path = args.log
    else:
        log_path = os.path.join(args.OUTDIR, 'viralcc.log')
    fh = logging.FileHandler(log_path, mode='a')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    root.addHandler(fh)

    # Add some environmental details
    logger.debug(mk_version())
    logger.debug(sys.version.replace('\n', ' '))
    logger.debug('Command line: {}'.format(' '.join(sys.argv)))

    try:
        if args.command == 'raw':
            if args.enzyme is not None:
                logger.info('Begin constructing raw contact matrix...')
                cm = RawContact(args.BAM,
                                args.enzyme,
                                args.FASTA,
                                args.OUTDIR,
                                min_mapq=ifelse(args.min_mapq, runtime_defaults['min_mapq']),
                                min_len=ifelse(args.min_len, runtime_defaults['min_len']),
                                min_match=ifelse(args.min_match, runtime_defaults['min_match']),
                                min_signal=ifelse(args.min_signal, runtime_defaults['min_signal']))

                logger.info('Raw contact matrix construction finished')
                
        if args.command == 'pipeline':
            start_time = time.time()
            cm = ContactMatrix(args.BAM,
                                args.FASTA,
                                args.VIRAL,
                                args.OUTDIR,
                                min_mapq=ifelse(args.min_mapq, runtime_defaults['min_mapq']),
                                min_len=ifelse(args.min_len, runtime_defaults['min_len']),
                                min_match=ifelse(args.min_match, runtime_defaults['min_match']),
                                min_k=ifelse(args.min_k, runtime_defaults['min_k']))


            cl = ClusterBin(args.OUTDIR,
                            cm.viral_info ,
                            cm.map_combine ,
                            random_seed = ifelse(args.random_seed, runtime_defaults['random_seed']))
                            
            logger.info('Clustering fininshed')
            logger.info('Writing bins...')
            gen_bins(args.FASTA , os.path.join(args.OUTDIR ,'cluster_viral_contig.txt') , os.path.join(args.OUTDIR ,'VIRAL_BIN'))
            
            end_time = time.time()
            logger.info('ViralCC consumes {} seconds in total'.format(str(end_time-start_time)))

        if args.command == 'link':
            logger.info('Begin finding MGE-host linkages...')
            
            import pandas as pd
            import numpy as np
            from scipy.sparse import load_npz
            from Bio import SeqIO
            import glob
            
            # Load contact matrix
            logger.info('Loading contact matrix...')
            contact_matrix = load_npz(args.contact).tocoo()
            
            # Read contig information to map indices to contig names
            logger.info('Reading contig information...')
            contig_names = []
            with open(args.contig_info, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        contig_names.append(line[1:].split()[0])
            
            # Get viral contigs from viral bin folder
            logger.info('Identifying viral contigs...')
            viral_contigs = set()
            viral_bin_files = glob.glob(os.path.join(args.viral_bin, '*.fa')) + \
                             glob.glob(os.path.join(args.viral_bin, '*.fasta'))
            for bin_file in viral_bin_files:
                for record in SeqIO.parse(bin_file, 'fasta'):
                    viral_contigs.add(record.id)
            
            logger.info('Found {} viral contigs'.format(len(viral_contigs)))
            
            # Get host contigs from host bin folder
            logger.info('Identifying host contigs...')
            host_contigs = set()
            host_bin_files = glob.glob(os.path.join(args.host_bin, '*.fa')) + \
                            glob.glob(os.path.join(args.host_bin, '*.fasta'))
            for bin_file in host_bin_files:
                for record in SeqIO.parse(bin_file, 'fasta'):
                    host_contigs.add(record.id)
            
            logger.info('Found {} host contigs'.format(len(host_contigs)))
            
            # Load optional annotations
            viral_annotations = {}
            host_annotations = {}
            
            if args.viral_annotation is not None:
                logger.info('Loading viral annotations...')
                # Assuming VIRGO output format - adjust as needed
                try:
                    virgo_df = pd.read_csv(args.viral_annotation, sep='\t')
                    # Adjust column names based on actual VIRGO output format
                    if 'contig_id' in virgo_df.columns and 'annotation' in virgo_df.columns:
                        viral_annotations = dict(zip(virgo_df['contig_id'], virgo_df['annotation']))
                    logger.info('Loaded {} viral annotations'.format(len(viral_annotations)))
                except Exception as e:
                    logger.warning('Failed to load viral annotations: {}'.format(str(e)))
            
            if args.host_annotation is not None:
                logger.info('Loading host annotations...')
                # Assuming GTDBTK output format - adjust as needed
                try:
                    gtdbtk_df = pd.read_csv(args.host_annotation, sep='\t')
                    # Adjust column names based on actual GTDBTK output format
                    if 'user_genome' in gtdbtk_df.columns and 'classification' in gtdbtk_df.columns:
                        host_annotations = dict(zip(gtdbtk_df['user_genome'], gtdbtk_df['classification']))
                    logger.info('Loaded {} host annotations'.format(len(host_annotations)))
                except Exception as e:
                    logger.warning('Failed to load host annotations: {}'.format(str(e)))
            
            # Extract virus-host linkages
            logger.info('Extracting virus-host linkages...')
            linkages = []
            
            for r, c, v in zip(contact_matrix.row, contact_matrix.col, contact_matrix.data):
                if v <= 0:
                    continue
                
                contig_a = contig_names[r]
                contig_b = contig_names[c]
                
                # Check if one is viral and one is host
                if contig_a in viral_contigs and contig_b in host_contigs:
                    linkage = {
                        'viral_contig': contig_a,
                        'host_contig': contig_b,
                        'contact_strength': int(v)
                    }
                    if viral_annotations:
                        linkage['viral_annotation'] = viral_annotations.get(contig_a, '')
                    if host_annotations:
                        linkage['host_annotation'] = host_annotations.get(contig_b, '')
                    linkages.append(linkage)
                    
                elif contig_b in viral_contigs and contig_a in host_contigs:
                    linkage = {
                        'viral_contig': contig_b,
                        'host_contig': contig_a,
                        'contact_strength': int(v)
                    }
                    if viral_annotations:
                        linkage['viral_annotation'] = viral_annotations.get(contig_b, '')
                    if host_annotations:
                        linkage['host_annotation'] = host_annotations.get(contig_a, '')
                    linkages.append(linkage)
            
            # Create DataFrame and save
            logger.info('Creating output table...')
            linkages_df = pd.DataFrame(linkages)
            
            if len(linkages_df) > 0:
                # Sort by contact strength
                linkages_df = linkages_df.sort_values('contact_strength', ascending=False)
                
                # Save to output
                output_file = os.path.join(args.OUTDIR, 'virus_host_linkages.tsv')
                linkages_df.to_csv(output_file, sep='\t', index=False)
                
                logger.info('Found {} virus-host linkages'.format(len(linkages_df)))
                logger.info('Unique viral contigs with hosts: {}'.format(linkages_df['viral_contig'].nunique()))
                logger.info('Unique host contigs with viruses: {}'.format(linkages_df['host_contig'].nunique()))
                logger.info('Output saved to: {}'.format(output_file))
            else:
                logger.warning('No virus-host linkages found')
            
            logger.info('MGE-host linkage analysis completed')


    except ApplicationException:
        logger.error('ApplicationException Error')
        sys.exit(1)