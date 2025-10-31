#########The structure of the main script is modified from bin3C########
from raw_contact import RawContact
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
    cm_raw = subparsers.add_parser('raw', parents=[global_parser],
                                      description='Raw contacts.')

    cmd_pl = subparsers.add_parser('pipeline', parents=[global_parser],
                                      description='Retrieve complete viral genomes.')


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




    except ApplicationException:
        logger.error('ApplicationException Error')
        sys.exit(1)

