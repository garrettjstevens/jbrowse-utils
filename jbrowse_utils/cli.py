import argparse
import json

from .prepare_refseqs import format_sequences


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(
        dest='COMMAND',
        # help='Choose one:',
        title='Commands',
        description='These are the available commands:')
    subparsers.required = True
    parser_prepare_refseqs = subparsers.add_parser(
        'prepare-refseqs',
        help='Formats reference sequence data for use with JBrowse',
        description='Formats reference sequence data for use with JBrowse')
    configure_prepare_refseqs_parser(parser_prepare_refseqs)
    args = parser.parse_args()
    args.func(args)


def configure_prepare_refseqs_parser(parser):
    parser.set_defaults(func=parse_prepare_refseqs)
    input_group = parser.add_argument_group(
        title='SEQUENCE FILE',
        description='Specify one of:')
    group = input_group.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '--gff',
        metavar='<GFF3 file>',
        help='Must be GFF version 3 with an embedded FASTA section')
    group.add_argument(
        '--fasta',
        dest='fastas',
        nargs='+',
        action='append',
        metavar='<FASTA file>',
        help=(
            'Can specify multiple FASTAs after one flag or use the flag '
            'multiple times. Can be a gzipped file (ending in .gz or .gzip). '
            'Can optionally supply --refs.'))
    group.add_argument(
        '--indexed_fasta',
        metavar='<FASTA file>',
        help='An index with the same name plus ".fai" must be present.')
    group.add_argument(
        '--twobit',
        metavar='<2BIT file>',
        help='A single .2bit file.')
    group.add_argument(
        '--conf',
        metavar='<biodb config file>',
        help=(
            'JBrowse config file pointing to a BioPerl database. Must also '
            'supply --refs.'))
    group.add_argument(
        '--sizes',
        nargs='+',
        action='append',
        metavar='<sizes file>',
        help=(
            'Can specify multiple sizes files after one flag or use the flag '
            'multiple times.'))
    group.add_argument(
        '--gff-sizes',
        nargs='+',
        action='append',
        metavar='<GFF file>',
        help=(
            'Must contain ##sequence-region lines as described in the '
            'GFF specs. Can specify multiple GFF files after one flag or use '
            'the flag multiple times.'))

    options_group = parser.add_argument_group(
        title='OPTIONS')
    options_group.add_argument(
        '--noSort',
        action='store_false',
        dest='sort',
        help=(
            'If using GFF or FASTA input, preserve the order of the reference '
            'sequences (sorts alphabetically by default).'))
    options_group.add_argument(
        '--out',
        metavar='<output directory>',
        default='data/',
        help='Optional directory to write to. Defaults to data/.')
    options_group.add_argument(
        '--noseq',
        action='store_false',
        dest='seq',
        help=(
            'Do not store the actual sequence bases, just the sequence metadata '
            '(name, length, and so forth).'))
    options_group.add_argument(
        '--refs',
        nargs='+',
        metavar='<refseq names>',
        help=(
            'Output only the sequences with the given (comma- or space-'
            'separated) names.'))
    options_group.add_argument(
        '--compress',
        action='store_true',
        help=(
            'If passed, compress the reference sequences with gzip, making the '
            'chunks be .txt.gz. NOTE: this requires a bit of additional web '
            'server configuration to be served correctly.'))
    options_group.add_argument(
        '--chunksize',
        type=int,
        metavar='<int>',
        default=20000,
        help=(
            'Size of sequence chunks to make, in base pairs. Default 20kb. '
            'This is multiplied by 4 if --compress is passed, so that the '
            'compressed sequence files are still approximately this size.'))
    options_group.add_argument(
        '--nohash',
        action='store_false',
        dest='hash',
        help=(
            'Store sequences in a flat seq/{seqname}/{chunk}.txt structure, '
            'instead of the new (more scalable) '
            '/seq/hash/hash/hash/{seqname}-{chunk}.txt structure.'))
    options_group.add_argument(
        '--trackLabel',
        metavar='<label>',
        help=(
            'The unique name of the sequence track, '
            'default is to set to the same as --seqType.'))
    options_group.add_argument(
        '--key',
        metavar='<string>',
        default='Reference sequence',
        help=(
            'The displayed name of the sequence track, defaults to '
            '"Reference sequence".'))
    options_group.add_argument(
        '--seqType',
        metavar='<string>',
        default='dna',
        help=('The name of the alphabet used for these reference sequences, '
              'usually either "dna", "rna", or "protein". (default "dna")'))
    options_group.add_argument(
        '--trackConfig',
        type=json.loads,
        metavar='{ JSON-format extra configuration for this track }',
        help=(
            'Additional top-level configuration for the client, in JSON syntax. '
            'Example: --trackConfig \'{ "glyph": "ProcessedTranscript" }\''))


def parse_prepare_refseqs(args):
    # Turn e.g. [['a.fa', 'b.fa'], ['c.fa']] into ['a.fa', 'b.fa', 'c.fa']
    # https://stackoverflow.com/questions/952914/making-a-flat-list-out-of-list-of-lists-in-python
    delattr(args, 'COMMAND')
    delattr(args, 'func')
    if args.fastas:
        args.fastas = [fasta for arg in args.fastas for fasta in arg]
    if args.sizes:
        args.sizes = [size for arg in args.sizes for size in arg]
    if args.gff_sizes:
        args.gff_sizes = [size for arg in args.gff_sizes for size in arg]
    if args.refs:
        args.refs = [ref for comma_list in args.refs for ref in comma_list.split(',')]
    format_sequences(**vars(args))
