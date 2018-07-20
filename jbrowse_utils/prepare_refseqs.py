import binascii
import gzip
import json
import os
import os.path
import re
import shutil
import struct

from .genome_db import GenomeDB
from .json_file_store import JsonFileStore


def format_sequences(
        # *,
        gff=None,
        fastas=None,
        indexed_fasta=None,
        twobit=None,
        conf=None,
        sizes=None,
        gff_sizes=None,

        sort=True,
        out='data/',
        seq=True,
        refs=None,
        compress=False,
        chunksize=20000,
        hash=True,
        trackLabel=None,
        key='Reference sequence',
        seqType='dna',
        trackConfig=None):
    """ Format reference sequence data for use with JBrowse

    **While documented, `conf` is currently not implemented (it's a work in
    progress) and is use will result in a NotImplementedError.**

    Only one of the input types should be specified:
    :param gff: Must be GFF version 3 with an embedded FASTA section
    :type fastas: list
    :param fastas: Can be a gzipped file (ending in .gz or .gzip).
    :param indexed_fasta: An index with the same name plus ".fai" must be
        present.
    :param twobit: A single .2bit file.
    :param conf: JBrowse config file pointing to a BioPerl database.
        Must also supply --refs.
    :type sizes: list
    :param sizes: Two-column "ref_name ref_size" file
    :type gff_sizes: list
    :param gff_sizes: Must contain ##sequence-region lines as described in the
        GFF specs.

    Other optional parameters
    :param sort: If using GFF or FASTA input, sorts reference sequences
        alphabetically (True by default). If false, preserve the order of the
        reference sequences.
    :param out: Optional directory to write to. Defaults to data/.
    :param seq: Whether to store the actual sequence bases, as opposed to just
        the sequence metadata (name, length, etc.). Defaults to True.
    :type refs: list
    :param refs: Output only the sequences with the given names.
    :param compress: Compress the reference sequences with gzip, making the
        chunks be .txt.gz. Defaults to False. NOTE: this requires a bit of
        additional web server configuration to be served correctly.
    :param chunksize: Size of sequence chunks to make, in base pairs.
        Default is 20,000. This is multiplied by 4 if compress is True, so that
        the compressed sequence files are still approximately this size.
    :param hash: Store sequences in a /seq/hash/hash/hash/{seqname}-{chunk}.txt
        structure. Defaults to True. If False, store in the old (less
        scalable) seq/{seqname}/{chunk}.txt structure.
    :param trackLabel: The unique name of the sequence track, default is to set
        to the same as seqType.
    :param key: The displayed name of the sequence track, defaults to
        "Reference sequence".
    :param seqType: The name of the alphabet used for these reference sequences,
        usually either "dna", "rna", or "protein". Defaults to "dna"
    :param trackConfig: dict
    :param trackConfig: Additional top-level configuration for the client,
        in JSON-serializable dict.
            Example: {'glyph': 'ProcessedTranscript'}
    """
    seq_types = [gff, fastas, indexed_fasta, twobit, conf, sizes, gff_sizes]
    if not sum(bool(seq_type) for seq_type in seq_types) == 1:
        raise ValueError(
            'Must specify one (and only one) of the following sequence types: '
            '{}'.format(', '.join(seq_types)))
    if conf:
        raise NotImplementedError(
            'The "--conf" option has not yet been implemented')
    if compress:
        chunksize = 4 * chunksize
    if trackConfig:
        if isinstance(trackConfig, str):
            trackConfig = json.loads(trackConfig)

    json_store = JsonFileStore(out, compress)

    opts = {
        'sort': sort,
        'out': out,
        'seq': seq,
        'refs': refs,
        'compress': compress,
        'chunk_size': chunksize,
        'hash': hash,
        'track_label': trackLabel,
        'key': key,
        'seq_type': seqType,
        'track_config': trackConfig}

    if indexed_fasta:
        export_fai(indexed_fasta, json_store, **opts)
        write_track_entry('indexed_fasta', json_store, src=indexed_fasta, **opts)
    elif twobit:
        export_twobit(twobit, json_store, **opts)
        write_track_entry('twobit', json_store, src=twobit, **opts)
    elif fastas:
        export_fastas(fastas, json_store, **opts)
        write_track_entry('fastas', json_store, **opts)
    elif gff:
        export_gff(gff, json_store, **opts)
        write_track_entry('gff', json_store, **opts)
    elif conf:
        export_conf(conf, json_store, **opts)
        write_track_entry('conf', json_store, **opts)
    elif sizes:
        export_sizes(sizes, json_store, **opts)
        write_track_entry('sizes', json_store, **opts)
    elif gff_sizes:
        export_gff_sizes(gff_sizes, json_store, **opts)
        write_track_entry('gff_sizes', json_store, **opts)


def export_fai(indexed_fasta, json_store, **opts):
    """opts are same as optional parameters for format_sequences()"""
    refseqs = {}
    original_order = []
    fai = indexed_fasta + '.fai'
    with open(fai) as infile:
        for line in infile:
            fai_match = re.match(
                r'^([^\t]+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)$', line.strip())
            if fai_match:
                refseqs[fai_match.group(1)] = {
                    'name': fai_match.group(1),
                    'start': 0,
                    'end': fai_match.group(2),
                    'offset': fai_match.group(3),
                    'line_length': fai_match.group(4),
                    'line_byte_length': fai_match.group(5)}
                if not opts['sort']:
                    original_order.append(fai_match.group(1))
            else:
                raise ValueError(
                    'Improperly-formatted line in fai file ({}):\n{}'.format(
                        fai, line.strip()))
    direc = os.path.join(opts['out'], 'seq')
    os.makedirs(direc, exist_ok=True)
    shutil.copy(fai, direc)
    shutil.copy(indexed_fasta, direc)
    write_refseqs_json(refseqs, json_store, original_order, **opts)


def export_twobit(twobit, json_store, **opts):
    """opts are same as optional parameters for format_sequences()"""
    toc = {}
    refseqs = {}
    with open(twobit, 'rb') as infile:
        raw = infile.read(4 * 4)
        for template in ('<', '>'):
            sig, ver, cnt, reserved = struct.unpack(template + '4L', raw)
            if sig == 0x1A412743:
                break
        else:
            raise ValueError('Invalid 2bit file: {}'.format(twobit))

        for i in range(cnt):
            raw = infile.read(1)
            # Read size of record name
            size, = struct.unpack('B', raw)
            # Read name of the record
            name = infile.read(size)
            # Read and store offset
            raw = infile.read(4)
            toc[name], = struct.unpack(template + 'L', raw)

        for name in sorted(toc):
            offset = toc[name]
            infile.seek(offset)
            raw = infile.read(4)
            size, = struct.unpack(template + 'L', raw)
            refseqs['name'] = {
                'name': name,
                'length': size,
                'start': 0,
                'end': size}
    direc = os.path.join(opts['out'], 'seq')
    os.makedirs(direc, exist_ok=True)
    shutil.copy(twobit, direc)
    write_refseqs_json(refseqs, json_store, **opts)


def export_fastas(fastas, json_store, **opts):
    """opts are same as optional parameters for format_sequences()"""
    accept_all_refs = False if opts['refs'] else True
    refseqs = {}
    original_order = []
    for fasta in fastas:
        if hasattr(fasta, 'read'):
            infile = fasta
        elif fasta.endswith('.gz') or fasta.endswith('.gzip'):
            infile = gzip.open(fasta, 'rb')
        else:
            infile = open(fasta, 'rb')

        curr_seq = {}
        curr_chunk = b''
        chunk_num = 0

        fasta_header = re.compile(rb'^\s*>\s*(\S+)\s*(.*)')
        for line in infile:
            header_match = re.match(fasta_header, line)
            if header_match:
                if curr_seq:
                    curr_seq, curr_chunk, chunk_num = _write_fasta_chunks(
                        curr_seq, curr_chunk, chunk_num, flush=True, **opts)
                if accept_all_refs or header_match.group(1) in opts['refs']:
                    chunk_num = 0
                    curr_chunk = b''
                    curr_seq = refseqs[header_match.group(1).decode()] = {
                        'name': header_match.group(1),
                        'start': 0,
                        'end': 0,
                        'seqChunkSize': opts['chunk_size']}
                    if header_match.group(2):
                        curr_seq['description'] = header_match.group(2)
                    if not opts['sort']:
                        original_order.append(header_match.group(1))
                else:
                    curr_seq = {}
            elif curr_seq and line.strip():
                line = re.sub(rb'[\s\r\n]', b'', line)
                curr_seq['end'] += len(line)
                if opts['seq']:
                    curr_chunk += line
                    curr_seq, curr_chunk, chunk_num = _write_fasta_chunks(
                        curr_seq, curr_chunk, chunk_num, **opts)
        _write_fasta_chunks(curr_seq, curr_chunk, chunk_num, flush=True, **opts)
        infile.close()
        write_refseqs_json(refseqs, json_store, original_order, **opts)


def export_gff(gff, json_store, **opts):
    """opts are same as optional parameters for format_sequences()"""
    with open(gff, 'rb') as infile:
        for line in infile:
            if re.match(rb'^##FASTA\s*$', line, flags=re.IGNORECASE):
                # start of the sequence block, pass the filehandle to our fasta database
                export_fastas([infile], json_store, **opts)
                break
            elif re.match(rb'^>', line):
                # beginning of implicit sequence block, need to seek back
                infile.seek(-len(line), 1)
                export_fastas([infile], json_store, **opts)
                break


def export_conf(conf, json_store, **opts):
    """ CURRENTLY NOT IMPLEMENTED

    opts are same as optional parameters for format_sequences()"""
    pass


def export_sizes(sizes, json_store, **opts):
    """opts are same as optional parameters for format_sequences()"""
    refseqs = {}
    for sizefile in sizes:
        with open(sizefile) as infile:
            for line in infile:
                if line.strip():
                    name, length = line.strip().split()
                    refseqs[name] = {
                        'name': name,
                        'start': 0,
                        'end': int(length),
                        'length': int(length)}
    write_refseqs_json(refseqs, json_store, **opts)


def export_gff_sizes(gff_sizes, json_store, **opts):
    """opts are same as optional parameters for format_sequences()"""
    refseqs = {}
    for sizefile in gff_sizes:
        with open(sizefile) as infile:
            for line in infile:
                if line.startswith('^##sequence-region'):
                    _, name, start, end = line.strip().split()
                    refseqs[name] = {
                        'name': name,
                        'start': int(start) - 1,
                        'end': int(end),
                        'length': int(end)}
    write_refseqs_json(refseqs, json_store, **opts)


def _write_fasta_chunks(curr_seq, curr_chunk, chunk_num, flush=False, **opts):
    if not opts['seq']:
        return
    while (flush and curr_chunk) or len(curr_chunk) >= opts['chunk_size']:
        shift = curr_chunk[:opts['chunk_size']]
        curr_chunk = curr_chunk[opts['chunk_size']:]
        with _open_chunk_file(curr_seq, chunk_num, **opts) as outfile:
            outfile.write(shift)
        chunk_num += 1
    return curr_seq, curr_chunk, chunk_num


def _open_chunk_file(ref_info, chunk_num, **opts):
    if opts['hash']:
        direc = os.path.join(opts['out'], 'seq', *_crc32_path(ref_info['name']))
        file = os.path.join(
            direc,
            '{}-{}.txt'.format(str(ref_info['name'].decode()), chunk_num))
    else:
        direc = os.path.join(opts['out'], 'seq', str(ref_info['name'].decode()))
        file = os.path.join(
            direc,
            '{}.txt'.format(chunk_num))
    os.makedirs(direc, exist_ok=True)
    if opts['compress']:
        file += 'z'
        outfile = gzip.open(file, 'wb')
    else:
        outfile = open(file, 'wb')
    return outfile


def _crc32_path(string):
    crc = (binascii.crc32(string) & 0xffffffff)
    crc_hex = '{:08x}'.format(crc)
    # e.g. '6092e02d' becomes ['609', '2e0', '2d']
    return [crc_hex[i:i + 3] for i in range(0, len(crc_hex), 3)]


def write_track_entry(seq_source, json_store, src=None, **opts):
    """ Write a a track entry for the prepared reference sequence

    :param seq_source: Type of sequence from which reference was prepared
    :param json_store: JsonFileStore object defining the output directory
    :param src: For "indexed_fasta" and "twobit", the name of the input file
    :param opts: Same as optional parameters for format_sequences()
    """
    if not opts['seq']:
        return
    if opts['track_label']:
        seq_track_name = opts['track_label']
    elif re.match(r'^[dr]na$', opts['seq_type'], flags=re.IGNORECASE):
        seq_track_name = opts['seq_type'].upper()
    else:
        seq_track_name = opts['seq_type'].lower()
    json_store.touch('tracks.conf')

    def add_track(data):
        if opts['hash']:
            seq_url_template = 'seq/{refseq_dirpath}/{refseq}-'
        else:
            seq_url_template = 'seq/{refseq}/'
        track = {
            'label': seq_track_name,
            'key': opts['key'],
            'type': 'SequenceTrack',
            'category': 'Reference sequence',
            'storeClass': 'JBrowse/Store/Sequence/StaticChunked',
            'chunkSize': opts['chunk_size'],
            'urlTemplate': seq_url_template,
            'seqType': opts['seq_type'].lower()
        }
        if opts['compress']:
            track['compress'] = 1
        if not opts['seq_type'].lower() == 'dna':
            track['showReverseStrand'] = 0
        if opts['seq_type'].lower() == 'protein':
            track['showTranslation'] = 0
        # Merge in any extra trackConfig supplied by the user.
        if opts['track_config']:
            track.update(opts['track_config'])
        if seq_source == 'indexed_fasta':
            del track['chunkSize']
            fasta_template = os.path.join('seq', os.path.basename(src))
            track.update({
                'storeClass': 'JBrowse/Store/Sequence/IndexedFasta',
                'urlTemplate': fasta_template,
                'faiUrlTemplate': fasta_template + '.fai',
                'useAsRefSeqStore': 1
            })
        if seq_source == 'twobit':
            del track['chunkSize']
            fasta_template = os.path.join('seq', os.path.basename(src))
            track.update({
                'storeClass': 'JBrowse/Store/Sequence/TwoBit',
                'urlTemplate': fasta_template,
                'useAsRefSeqStore': 1
            })
        if not data:
            data = {
                'formatVersion': 1,
                'tracks': [track]
            }
            return data
        for idx, data_track in enumerate(data['tracks'][:]):
            if data_track['label'] == seq_track_name:
                data['tracks'][idx] = track
                break
        return data
    json_store.modify('trackList.json', add_track)


def write_refseqs_json(refseqs, json_store, ref_order=None, **opts):
    """ Add JSON description of reference sequences to the seq/ directory

    :type refseqs: dict
    :param refseqs: Key-value pairs to be written to the JSON
    :param json_store: JsonFileStore object defining the output directory
    :param ref_order: Optional order for the reference sequences
    :param opts: Same as optional parameters for format_sequences()
    """
    os.makedirs(os.path.join(opts['out'], 'seq'), exist_ok=True)

    def add_refs(data):
        if data:
            for seq in data[:]:
                if seq['name'] in refseqs:
                    data.remove(seq)
        else:
            data = []
        for name in ref_order if ref_order else sorted(refseqs):
            data.append(refseqs[name])
        return data
    json_store.modify('seq/refSeqs.json', add_refs)
    if opts['compress']:
        with open(os.path.join(opts['out'], 'seq', '.htaccess'), 'w') as outfile:
            gdb = GenomeDB()
            outfile.write(gdb.precompression_htaccess('.txtz','.jsonz'))
