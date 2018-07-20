# jbrowse-utils
A set of Python utilities for formatting data to be used with JBrowse

## Intro

1. What is JBrowse?

    JBrowse is an awesem modern genome browser built with JavaScript and HTML5. (See https://github.com/GMOD/jbrowse and
     http://jbrowse.org/)

2. What is this project for?

    JBrowse has a number of tools for working with JBrowse data and configurations, for example to format sequence data
    so that it can be viewed, or to add new data sources, etc. Those tools are written in Perl, and this is basically a
    rewrite of those tools in Python.

3. Why rewrite the tools in Python?

    A few reasons. It started because I had a terrible time trying to get some of the Perl dependencies that the JBrowse
    utils use to install. It sounded like a fun project to dive into the code and see if I could port some of it to
    Python, mostly for my own use, but also in case it was useful to anyone else who was sick of fighting the Perl
    installations.

 4. What's in this project?
 
    So far this is just a rewrite of `prepare-refseqs.pl`, which is used for preparing the reference sequence data for
    use in JBrowse. I think I've duplicated all the functionality except for the `--conf` option to import from a BioDB
    config file.

## Installation

The only requirement so far is Python 3+, with no other dependencies. Just clone and run.
```bash
git clone https://github.com/garrettjstevens/jbrowse-utils.git
cd jbrowse_utils
bin/jbrowse-utils --help
bin/jbrowse-utils prepare-refseqs --help
```

## Usage
```bash
$ bin/jbrowse-utils prepare-refseqs --help
usage: jbrowse-utils prepare-refseqs [-h]
                                     (--gff <GFF3 file> | --fasta <FASTA file> [<FASTA file> ...] | --indexed_fasta <FASTA file> | --twobit <2BIT file> | --conf <biodb config file> | --sizes <sizes file> [<sizes file> ...] | --gff-sizes <GFF file>)
                                     [--noSort] [--out <output directory>]
                                     [--noseq]
                                     [--refs <refseq names> [<refseq names> ...]]
                                     [--compress] [--chunksize <int>]
                                     [--nohash] [--trackLabel <label>]
                                     [--key <string>] [--seqType <string>]
                                     [--trackConfig { JSON-format extra configuration for this track }]

Formats reference sequence data for use with JBrowse

optional arguments:
  -h, --help            show this help message and exit

SEQUENCE FILE:
  Specify one of:

  --gff <GFF3 file>     Must be GFF version 3 with an embedded FASTA section
  --fasta <FASTA file> [<FASTA file> ...]
                        Can specify multiple FASTAs after one flag or use the
                        flag multiple times. Can be a gzipped file (ending in
                        .gz or .gzip). Can optionally supply --refs.
  --indexed_fasta <FASTA file>
                        An index with the same name plus ".fai" must be
                        present.
  --twobit <2BIT file>  A single .2bit file.
  --conf <biodb config file>
                        JBrowse config file pointing to a BioPerl database.
                        Must also supply --refs.
  --sizes <sizes file> [<sizes file> ...]
                        Can specify multiple sizes files after one flag or use
                        the flag multiple times.
  --gff-sizes <GFF file> [<GFF file> ...]
                        Must contain ##sequence-region lines as described in
                        the GFF specs. Can specify multiple GFF files after
                        one flag or use the flag multiple times.

OPTIONS:
  --noSort              If using GFF or FASTA input, preserve the order of the
                        reference sequences (sorts alphabetically by default).
  --out <output directory>
                        Optional directory to write to. Defaults to data/.
  --noseq               Do not store the actual sequence bases, just the
                        sequence metadata (name, length, and so forth).
  --refs <refseq names> [<refseq names> ...]
                        Output only the sequences with the given (comma- or
                        space-separated) names.
  --compress            If passed, compress the reference sequences with gzip,
                        making the chunks be .txt.gz. NOTE: this requires a
                        bit of additional web server configuration to be
                        served correctly.
  --chunksize <int>     Size of sequence chunks to make, in base pairs.
                        Default 20kb. This is multiplied by 4 if --compress is
                        passed, so that the compressed sequence files are
                        still approximately this size.
  --nohash              Store sequences in a flat seq/{seqname}/{chunk}.txt
                        structure, instead of the new (more scalable)
                        /seq/hash/hash/hash/{seqname}-{chunk}.txt structure.
  --trackLabel <label>  The unique name of the sequence track, default is to
                        set to the same as --seqType.
  --key <string>        The displayed name of the sequence track, defaults to
                        "Reference sequence".
  --seqType <string>    The name of the alphabet used for these reference
                        sequences, usually either "dna", "rna", or "protein".
                        (default "dna")
  --trackConfig { JSON-format extra configuration for this track }
                        Additional top-level configuration for the client, in
                        JSON syntax. Example: --trackConfig '{ "glyph":
                        "ProcessedTranscript" }'
```

## Credits

This is a derivative of the work in JBrowse, so all credit for the basic ideas go to them.

## License

This is released under the same license as JBrowse: GNU Lesser General Public License v2.1

Under the terms of that license, I give full credit to the original authors ()https://github.com/GMOD/jbrowse) for the
ideas used in this work. 
