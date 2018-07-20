import unittest
import os.path
import tempfile

from jbrowse_utils import prepare_refseqs


class TestPrepareRefseqs(unittest.TestCase):

    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), 'data')

    def test_gff(self):
        with tempfile.TemporaryDirectory() as tmpdirname:
            gff = os.path.join(self.data_dir, 'embed_fasta.gff3')
            prepare_refseqs.format_sequences(gff=gff, out=tmpdirname)

    def test_fastas(self):
        with tempfile.TemporaryDirectory() as tmpdirname:
            fastas = [os.path.join(self.data_dir, fa) for fa in ('a.fa', 'b.fa')]
            prepare_refseqs.format_sequences(fastas=fastas, out=tmpdirname)

    def test_indexed_fasta(self):
        with tempfile.TemporaryDirectory() as tmpdirname:
            fasta = os.path.join(self.data_dir, 'a.fa')
            prepare_refseqs.format_sequences(indexed_fasta=fasta, out=tmpdirname)

    def test_twobit(self):
        with tempfile.TemporaryDirectory() as tmpdirname:
            twobit = os.path.join(self.data_dir, 'a.2bit')
            prepare_refseqs.format_sequences(twobit=twobit, out=tmpdirname)

    def test_sizes(self):
        with tempfile.TemporaryDirectory() as tmpdirname:
            sizes = os.path.join(self.data_dir, 'a.sizes')
            prepare_refseqs.format_sequences(sizes=[sizes], out=tmpdirname)

    def test_gff_sizes(self):
        with tempfile.TemporaryDirectory() as tmpdirname:
            gff_sizes = os.path.join(self.data_dir, 'sizes.gff3')
            prepare_refseqs.format_sequences(gff_sizes=[gff_sizes], out=tmpdirname)


if __name__ == '__main__':
    unittest.main()
