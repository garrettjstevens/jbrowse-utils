import collections
import json
import os
import os.path

from .genome_db import GenomeDB


class JsonFileStore(object):
    """Manage a directory structure of .json or .jsonz files"""

    def __init__(self, outdir, compress, **kwargs):
        """Initialize a directory to store JSON files

        :param outdir: Directory in which JSON files will be placed
        :type compress: bool
        :param compress: Whether to compress files to .jsonz
        :param kwargs: Additional optional arguments to pass to the
            json.dump() call
        """
        if 'separators' in kwargs:
            self.separators = kwargs.pop('separators')
        else:
            self.separators = (',', ':')
        if 'cls' in kwargs:
            self.cls = kwargs.pop('cls')
        else:
            self.cls = BytesEncoder
        self.kwargs = kwargs
        self.outdir = outdir
        self.compress = compress
        self.ext = '.jsonz' if compress else '.json'
        self.htaccess_written = False
        os.makedirs(outdir, exist_ok=True)

    def modify(self, filename, callback):
        """Modify an existing file in the directory

        :param filename: Name of the file to modify
        :param callback: function
        :param callback: A function that defines the modification to be done.
            The function should take the loaded JSON data and return the
            data to be dumped back into the JSON.
        """
        self._write_htaccess()

        data = None
        fn = self.full_path(filename)
        if os.path.isfile(fn) and os.stat(fn).st_size > 0:
            with open(fn) as infile:
                data = json.load(
                    infile, object_pairs_hook=collections.OrderedDict)
        data = callback(data)
        with open(fn, 'w') as outfile:
            json.dump(
                data, outfile,
                separators=self.separators, cls=self.cls, **self.kwargs)

    def touch(self, filename):
        """Create an empty JSON file or do nothing if it exists

        :param filename: Name of the JSON to create
        """
        with open(self.full_path(filename), 'a'):
            pass

    def full_path(self, path):
        """ Join the output directory with the path in the directory

        :param path: Path within the output directory
        """
        return os.path.join(self.outdir, path)

    def _write_htaccess(self):
        """ If JSON files are compressed, write a ".htaccess" file to the
        directory that will automatically configure users with
        Apache (and AllowOverride on) to serve the compressed files correctly
        """
        if self.compress and not self.htaccess_written:
            with open(os.path.join(self.outdir, '.htaccess'), 'w') as outfile:
                gdb = GenomeDB()
                outfile.write(
                    gdb.precompression_htaccess('.jsonz', '.txtz', '.txt.gz'))
            self.htaccess_written = True


class BytesEncoder(json.JSONEncoder):
    """Convert bytes to strings to allow dumping to JSON"""
    def default(self, o):
        if isinstance(o, bytes):
            return o.decode()
        return json.JSONEncoder.default(self, o)
