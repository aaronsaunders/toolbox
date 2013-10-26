"""
Mappings from file extensions to biopython types.

Copied from Erik Matsens seqmagick https://github.com/fhcrc/seqmagick/

"""
import argparse
import contextlib
import copy
import functools
import os
import os.path
import signal
import sys
import tempfile

import bz2
import gzip
import os.path
import sys

# Define mappings in a dictionary with extension : BioPython_file_type.
EXTENSION_TO_TYPE = {'.aln': 'clustal',
                     '.afa': 'fasta',
                     '.fa': 'fasta',
                     '.faa': 'fasta',
                     '.fas': 'fasta',
                     '.fasta': 'fasta',
                     '.fastq': 'fastq',
                     '.fq': 'fastq',
                     '.ffn': 'fasta',
                     '.fna': 'fasta',
                     '.frn': 'fasta',
                     '.gb': 'genbank',
                     '.gbk': 'genbank',
                     '.needle': 'emboss',
                     '.nex': 'nexus',
                     '.phy': 'phylip',
                     '.phylip': 'phylip',
                     '.phyx': 'phylip-relaxed',
                     '.qual': 'qual',
                     '.sff': 'sff-trim',
                     '.sth': 'stockholm',
                     '.sto': 'stockholm',}

COMPRESS_EXT = {'.bz2': bz2.BZ2File, '.gz': gzip.open, '.bz': bz2.BZ2File}


class UnknownExtensionError(ValueError):
    pass


def from_extension(extension):
    """
    Look up the BioPython file type corresponding with input extension.

    Look up is case insensitive.
    """
    if not extension.startswith('.'):
        raise ValueError("Extensions must begin with a period.")
    try:
        return EXTENSION_TO_TYPE[extension.lower()]
    except KeyError:
        raise UnknownExtensionError("seqmagick does not know how to handle " +
                "files with extensions like this: " + extension)


def from_filename(file_name):
    """
    Look up the BioPython file type corresponding to an input file name.
    """
    base, extension = os.path.splitext(file_name)
    if extension in COMPRESS_EXT:
        # Compressed file
        extension = os.path.splitext(base)[1]
    return from_extension(extension)

def from_handle(fh, stream_default='fasta'):
    """
    Look up the BioPython file type corresponding to a file-like object.

    For stdin, stdout, and stderr, ``stream_default`` is used.
    """
    if fh in (sys.stdin, sys.stdout, sys.stderr):
        return stream_default
    return from_filename(fh.name)


class FileType(object):
    """
    Near clone of argparse.FileType, supporting gzip and bzip
    """
    def __init__(self, mode='r'):
        self.mode = mode
        self.ext_map = COMPRESS_EXT.copy()

    def _get_handle(self, file_path):
        ext = os.path.splitext(file_path)[1].lower()
        return self.ext_map.get(ext, open)(file_path, self.mode)

    def __call__(self, string):
        if string == '-':
            if 'r' in self.mode:
                return sys.stdin
            elif 'w' in self.mode:
                return sys.stdout
            else:
                raise ValueError("Invalid mode: {0}".format(string))
        else:
            return self._get_handle(string)
