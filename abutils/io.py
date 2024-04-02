#!/usr/bin/env python
# filename: io.py


#
# Copyright (c) 2020 Bryan Briney
# License: The MIT license (http://opensource.org/licenses/MIT)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or
# substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
# BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


import glob
import os
import sys
from typing import Iterable, Union

from .core.sequence import (
    from_mongodb,
    read_airr,
    read_csv,
    read_fasta,
    read_fastq,
    read_json,
    to_csv,
    to_fasta,
    to_fastq,
)
from .utils.convert import abi_to_fasta

# from .utils.pipeline import list_files, make_dir


def make_dir(directory: str) -> None:
    """
    Makes a directory, if it doesn't already exist.

    Parameters
    ----------
    directory : str
        Path to a directory.

    """
    if not os.path.exists(directory):
        os.makedirs(os.path.abspath(directory))


def list_files(
    directory: str, extension: Union[str, Iterable, None] = None
) -> Iterable[str]:
    """
    Lists files in a given directory.

    Parameters
    ----------
    directory : str
        Path to a directory.

    extension : str
        If supplied, only files that end with the specificied extension(s) will be returned. Can be either
        a string or a list of strings. Extension evaluation is case-insensitive and can match complex
        extensions (e.g. '.fastq.gz'). Default is ``None``, which returns all files in the directory,
        regardless of extension.

    Returns
    -------
    Iterable[str]

    """
    if os.path.isdir(directory):
        expanded_dir = os.path.expanduser(directory)
        files = sorted(glob.glob(expanded_dir + "/*"))
    else:
        files = [
            directory,
        ]
    if extension is not None:
        if isinstance(extension, str):
            extension = [
                extension,
            ]
        files = [
            f
            for f in files
            if any(
                [
                    any([f.lower().endswith(e.lower()) for e in extension]),
                    any([f.endswith(e.upper()) for e in extension]),
                    any([f.endswith(e.lower()) for e in extension]),
                ]
            )
        ]
    return files


def rename_file(file: str, new_name: str) -> None:
    """
    Renames a file.

    Parameters
    ----------
    file : str
        Path to the file to be renamed.

    new_name : str
        New name for the file.

    """
    os.rename(file, new_name)


def delete_files(files: Union[str, Iterable]) -> None:
    """
    Deletes files.

    Parameters
    ----------
    files : Union[str, Iterable]
        Path to a file or an iterable of file paths.

    """
    if isinstance(files, str):
        files = [
            files,
        ]
    for f in files:
        if os.path.exists(f):
            os.remove(f)


def concatenate_files(files: Iterable[str], output_file: str) -> None:
    """
    Concatenates multiple files into a single file.

    Parameters
    ----------
    files : Iterable[str]
        Iterable of file paths.

    output_file : str
        Path to the output file.

    """
    with open(output_file, "w") as outfile:
        for fname in files:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)


def read_sequences(
    file=None,
    format="airr",
    sep="\t",
    fields=None,
    match=None,
    id_key="sequence_id",
    sequence_key="sequence",
    db=None,
    collection=None,
    mongodb_kwargs=None,
):
    """
    Reads sequence data from a file and returns ``abutils.Seqeunce`` objects.

    Args:
    -----

    file (str): path to a file containing sequence data in any of the supported formats.

    format (str): format of the sequence file. Supported formats are: ``'airr'``, ``'tabular'``, ``'fasta'``,
        ``'fastq'``, ``'json'`` and ``'mongodb'``. Default is ``'airr'``.

    sep (str): character used to separate fields in ``'tabular'`` input files. This option is
        only used when ``format`` is ``'tabular'``. Default is ``'\t'``, which conforms with the
        default format for AIRR-compatible sequence annotation.

    id_key (str): name of the field containing the sequence ID. Default is ``'sequence_id'``.

    sequence_key (str): name of the field containing the sequence. Default is ``'sequence'``.

    db (str): mongodb database to query for sequence information. Required if ``format`` is ``'mongodb'``.

    collection (str): mongodb collection to query for sequence information. Required if ``format`` is ``'mongodb'``.

    mongodb_kwargs (dict): dictionary containing additional keyword arguments that will be passed to
        ``abutils.io.from_mongodb``.


    Returns:
    --------

    A list of ``abutils.Sequence`` objects.
    """
    format = format.lower()
    if format == "json":
        return read_json(
            file, id_key=id_key, sequence_key=sequence_key, fields=fields, match=match
        )
    elif format == "fasta":
        return read_fasta(file)
    elif format == "fastq":
        return read_fastq(file)
    elif format == "airr":
        return read_airr(file, fields=fields, match=match)
    elif format == "tabular":
        return read_csv(
            file,
            delimiter=sep,
            id_key=id_key,
            sequence_key=sequence_key,
            fields=fields,
            match=match,
        )
    elif format == "mongodb":
        if any([db is None, collection is None]):
            error = '<db> and <collection> are required arguments if the data type is "mongodb".'
            raise ValueError(error)
        return from_mongodb(db, collection, **mongodb_kwargs)
    else:
        error = f'Format type "{format}"" is not supported. '
        error += 'Supported file types are "airr", "fasta", "json", and "tabular".'
        raise ValueError(error)
