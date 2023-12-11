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


import sys

from .core.sequence import (
    read_csv,
    read_airr,
    read_fasta,
    read_json,
    from_mongodb,
    to_fasta,
    to_fastq,
)

from .utils.convert import abi_to_fasta
from .utils.pipeline import list_files, make_dir


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
        ``'json'`` and ``'mongodb'``. Default is ``'airr'``.

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

