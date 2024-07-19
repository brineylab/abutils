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
import gzip
import os
from typing import Iterable, Optional, Union

import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq
from natsort import natsorted

from .core.pair import Pair, pairs_to_csv
from .core.sequence import (
    Sequence,
    determine_fastx_format,
    from_mongodb,
    parse_fasta,
    parse_fastq,
    parse_fastx,
    read_airr,
    read_csv,
    read_fasta,
    read_fastq,
    read_fastx,
    read_json,
    sequences_to_csv,
    to_fasta,
    to_fastq,
)
from .utils.convert import abi_to_fasta

# -----------------------
#   General file I/O
# -----------------------


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
    directory: str,
    extension: Union[str, Iterable, None] = None,
    recursive: bool = False,
) -> Iterable[str]:
    """
    Lists files in a given directory.

    Parameters
    ----------
    directory : str
        Path to a directory. If a file path is passed instead, the returned list of files will contain
        only that file path.

    extension : str
        If supplied, only files that end with the specificied extension(s) will be returned. Can be either
        a string or a list of strings. Extension evaluation is case-insensitive and can match complex
        extensions (e.g. '.fastq.gz'). Default is ``None``, which returns all files in the directory,
        regardless of extension.

    Returns
    -------
    Iterable[str]

    """
    directory = os.path.abspath(directory)
    if os.path.exists(directory):
        if os.path.isdir(directory):
            if recursive:
                files = []
                for root, dirs, _files in os.walk(directory):
                    for f in _files:
                        file_path = os.path.join(root, f)
                        files.append(file_path)
            else:
                files = natsorted(glob.glob(directory + "/*"))
        else:
            files = [directory]
    else:
        raise ValueError(f"Directory {directory} does not exist.")
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


def split_parquet(
    parquet_file: str,
    output_directory: str,
    num_rows: int = 500,
    num_splits: Optional[int] = None,
    split_prefix: str = "chunk_",
    start_numbering_at: int = 0,
) -> Iterable[str]:
    """
    Splits a parquet file into multiple files.

    Parameters
    ----------
    parquet_file : str
        Path to the parquet file to be split.

    output_directory : str
        Path to the directory where the split files will be saved. If the directory does not exist, it will be created.

    num_rows : int, optional
        Number of rows per split file. Default is 500. If ``num_splits`` is supplied, this argument is ignored.

    num_splits : int, optional
        Number of split files to create. If not supplied, ``num_rows`` is used to determine the number of split files.

    split_prefix : str, optional
        Prefix for the split files, which is followed directly by the file number. Default is "chunk_".

    start_numbering_at : int, optional
        Start numbering the split files at this number. Default is 0.

    Returns
    -------
    Iterable[str]
        Iterable of file paths for the split files.

    """
    # outputs
    output_files = []
    output_directory = os.path.abspath(output_directory)
    make_dir(output_directory)
    # read parquet file
    pqfile = pq.ParquetFile(parquet_file)
    # number of rows per split file
    if num_splits is not None:
        total_rows = pqfile.metadata.num_rows
        num_rows = total_rows // num_splits
        if total_rows % num_splits != 0:
            num_rows += 1
    # split into batches and process
    batches = pqfile.iter_batches(batch_size=num_rows)
    for i, batch in enumerate(batches, start=start_numbering_at):
        table = pa.Table.from_batches([batch])
        output_file = os.path.join(output_directory, f"{split_prefix}{i}.parquet")
        pq.write_table(table, output_file)
        output_files.append(output_file)
    return output_files


# -----------------------
#     Sequence I/O
# -----------------------


def to_csv(
    sequences: Iterable[Union[Sequence, Pair]],
    csv_file: Optional[str] = None,
    sep: str = ",",
    header: bool = True,
    columns: Optional[Iterable] = None,
    index: bool = False,
    properties: Optional[Iterable[str]] = None,
    sequence_properties: Optional[Iterable[str]] = None,
    drop_na_columns: bool = True,
    order: Optional[Iterable[str]] = None,
    exclude: Optional[Union[str, Iterable[str]]] = None,
    leading: Optional[Union[str, Iterable[str]]] = None,
) -> Optional[pl.DataFrame]:
    """
    Saves a list of ``Pair`` objects to a CSV file.

    Parameters
    ----------
    sequences : Iterable[Sequence, Pair]
        List of ``Sequence`` or ``Pair`` objects to be saved to a CSV file. Required.

    csv_file : str
        Path to the output CSV file. If ``None``, the CSV will not be saved to a file and the
        function will return a ``polars.DataFrame`` object.

    sep : str, default=","
        Column delimiter. Default is ``","``.

    header : bool, default=True
        If ``True``, the CSV file will contain a header row. Default is ``True``.

    columns : list, default=None
        A list of fields to be retained in the output CSV file. Fields must be column
        names in the input file.

    index : bool, default=False
        If ``True``, the CSV file will contain an index column. Default is ``False``.

    properties : list, default=None
        A list of properties to be included in the CSV file. If not provided, everything
        in the ``annotations`` field of each heavy/light chain will be included.

    sequence_properties : list, default=None
        A list of sequence properties to be included. Differs from ``properties``, which
        refers to properties of the ``Pair`` object. These properties are those of the
        heavy/light ``Sequence`` objects. Ignored if the input is a list of ``Sequence`` objects.

    drop_na_columns : bool, default=True
        If ``True``, columns with all ``NaN`` values will be dropped from the CSV file.
        Default is ``True``.

    order : list, default=None
        A list of fields in the order they should appear in the CSV file.

    exclude : str or list, default=None
        Field or list of fields to be excluded from the CSV file.

    leading : str or list, default=None
        Field or list of fields to appear first in the CSV file. Supercedes ``order``, so
        if both are provided, fields in ``leading`` will appear first in the CSV file and
        remaining fields will appear in the order provided in ``order``.

    Returns
    -------
    Optional[pl.DataFrame]
        A ``polars.DataFrame`` object if ``csv_file`` is ``None``.

    """
    if all([isinstance(s, Sequence) for s in sequences]):
        return sequences_to_csv(
            sequences,
            csv_file=csv_file,
            sep=sep,
            header=header,
            columns=columns,
            index=index,
            properties=properties,
            drop_na_columns=drop_na_columns,
            order=order,
            exclude=exclude,
            leading=leading,
        )
    elif all([isinstance(s, Pair) for s in sequences]):
        return pairs_to_csv(
            sequences,
            csv_file=csv_file,
            sep=sep,
            header=header,
            columns=columns,
            index=index,
            properties=properties,
            sequence_properties=sequence_properties,
            drop_na_columns=drop_na_columns,
            order=order,
            exclude=exclude,
            leading=leading,
        )
    else:
        raise ValueError(
            "All elements in the input list must be of the same type (either Sequence or Pair)."
        )


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
