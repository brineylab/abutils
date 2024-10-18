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

import pandas as pd
import polars as pl

# import pyarrow as pa
# import pyarrow.parquet as pq
from natsort import natsorted

from .core.pair import (
    Pair,
    pairs_from_pandas,
    pairs_from_polars,
    pairs_to_csv,
    pairs_to_pandas,
    pairs_to_parquet,
    pairs_to_polars,
)
from .core.sequence import (
    Sequence,
    determine_fastx_format,
    from_mongodb,
    parse_fasta,
    parse_fastq,
    parse_fastx,
    read_airr,
    # read_csv,
    read_fasta,
    read_fastq,
    read_fastx,
    read_json,
    # read_parquet,
    sequences_from_pandas,
    sequences_from_polars,
    sequences_to_csv,
    sequences_to_pandas,
    sequences_to_parquet,
    sequences_to_polars,
    to_airr,
    to_fasta,
    to_fastq,
)
from .utils.convert import abi_to_fasta

# to_airr.__module__ = "abutils.io"

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
    ignore_dot_files: bool = True,
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
    if ignore_dot_files:
        files = [f for f in files if not os.path.basename(f).startswith(".")]
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


# --------------------------------------------
#      File splitting and concatenation
# --------------------------------------------


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
    import pyarrow as pa  # for some reason, pyarrow breaks readthedocs if imported up top
    import pyarrow.parquet as pq

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


def split_fastx(
    fastx_file: str,
    output_directory: str,
    chunksize: int = 500,
    start_numbering_at: int = 0,
    fmt: Optional[str] = None,
) -> Iterable[str]:
    """
    Splits a FASTA or FASTQ file into multiple files, each containing a specified number of sequences.

    Parameters
    ----------
    fastx_file : str
        Path to the FASTA or FASTQ file to be split.

    output_directory : str
        Path to the directory where the split files will be saved. If the directory does not exist, it will be created.

    chunksize : int, optional
        Number of sequences per split file. Default is 500. The last file may contain fewer sequences than this number.

    start_numbering_at : int, optional
        Start numbering the split files at this number. Default is 0.

    fmt : str, optional
        Format of the input file. If not supplied, the format will be determined automatically.

    Returns
    -------
    Iterable[str]
        Iterable of file paths for the split files.

    """
    # outputs
    output_files = []
    output_directory = os.path.abspath(output_directory)
    make_dir(output_directory)
    output_basename = ".".join(
        os.path.basename(fastx_file).rstrip(".gz").split(".")[:-1]
    )

    # get input file format
    if fmt is None:
        fmt = determine_fastx_format(fastx_file)

    # split sequences into chunks, write the chunks to files
    file_num = start_numbering_at
    chunk = []
    for s in parse_fastx(fastx_file):
        chunk.append(s.fastq if fmt == "fastq" else s.fasta)
        if len(chunk) == chunksize:
            output_file = os.path.join(
                output_directory, f"{output_basename}_{file_num}.{fmt}"
            )
            with open(output_file, "w") as f:
                f.write("\n".join(chunk))
            chunk = []
            file_num += 1
            output_files.append(output_file)

    # write the last chunk to a file
    if chunk:
        output_file = os.path.join(
            output_directory, f"{output_basename}_{file_num}.{fmt}"
        )
        with open(output_file, "w") as f:
            f.write("\n".join(chunk))
        output_files.append(output_file)

    return output_files


def split_fasta(
    fasta_file: str,
    output_directory: str,
    chunksize: int = 500,
    start_numbering_at: int = 0,
) -> Iterable[str]:
    """
    Splits a FASTA or FASTQ file into multiple files, each containing a specified number of sequences.

    Parameters
    ----------
    fasta_file : str
        Path to the FASTA file to be split.

    output_directory : str
        Path to the directory where the split files will be saved. If the directory does not exist, it will be created.

    chunksize : int, optional
        Number of sequences per split file. Default is 500. The last file may contain fewer sequences than this number.

    start_numbering_at : int, optional
        Start numbering the split files at this number. Default is 0.

    Returns
    -------
    Iterable[str]
        Iterable of file paths for the split files.

    """
    return split_fastx(
        fastx_file=fasta_file,
        output_directory=output_directory,
        chunksize=chunksize,
        start_numbering_at=start_numbering_at,
        fmt="fasta",
    )


def split_fastq(
    fastq_file: str,
    output_directory: str,
    chunksize: int = 500,
    start_numbering_at: int = 0,
) -> Iterable[str]:
    """
    Splits a FASTA or FASTQ file into multiple files, each containing a specified number of sequences.

    Parameters
    ----------
    fastq_file : str
        Path to the FASTQ file to be split.

    output_directory : str
        Path to the directory where the split files will be saved. If the directory does not exist, it will be created.

    chunksize : int, optional
        Number of sequences per split file. Default is 500. The last file may contain fewer sequences than this number.

    start_numbering_at : int, optional
        Start numbering the split files at this number. Default is 0.

    Returns
    -------
    Iterable[str]
        Iterable of file paths for the split files.

    """
    return split_fastx(
        fastx_file=fastq_file,
        output_directory=output_directory,
        chunksize=chunksize,
        start_numbering_at=start_numbering_at,
        fmt="fastq",
    )


# ------------------------------
#     Sequence and Pair I/O
# ------------------------------


def read_csv(
    csv_file: str,
    separator: str = ",",
    match: Optional[dict] = None,
    fields: Optional[Iterable] = None,
    id_key: str = "sequence_id",
    sequence_key: str = "sequence",
) -> Iterable[Union[Sequence, Pair]]:
    """
    Reads a CSV file and returns a list of ``Sequence`` or ``Pair`` objects.

    Parameters
    ----------
    csv_file : str
        Path to the CSV file to be read.

    separator : str, default=","
        Column separator. Default is ``","``.

    match : dict, optional
        A ``dict`` for filtering sequences from the input file.
        Sequences must match all conditions to be returned. For example,
        the following ``dict`` will filter out all sequences for which
        the ``'v_gene:0'`` field is not ``'IGHV1-2'``:

        .. code-block:: python

            {'v_gene:0': 'IGHV1-2'}

    fields : list, optional
        A list of fields to be read from the input file. If not provided, all fields will be read.

    id_key : str, default="sequence_id"
        The name of the column that contains the sequence IDs. Default is ``"sequence_id"``.

    sequence_key : str, default="sequence"
        The name of the column that contains the sequence data. Default is ``"sequence"``.

    Returns
    -------
    Iterable[Union[Sequence, Pair]]
        A list of ``Sequence`` or ``Pair`` objects.

    """
    df = pl.read_csv(csv_file, separator=separator, infer_schema_length=None)
    if any([c in df.columns for c in ["sequence:0", "sequence_id:0"]]):
        return pairs_from_polars(
            df, match=match, fields=fields, id_key=id_key, sequence_key=sequence_key
        )
    else:
        return sequences_from_polars(
            df, match=match, fields=fields, id_key=id_key, sequence_key=sequence_key
        )


def to_csv(
    sequences: Iterable[Union[Sequence, Pair]],
    csv_file: Optional[str] = None,
    separator: str = ",",
    header: bool = True,
    columns: Optional[Iterable] = None,
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

    separator : str, default=","
        Column separator. Default is ``","``.

    header : bool, default=True
        If ``True``, the CSV file will contain a header row. Default is ``True``.

    columns : list, default=None
        A list of fields to be retained in the output CSV file. Fields must be column
        names in the input file.

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
            separator=separator,
            header=header,
            columns=columns,
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
            separator=separator,
            header=header,
            columns=columns,
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


def read_parquet(
    parquet_file: str,
    match: Optional[dict] = None,
    fields: Optional[Iterable] = None,
    id_key: str = "sequence_id",
    sequence_key: str = "sequence",
) -> Iterable[Union[Sequence, Pair]]:
    """
    Reads a Parquet file and returns a list of ``Sequence`` or ``Pair`` objects.

    Parameters
    ----------
    parquet_file : str
        Path to the Parquet file to be read.

    match : dict, optional
        A ``dict`` for filtering sequences from the input file.
        Sequences must match all conditions to be returned. For example,
        the following ``dict`` will filter out all sequences for which
        the ``'v_gene:0'`` field is not ``'IGHV1-2'``:

        .. code-block:: python

            {'v_gene:0': 'IGHV1-2'}

    fields : list, optional
        A list of fields to be read from the input file. If not provided, all fields will be read.

    id_key : str, default="sequence_id"
        The name of the column that contains the sequence IDs. Default is ``"sequence_id"``.

    sequence_key : str, default="sequence"
        The name of the column that contains the sequence data. Default is ``"sequence"``.

    Returns
    -------
    Iterable[Union[Sequence, Pair]]
        A list of ``Sequence`` or ``Pair`` objects.

    """
    df = pl.read_parquet(parquet_file)
    if any([c in df.columns for c in ["sequence:0", "sequence_id:0"]]):
        return pairs_from_polars(
            df, match=match, fields=fields, id_key=id_key, sequence_key=sequence_key
        )
    else:
        return sequences_from_polars(
            df, match=match, fields=fields, id_key=id_key, sequence_key=sequence_key
        )


def to_parquet(
    sequences: Iterable[Union[Sequence, Pair]],
    parquet_file: Optional[str] = None,
    columns: Optional[Iterable] = None,
    properties: Optional[Iterable[str]] = None,
    sequence_properties: Optional[Iterable[str]] = None,
    drop_na_columns: bool = True,
    order: Optional[Iterable[str]] = None,
    exclude: Optional[Union[str, Iterable[str]]] = None,
    leading: Optional[Union[str, Iterable[str]]] = None,
) -> Optional[pl.DataFrame]:
    """
    Saves a list of ``Pair`` objects to a Parquet file.

    Parameters
    ----------
    sequences : Iterable[Sequence, Pair]
        List of ``Sequence`` or ``Pair`` objects to be saved to a Parquet file. Required.

    parquet_file : str
        Path to the output Parquet file. If ``None``, the Parquet file will not be saved to a file and the
        function will return a ``polars.DataFrame`` object.

    columns : list, default=None
        A list of fields to be retained in the output Parquet file. Fields must be column
        names in the input file.

    properties : list, default=None
        A list of properties to be included in the Parquet file. If not provided, everything
        in the ``annotations`` field of each heavy/light chain will be included.

    sequence_properties : list, default=None
        A list of sequence properties to be included. Differs from ``properties``, which
        refers to properties of the ``Pair`` object. These properties are those of the
        heavy/light ``Sequence`` objects. Ignored if the input is a list of ``Sequence`` objects.

    drop_na_columns : bool, default=True
        If ``True``, columns with all ``NaN`` values will be dropped from the Parquet file.
        Default is ``True``.

    order : list, default=None
        A list of fields in the order they should appear in the Parquet file.

    exclude : str or list, default=None
        Field or list of fields to be excluded from the Parquet file.

    leading : str or list, default=None
        Field or list of fields to appear first in the Parquet file. Supercedes ``order``, so
        if both are provided, fields in ``leading`` will appear first in the Parquet file and
        remaining fields will appear in the order provided in ``order``.

    Returns
    -------
    Optional[pl.DataFrame]
        A ``polars.DataFrame`` object if ``csv_file`` is ``None``.

    """
    if all([isinstance(s, Sequence) for s in sequences]):
        return sequences_to_parquet(
            sequences,
            parquet_file=parquet_file,
            columns=columns,
            properties=properties,
            drop_na_columns=drop_na_columns,
            order=order,
            exclude=exclude,
            leading=leading,
        )
    elif all([isinstance(s, Pair) for s in sequences]):
        return pairs_to_parquet(
            sequences,
            parquet_file=parquet_file,
            columns=columns,
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


def from_polars(
    df: Union[pl.DataFrame, pl.LazyFrame],
    match: Optional[dict] = None,
    fields: Optional[Iterable] = None,
    id_key: str = "sequence_id",
    sequence_key: str = "sequence",
) -> Iterable[Union[Sequence, Pair]]:
    """
    Reads a Polars DataFrame and returns a list of ``Sequence`` or ``Pair`` objects.

    Parameters
    ----------
    df : polars.DataFrame or polars.LazyFrame
        The input Polars DataFrame or LazyFrame.

    match : dict, optional
        A ``dict`` for filtering sequences from the input file.
        Sequences must match all conditions to be returned. For example,
        the following ``dict`` will filter out all sequences for which
        the ``'v_gene:0'`` field is not ``'IGHV1-2'``:

        .. code-block:: python

            {'v_gene:0': 'IGHV1-2'}

    fields : list, optional
        A list of fields to be read from the input file. If not provided, all fields will be read.

    id_key : str, default="sequence_id"
        The name of the column that contains the sequence IDs. Default is ``"sequence_id"``.

    sequence_key : str, default="sequence"
        The name of the column that contains the sequence data. Default is ``"sequence"``.

    Returns
    -------
    Iterable[Union[Sequence, Pair]]
        A list of ``Sequence`` or ``Pair`` objects.

    """
    if isinstance(df, pl.LazyFrame):
        df = df.collect()
    if any([c in df.columns for c in ["sequence:0", "sequence_id:0"]]):
        return pairs_from_polars(
            df, match=match, fields=fields, id_key=id_key, sequence_key=sequence_key
        )
    else:
        return sequences_from_polars(
            df, match=match, fields=fields, id_key=id_key, sequence_key=sequence_key
        )


def to_polars(
    sequences: Iterable[Union[Sequence, Pair]],
    columns: Optional[Iterable] = None,
    properties: Optional[Iterable[str]] = None,
    sequence_properties: Optional[Iterable[str]] = None,
    drop_na_columns: bool = True,
    order: Optional[Iterable[str]] = None,
    exclude: Optional[Union[str, Iterable[str]]] = None,
    leading: Optional[Union[str, Iterable[str]]] = None,
) -> Optional[pl.DataFrame]:
    """
    Saves a list of ``Pair`` objects to a Polars DataFrame.

    Parameters
    ----------
    sequences : Iterable[Sequence, Pair]
        List of ``Sequence`` or ``Pair`` objects to be saved to a Polars DataFrame. Required.

    columns : list, default=None
        A list of fields to be retained in the output Polars DataFrame. Fields must be column
        names in the input file.

    properties : list, default=None
        A list of properties to be included in the Polars DataFrame. If not provided, everything
        in the ``annotations`` field of each heavy/light chain will be included.

    sequence_properties : list, default=None
        A list of sequence properties to be included. Differs from ``properties``, which
        refers to properties of the ``Pair`` object. These properties are those of the
        heavy/light ``Sequence`` objects. Ignored if the input is a list of ``Sequence`` objects.

    drop_na_columns : bool, default=True
        If ``True``, columns with all ``NaN`` values will be dropped from the Polars DataFrame.
        Default is ``True``.

    order : list, default=None
        A list of fields in the order they should appear in the Polars DataFrame.

    exclude : str or list, default=None
        Field or list of fields to be excluded from the Polars DataFrame.

    leading : str or list, default=None
        Field or list of fields to appear first in the Polars DataFrame. Supercedes ``order``, so
        if both are provided, fields in ``leading`` will appear first in the Polars DataFrame and
        remaining fields will appear in the order provided in ``order``.

    Returns
    -------
    pl.DataFrame
        A ``polars.DataFrame`` object.

    """
    if all([isinstance(s, Sequence) for s in sequences]):
        return sequences_to_polars(
            sequences,
            columns=columns,
            properties=properties,
            drop_na_columns=drop_na_columns,
            order=order,
            exclude=exclude,
            leading=leading,
        )
    elif all([isinstance(s, Pair) for s in sequences]):
        return pairs_to_polars(
            sequences,
            columns=columns,
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


def from_pandas(
    df: pd.DataFrame,
    match: Optional[dict] = None,
    fields: Optional[Iterable] = None,
    id_key: str = "sequence_id",
    sequence_key: str = "sequence",
) -> Iterable[Union[Sequence, Pair]]:
    """
    Reads a Pandas DataFrame and returns a list of ``Sequence`` or ``Pair`` objects.

    Parameters
    ----------
    df : pd.DataFrame
        The input Pandas DataFrame.

    match : dict, optional
        A ``dict`` for filtering sequences from the input file.
        Sequences must match all conditions to be returned. For example,
        the following ``dict`` will filter out all sequences for which
        the ``'v_gene:0'`` field is not ``'IGHV1-2'``:

        .. code-block:: python

            {'v_gene:0': 'IGHV1-2'}

    fields : list, optional
        A list of fields to be read from the input file. If not provided, all fields will be read.

    id_key : str, default="sequence_id"
        The name of the column that contains the sequence IDs. Default is ``"sequence_id"``.

    sequence_key : str, default="sequence"
        The name of the column that contains the sequence data. Default is ``"sequence"``.

    Returns
    -------
    Iterable[Union[Sequence, Pair]]
        A list of ``Sequence`` or ``Pair`` objects.

    """
    df = pl.from_pandas(df)
    if any([c in df.columns for c in ["sequence:0", "sequence_id:0"]]):
        return pairs_from_polars(
            df, match=match, fields=fields, id_key=id_key, sequence_key=sequence_key
        )
    else:
        return sequences_from_polars(
            df, match=match, fields=fields, id_key=id_key, sequence_key=sequence_key
        )


def to_pandas(
    sequences: Iterable[Union[Sequence, Pair]],
    columns: Optional[Iterable] = None,
    properties: Optional[Iterable[str]] = None,
    sequence_properties: Optional[Iterable[str]] = None,
    drop_na_columns: bool = True,
    order: Optional[Iterable[str]] = None,
    exclude: Optional[Union[str, Iterable[str]]] = None,
    leading: Optional[Union[str, Iterable[str]]] = None,
) -> Optional[pl.DataFrame]:
    """
    Saves a list of ``Pair`` objects to a Pandas DataFrame.

    Parameters
    ----------
    sequences : Iterable[Sequence, Pair]
        List of ``Sequence`` or ``Pair`` objects to be saved to a Pandas DataFrame. Required.

    columns : list, default=None
        A list of fields to be retained in the output Pandas DataFrame. Fields must be column
        names in the input file.

    properties : list, default=None
        A list of properties to be included in the Pandas DataFrame. If not provided, everything
        in the ``annotations`` field of each heavy/light chain will be included.

    sequence_properties : list, default=None
        A list of sequence properties to be included. Differs from ``properties``, which
        refers to properties of the ``Pair`` object. These properties are those of the
        heavy/light ``Sequence`` objects. Ignored if the input is a list of ``Sequence`` objects.

    drop_na_columns : bool, default=True
        If ``True``, columns with all ``NaN`` values will be dropped from the Pandas DataFrame.
        Default is ``True``.

    order : list, default=None
        A list of fields in the order they should appear in the Pandas DataFrame.

    exclude : str or list, default=None
        Field or list of fields to be excluded from the Pandas DataFrame.

    leading : str or list, default=None
        Field or list of fields to appear first in the Pandas DataFrame. Supercedes ``order``, so
        if both are provided, fields in ``leading`` will appear first in the Pandas DataFrame and
        remaining fields will appear in the order provided in ``order``.

    Returns
    -------
    pl.DataFrame
        A ``polars.DataFrame`` object.

    """
    if all([isinstance(s, Sequence) for s in sequences]):
        return sequences_to_pandas(
            sequences,
            columns=columns,
            properties=properties,
            drop_na_columns=drop_na_columns,
            order=order,
            exclude=exclude,
            leading=leading,
        )
    elif all([isinstance(s, Pair) for s in sequences]):
        return pairs_to_pandas(
            sequences,
            columns=columns,
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
