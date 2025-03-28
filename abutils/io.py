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
import re
import tempfile
from typing import Iterable, Optional, Union

import pandas as pd
import polars as pl
from Bio.SeqRecord import SeqRecord

# import pyarrow as pa
# import pyarrow.parquet as pq
from natsort import natsorted

from .core.pair import (
    Pair,
    # pairs_from_pandas,
    pairs_from_polars,
    # pairs_to_csv,
    # pairs_to_pandas,
    # pairs_to_parquet,
    pairs_to_polars,
)
from .core.sequence import (
    Sequence,
    determine_fastx_format,
    from_mongodb,
    parse_fasta,
    parse_fastq,
    parse_fastx,
    # read_airr,
    # read_csv,
    read_fasta,
    read_fastq,
    read_fastx,
    read_json,
    # read_parquet,
    # sequences_from_pandas,
    sequences_from_polars,
    # sequences_to_csv,
    # sequences_to_pandas,
    # sequences_to_parquet,
    sequences_to_polars,
    # to_airr,
    # to_fasta,
    to_fastq,
)
from .utils.path import delete_files, list_files, make_dir, rename_file

# =======================
#   General file I/O
# =======================


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
    match: Optional[str] = None,
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

    recursive : bool, default=False
        If ``True``, the directory will be searched recursively, and all files in all subdirectories will be returned.

    match : str, optional
        If supplied, only files that match the specified pattern will be returned. Regular expressions are supported.

    ignore_dot_files : bool, default=True
        If ``True``, dot files (hidden files) will be ignored.

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
    if match is not None:
        files = [f for f in files if re.match(match, os.path.basename(f)) is not None]
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


# ============================================
#      File splitting and concatenation
# ============================================


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


def concatenate_parquet(parquet_files: Iterable[str], output_file: str) -> str:
    """
    Concatenates multiple parquet files into a single file.

    Parameters
    ----------
    parquet_files : Iterable[str]
        Iterable of file paths to the parquet files to be concatenated.

    output_file : str
        Path to the output file.

    Returns
    -------
    str
        Path to the output file.

    """
    make_dir(os.path.dirname(output_file))

    # do all parquet operations lazily to save memory
    dfs = [pl.scan_parquet(f) for f in parquet_files]
    df = pl.concat(dfs)
    df.sink_parquet(output_file)
    return output_file


def split_csv(
    csv_file: str, output_directory: str, separator: str = ",", chunksize: int = 500
) -> Iterable[str]:
    """
    Splits a CSV file into multiple files, each containing a specified number of rows, including the header in each split.

    Parameters
    ----------
    csv_file : str
        Path to the CSV file to be split.

    output_directory : str
        Path to the directory where the split files will be saved. If the directory does not exist, it will be created.

    separator : str, default=","
        Column separator. Default is `","`.

    chunksize : int, default=500
        Number of rows per split file. Default is 500. The last file may contain fewer rows than this number.

    Returns
    -------
    Iterable[str]
        Iterable of file paths for the split files.

    """
    make_dir(output_directory)

    output_files = []
    output_basename = os.path.splitext(os.path.basename(csv_file))[0]

    for i, chunk in enumerate(
        pd.read_csv(csv_file, sep=separator, chunksize=chunksize)
    ):
        output_file = os.path.join(output_directory, f"{output_basename}_part_{i}.csv")
        chunk.to_csv(output_file, sep=separator, index=False, mode="w", header=True)
        output_files.append(output_file)

    return output_files


def concatenate_csv(
    csv_files: Iterable[str], output_file: str, separator: str = ","
) -> str:
    """
    Concatenates multiple CSV files into a single file.

    Parameters
    ----------
    csv_files : Iterable[str]
        Iterable of file paths to the CSV files to be concatenated.

    output_file : str
        Path to the output file.

    separator : str, default=","
        Column separator. Default is ``","``.

    Returns
    -------
    str
        Path to the output file.

    """
    make_dir(os.path.dirname(output_file))

    # do all csv operations lazily to save memory
    dfs = [pl.scan_csv(f, separator=separator) for f in csv_files]
    df = pl.concat(dfs)
    df.sink_csv(output_file, separator=separator)
    return output_file


def split_airr(
    airr_file: str, output_directory: str, chunksize: int = 500
) -> Iterable[str]:
    """
    Splits an AIRR-formatted file into multiple files.

    Parameters
    ----------
    airr_file : str
        Path to the AIRR-formatted file to be split.

    output_directory : str
        Path to the directory where the split files will be saved.

    chunksize : int, default=500
        Number of rows per split file. Default is 500. The last file may contain fewer rows than this number.

    Returns
    -------
    Iterable[str]
        Iterable of file paths for the split files.

    """
    return split_csv(airr_file, output_directory, separator="\t", chunksize=chunksize)


def concatenate_airr(airr_files: Iterable[str], output_file: str) -> str:
    """
    Concatenates multiple AIRR-formatted files into a single file.

    Parameters
    ----------
    airr_files : Iterable[str]
        Iterable of file paths to the AIRR-formatted files to be concatenated.

    output_file : str
        Path to the output file.

    Returns
    -------
    str
        Path to the output file.

    """
    return concatenate_csv(airr_files, output_file, separator="\t")


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


# ===============================
#     Sequence and Pair I/O
# ===============================


def to_fasta(
    sequences: Union[str, Iterable],
    fasta_file: Optional[str] = None,
    id_key: Optional[str] = None,
    sequence_key: Optional[str] = None,
    tempfile_dir: Optional[str] = None,
    append_chain: bool = True,
    as_string: bool = False,
) -> str:
    """
    Writes sequences to a FASTA-formatted file or returns a FASTA-formatted string.

    Parameters
    ----------
    sequences : str or Iterable
        Accepts any of the following:
            1. list of abutils ``Sequence`` and/or ``Pair`` objects
            2. FASTA/Q-formatted string
            3. path to a FASTA/Q-formatted file
            4. list of BioPython ``SeqRecord`` objects
            5. list of lists/tuples, of the format ``[sequence_id, sequence]``
        Required.

        .. note::
            Processing a list containing a mixture of ``Sequence`` and/or ``Pair`` objects is supported.

    fasta_file : str, default=None
        Path to the output FASTA file. If neither `fasta_file` nor `tempfile_dir`
        are provided, a FASTA-formatted string will be returned.

    id_key : str, default=None
        Name of the annotation field containing the sequence ID. If not provided,
        ``sequence.id`` is used.

    sequence_key : str, default=None
        Name of the annotation field containg the sequence. If not provided,
        ``sequence.sequence`` is used.

    tempfile_dir : str, optional
        If `fasta_file` is not provided, directory into which the tempfile
        should be created. If the directory does not exist, it will be
        created.

    append_chain : bool, default=True
        If ``True``, the chain (heavy or light) will be appended to the sequence name:
        ``>MySequence_heavy``.

        .. note::
            This option is ignored unless a list containing ``Pair`` objects is provided.

    as_string : bool, default=False
        Deprecated. Kept for backwards compatibility.


    Returns
    --------
    fasta : str
        Path to a FASTA file or a FASTA-formatted string

    """
    if isinstance(sequences, str):
        # if sequences is already a FASTA/Q-formatted file
        if os.path.isfile(sequences):
            fasta_string = "\n".join([s.fasta for s in parse_fasta(sequences)])
        # if sequences is a FASTA-formatted string
        else:
            fasta_string = sequences
    # if sequences is a list of Sequences
    elif all([isinstance(s, (Sequence, Pair)) for s in sequences]):
        fasta_strings = []
        for s in sequences:
            if isinstance(s, Pair):
                fasta_strings.append(
                    s.fasta(
                        name_field=id_key,
                        sequence_field=sequence_key,
                        append_chain=append_chain,
                    )
                )
            else:
                fasta_strings.append(
                    s.as_fasta(name_field=id_key, seq_field=sequence_key)
                )
        fasta_string = "\n".join(fasta_strings)
    # anything else..
    else:
        fasta_string = "\n".join(
            [Sequence(s, id_key=id_key, seq_key=sequence_key).fasta for s in sequences]
        )
    # output
    if fasta_file is not None:
        make_dir(os.path.dirname(fasta_file))
        with open(fasta_file, "w") as f:
            f.write(fasta_string)
        return fasta_file
    elif tempfile_dir is not None:
        make_dir(tempfile_dir)
        ff = tempfile.NamedTemporaryFile(dir=tempfile_dir, delete=False)
        ff.close()
        with open(ff.name, "w") as f:
            f.write(fasta_string)
        return ff.name
    return fasta_string


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
    annotations: Optional[Iterable[str]] = None,
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

    annotations : list, default=None
        A list of annotation fields to be included in the Polars DataFrame. For ``Sequence``
        objects, this refers to fields in the ``annotations`` field. For ``Pair`` objects,
        this refers to fields in the heavy and light chain annotations.

    columns : list, default=None
        Used only for ``Pair`` objects. A list of fields to be retained in the output Polars
        DataFrame. Fields should be column names in the input file, such as ``"sequence:0"``,
        ``"sequence:1"``, ``"name"``, etc. This option is provided to allow differential
        selection of heavy and light chain fields. For cases in which the same fields will
        be selected for both chains, it is recommended to use ``annotations`` instead.

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
            annotations=annotations,
            properties=properties,
            drop_na_columns=drop_na_columns,
            order=order,
            exclude=exclude,
            leading=leading,
        )
    elif all([isinstance(s, Pair) for s in sequences]):
        return pairs_to_polars(
            sequences,
            annotations=annotations,
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
    return from_polars(
        df, match=match, fields=fields, id_key=id_key, sequence_key=sequence_key
    )


def to_pandas(
    sequences: Iterable[Union[Sequence, Pair]],
    annotations: Optional[Iterable[str]] = None,
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

    annotations : list, default=None
        A list of annotation fields to be included in the Pandas DataFrame. For ``Sequence``
        objects, this refers to fields in the ``annotations`` field. For ``Pair`` objects,
        this refers to fields in the heavy and light chain annotations.

    columns : list, default=None
        Used only for ``Pair`` objects. A list of fields to be retained in the output Pandas
        DataFrame. Fields should be column names in the input file, such as ``"sequence:0"``,
        ``"sequence:1"``, ``"name"``, etc. This option is provided to allow differential
        selection of heavy and light chain fields. For cases in which the same fields will
        be selected for both chains, it is recommended to use ``annotations`` instead.

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
        df = sequences_to_polars(
            sequences,
            annotations=annotations,
            columns=columns,
            properties=properties,
            drop_na_columns=drop_na_columns,
            order=order,
            exclude=exclude,
            leading=leading,
        )
    elif all([isinstance(s, Pair) for s in sequences]):
        df = pairs_to_polars(
            sequences,
            annotations=annotations,
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

    return df.to_pandas()


def read_csv(
    csv_file: str,
    separator: str = ",",
    match: Optional[dict] = None,
    fields: Optional[Iterable] = None,
    id_key: str = "sequence_id",
    sequence_key: str = "sequence",
    as_dataframe: bool = False,
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

    as_dataframe : bool, default=False
        If ``True``, the function will return a ``polars.DataFrame`` object.

        .. note::
            If ``True``, ``fields`` will be used to select columns from the input file, but
            ``match`` will be ignored.

    Returns
    -------
    Iterable[Union[Sequence, Pair]] or polars.DataFrame
        A list of ``Sequence`` or ``Pair`` objects or a ``polars.DataFrame`` object.

    """
    df = pl.read_csv(csv_file, separator=separator, infer_schema_length=None)
    if as_dataframe:
        if fields is not None:
            df = df.select(fields)
        return df
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
    csv_file: str,
    separator: str = ",",
    header: bool = True,
    annotations: Optional[Iterable[str]] = None,
    columns: Optional[Iterable] = None,
    properties: Optional[Iterable[str]] = None,
    sequence_properties: Optional[Iterable[str]] = None,
    drop_na_columns: bool = True,
    order: Optional[Iterable[str]] = None,
    exclude: Optional[Union[str, Iterable[str]]] = None,
    leading: Optional[Union[str, Iterable[str]]] = None,
) -> None:
    """
    Saves a list of ``Pair`` objects to a CSV file.

    Parameters
    ----------
    sequences : Iterable[Sequence, Pair]
        List of ``Sequence`` or ``Pair`` objects to be saved to a CSV file. Required.

    csv_file : str
        Path to the output CSV file. Required.

    separator : str, default=","
        Column separator. Default is ``","``.

    header : bool, default=True
        If ``True``, the CSV file will contain a header row. Default is ``True``.

    annotations : list, default=None
        A list of annotation fields to be included in the CSV file. For ``Sequence``
        objects, this refers to fields in the ``annotations`` field. For ``Pair`` objects,
        this refers to fields in the heavy and light chain annotations.

    columns : list, default=None
        Used only for ``Pair`` objects. A list of fields to be retained in the output CSV
        file. Fields should be column names in the input file, such as ``"sequence:0"``,
        ``"sequence:1"``, ``"name"``, etc. This option is provided to allow differential
        selection of heavy and light chain fields. For cases in which the same fields will
        be selected for both chains, it is recommended to use ``annotations`` instead.

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
    None

    """
    if all([isinstance(s, Sequence) for s in sequences]):
        df = sequences_to_polars(
            sequences,
            annotations=annotations,
            columns=columns,
            properties=properties,
            drop_na_columns=drop_na_columns,
            order=order,
            exclude=exclude,
            leading=leading,
        )
    elif all([isinstance(s, Pair) for s in sequences]):
        df = pairs_to_polars(
            sequences,
            annotations=annotations,
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

    df.write_csv(csv_file, separator=separator, include_header=header)


def read_airr(
    airr_file: str,
    match: Optional[dict] = None,
    fields: Optional[Iterable] = None,
    id_key: str = "sequence_id",
    sequence_key: str = "sequence",
    as_dataframe: bool = False,
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

    as_dataframe : bool, default=False
        If ``True``, the function will return a ``polars.DataFrame`` object.

        .. note::
            If ``True``, ``fields`` will be used to select columns from the input file, but
            ``match`` will be ignored.

    Returns
    -------
    Iterable[Union[Sequence, Pair]] or polars.DataFrame
        A list of ``Sequence`` or ``Pair`` objects or a ``polars.DataFrame`` object.

    """
    return read_csv(
        airr_file,
        separator="\t",
        match=match,
        fields=fields,
        id_key=id_key,
        sequence_key=sequence_key,
        as_dataframe=as_dataframe,
    )


def to_airr(
    sequences: Iterable[Union[Sequence, Pair]],
    airr_file: str,
    annotations: Optional[Iterable[str]] = None,
    columns: Optional[Iterable] = None,
    properties: Optional[Iterable[str]] = None,
    sequence_properties: Optional[Iterable[str]] = None,
    drop_na_columns: bool = True,
    order: Optional[Iterable[str]] = None,
    exclude: Optional[Union[str, Iterable[str]]] = None,
    leading: Optional[Union[str, Iterable[str]]] = None,
) -> None:
    """
    Saves a list of ``Pair`` objects to a CSV file.

    Parameters
    ----------
    sequences : Iterable[Sequence, Pair]
        List of ``Sequence`` or ``Pair`` objects to be saved to a CSV file. Required.

    airr_file : str
        Path to the output AIRR file. Required.

    annotations : list, default=None
        A list of annotation fields to be included in the CSV file. For ``Sequence``
        objects, this refers to fields in the ``annotations`` field. For ``Pair`` objects,
        this refers to fields in the heavy and light chain annotations.

    columns : list, default=None
        Used only for ``Pair`` objects. A list of fields to be retained in the output CSV
        file. Fields should be column names in the input file, such as ``"sequence:0"``,
        ``"sequence:1"``, ``"name"``, etc. This option is provided to allow differential
        selection of heavy and light chain fields. For cases in which the same fields will
        be selected for both chains, it is recommended to use ``annotations`` instead.

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
    None

    """
    to_csv(
        sequences,
        airr_file,
        separator="\t",
        header=True,
        annotations=annotations,
        columns=columns,
        properties=properties,
        sequence_properties=sequence_properties,
        drop_na_columns=drop_na_columns,
        order=order,
        exclude=exclude,
        leading=leading,
    )


def read_parquet(
    parquet_file: str,
    match: Optional[dict] = None,
    fields: Optional[Iterable] = None,
    id_key: str = "sequence_id",
    sequence_key: str = "sequence",
    as_dataframe: bool = False,
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

    as_dataframe : bool, default=False
        If ``True``, the function will return a ``polars.DataFrame`` object.

        .. note::
            If ``True``, ``fields`` will be used to select columns from the input file, but
            ``match`` will be ignored.

    Returns
    -------
    Iterable[Union[Sequence, Pair]] or polars.DataFrame
        A list of ``Sequence`` or ``Pair`` objects or a ``polars.DataFrame`` object.

    """
    df = pl.read_parquet(parquet_file)
    if as_dataframe:
        if fields is not None:
            df = df.select(fields)
        return df
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
    parquet_file: str,
    annotations: Optional[Iterable[str]] = None,
    columns: Optional[Iterable] = None,
    properties: Optional[Iterable[str]] = None,
    sequence_properties: Optional[Iterable[str]] = None,
    drop_na_columns: bool = True,
    order: Optional[Iterable[str]] = None,
    exclude: Optional[Union[str, Iterable[str]]] = None,
    leading: Optional[Union[str, Iterable[str]]] = None,
) -> None:
    """
    Saves a list of ``Pair`` objects to a Parquet file.

    Parameters
    ----------
    sequences : Iterable[Sequence, Pair]
        List of ``Sequence`` or ``Pair`` objects to be saved to a Parquet file. Required.

    parquet_file : str
        Path to the output Parquet file. Required.

    annotations : list, default=None
        A list of annotation fields to be included in the Parquet file. For ``Sequence``
        objects, this refers to fields in the ``annotations`` field. For ``Pair`` objects,
        this refers to fields in the heavy and light chain annotations.

    columns : list, default=None
        Used only for ``Pair`` objects. A list of fields to be retained in the output Polars
        DataFrame. Fields should be column names in the input file, such as ``"sequence:0"``,
        ``"sequence:1"``, ``"name"``, etc. This option is provided to allow differential
        selection of heavy and light chain fields. For cases in which the same fields will
        be selected for both chains, it is recommended to use ``annotations`` instead.

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
    None

    """
    if all([isinstance(s, Sequence) for s in sequences]):
        df = sequences_to_polars(
            sequences,
            annotations=annotations,
            columns=columns,
            properties=properties,
            drop_na_columns=drop_na_columns,
            order=order,
            exclude=exclude,
            leading=leading,
        )
    elif all([isinstance(s, Pair) for s in sequences]):
        df = pairs_to_polars(
            sequences,
            annotations=annotations,
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

    df.write_parquet(parquet_file)


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
