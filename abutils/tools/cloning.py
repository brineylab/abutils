#!/usr/bin/env python
# filename: cloning.py
#
#
# Copyright (c) 2023 Bryan Briney
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
#


from typing import Callable, Iterable, Optional, Union

from tqdm.auto import tqdm

from ..core.sequence import Sequence, codon_optimize

__all__ = ["build_synthesis_constructs"]


def build_synthesis_constructs(
    sequences: Union[Sequence, Iterable[Sequence], str],
    overhang_5: Optional[dict] = None,
    overhang_3: Optional[dict] = None,
    sequence_key: Optional[str] = None,
    id_key: Optional[str] = None,
    locus_key: Optional[str] = None,
    frame: Optional[int] = None,
    add_locus_to_name: bool = True,
    group_by_chain: bool = False,
    show_progress: bool = True,
    sort_func: Optional[Callable] = None,
):
    """
    Builds codon-optimized synthesis constructs, including Gibson overhangs suitable
    for cloning IGH, IGK and IGL variable region constructs into antibody expression
    vectors.


    .. seealso::
        | Thomas Tiller, Eric Meffre, Sergey Yurasov, Makoto Tsuiji, Michel C Nussenzweig, Hedda Wardemann
        | Efficient generation of monoclonal antibodies from single human B cells by single cell RT-PCR and expression vector cloning
        | *Journal of Immunological Methods* 2008, doi: 10.1016/j.jim.2007.09.017


    Parameters
    ----------
    sequences : ``Sequence``, str, or ``Iterable[Sequence]``
        A ``Sequence`` object, a sequence as a ``str``, or an ``Iterable`` of ``Sequence`` objects.

        .. note::
            If the provided sequences are nucleotide sequences, ``frame`` must also be provided
            so that codon optimization can be performed on the correct reading frame.

    overhang_5 : dict, optional
        A ``dict`` mapping the locus name to 5' Gibson overhangs. By default, Gibson
        overhangs corresponding to the expression vectors in Tiller et al, 2008:

            | **IGH:** ``catcctttttctagtagcaactgcaaccggtgtacac``
            | **IGK:** ``atcctttttctagtagcaactgcaaccggtgtacac``
            | **IGL:** ``atcctttttctagtagcaactgcaaccggtgtacac``

        To produce constructs without 5' Gibson overhangs, provide an empty dictionary.

    overhang_3 : dict, optional
        A ``dict`` mapping the locus name to 3' Gibson overhangs. By default, Gibson
        overhangs corresponding to the expression vectors in Tiller et al, 2008:

            | **IGH:** ``gcgtcgaccaagggcccatcggtcttcc``
            | **IGK:** ``cgtacggtggctgcaccatctgtcttcatc``
            | **IGL:** ``ggtcagcccaaggctgccccctcggtcactctgttcccgccctcgagtgaggagcttcaagccaacaaggcc``

        To produce constructs without 3' Gibson overhangs, provide an empty dictionary.

    sequence_key : str, default='sequence_aa'
        Field containing the sequence to be codon optimized. Default is ``'sequence_aa'`` if
        ``annotation_format == 'airr'`` or ``'vdj_aa'`` if ``annotation_format == 'json'``.
        Either nucleotide or amino acid sequences are acceptable.

    id_key : str, optional
        Field containing the name of the sequence. If not provided, ``sequence.id`` is used.

    locus_key : str, optional
        Field containing the name of the locus. If not provided, ``sequence["locus"]`` is used, which
        conforms with AIRR naming conventions.

        .. note::
            The data in ``locus_key`` should return one of the following:
                * IGH, IGK, or IGL (which are the names of the loci in AIRR naming conventions)
                * heavy, kappa, or lambda

    frame : int, default=1
        Reading frame to translate. Default is ``1``.

        .. note::
            This parameter is ignored if the input sequences are amino acid sequences.

    add_locus_to_name : bool, default=True
        If ``True``, the name of the sequence will be appended with the name of the locus, separated by an underscore.
        If ``False``, the name of the sequence will not be modified.

    group_by_chain : bool, default=False
        If ``True``, the output ``list`` of ``Sequence`` objects will be separated by chain (heavy first, light second).
        Within each chain group, sequence ordering will retained (either input order, or sorted by ``sort_func``).
        If ``False``, the output ``list`` of ``Sequence`` objects will not be grouped.

    in_place : bool, default=False
        If ``True``, the input ``Sequence`` objects will be returned with the optimized sequence
        populating the ``sequence.codon_optimized`` property. If ``False``, the optimized sequences
        will be returned as new ``Sequence`` objects.

    sort_func : Callable, optional
        A function to sort the output ``list`` of ``Sequence`` objects. By default,
        the output ``list`` is is in the same order as the input.


    Returns
    -------
    sequences : ``list`` of ``Sequence`` objects
        A ``list`` of ``Sequence`` objects. Each ``Sequence`` object has the following properties:

            | *id*: The sequence ID (either ``sequence.id`` or ``id_key``), to which the locus (``locus_key``) is appended if ``add_locus_to_name`` is ``True``.
            | *sequence*: The codon-optimized sequence, including Gibson overhangs.

        If ``sort_func`` is provided, the output ``list``  will be sorted by ``sort_func``, otherwise the output is in the same order as the input.

    """
    # get overhangs
    overhang_3 = overhang_3 if overhang_3 is not None else GIBSON3
    overhang_5 = overhang_5 if overhang_5 is not None else GIBSON5

    # process input
    if isinstance(sequences, (Sequence, str)):
        sequences = [sequences]
    sequences = [Sequence(s) for s in sequences]

    # sort
    if sort_func is not None:
        sequences = sorted(sequences, key=sort_func)
    if group_by_chain:
        locus_key = locus_key if locus_key is not None else "locus"
        heavies = [s for s in sequences if s[locus_key] in ["heavy", "IGH"]]
        lights = [s for s in sequences if s[locus_key] not in ["heavy", "IGH"]]
        sequences = heavies + lights

    # make codon-optimized synthesis constructs
    if show_progress:
        sequences = tqdm(sequences)
    optimized = []
    for sequence in sequences:
        seq = sequence[sequence_key] if sequence_key is not None else sequence.sequence
        name = sequence[id_key] if id_key is not None else sequence.id
        locus = sequence[locus_key] if locus_key is not None else sequence["locus"]
        if add_locus_to_name and locus is not None:
            name += f"_{locus}"
        opt_seq = codon_optimize(seq, frame=frame, as_string=True)
        synthesis_construct = (
            overhang_5.get(locus, "") + opt_seq + overhang_3.get(locus, "")
        )
        optimized.append(Sequence(synthesis_construct, id=name))
    return optimized


GIBSON5 = {
    "IGH": "catcctttttctagtagcaactgcaaccggtgtacac",
    "IGK": "atcctttttctagtagcaactgcaaccggtgtacac",
    "IGL": "atcctttttctagtagcaactgcaaccggtgtacac",
    "heavy": "catcctttttctagtagcaactgcaaccggtgtacac",
    "kappa": "atcctttttctagtagcaactgcaaccggtgtacac",
    "lambda": "atcctttttctagtagcaactgcaaccggtgtacac",
}

GIBSON3 = {
    "IGH": "gcgtcgaccaagggcccatcggtcttcc",
    "IGK": "cgtacggtggctgcaccatctgtcttcatc",
    "IGL": "ggtcagcccaaggctgccccctcggtcactctgttcccgccctcgagtgaggagcttcaagccaacaaggcc",
    "heavy": "gcgtcgaccaagggcccatcggtcttcc",
    "kappa": "cgtacggtggctgcaccatctgtcttcatc",
    "lambda": "ggtcagcccaaggctgccccctcggtcactctgttcccgccctcgagtgaggagcttcaagccaacaaggcc",
}


LOCUS_MAP = {
    "IGH": "heavy",
    "IGK": "kappa",
    "IGL": "lambda",
    "TRA": "alpha",
    "TRB": "beta",
    "TRD": "delta",
    "TRG": "gamma",
}
