#!/usr/bin/env python
# filename: phylo.py


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

from ..bin import get_path as get_binary_path
import subprocess as sp


def numbering(ab, debug: bool = False):
    """
    Numbers all AA residues in the Sequence object with different numbering scheme available.

    Parameters
    ----------
    ab : Sequence
        Sequence object processed with AbStar (the presence of the sequence_aa field 
        is mandatory). Required.

    debug : bool, default=False
        If ``True``, verbose output is printed. Default is False.

    Returns
    -------
    ab: Sequence
        Sequence object with added annotations. Numbering can be found in ab['numbering']

    """

    ab['numbering'] = {}
    ab['numbering']['kabat'] = anarci_wrap(ab, 'kabat', debug)
    ab['numbering']['IMGT'] = anarci_wrap(ab, 'IMGT', debug)
    ab['numbering']['chothia'] = anarci_wrap(ab, 'chothia', debug)
    ab['numbering']['martin'] = anarci_wrap(ab, 'Martin', debug)
    ab['numbering']['Aho'] = anarci_wrap(ab, 'Aho', debug)
    ab['numbering']['wolfguy'] = anarci_wrap(ab, 'Wolfguy', debug)
    return ab


def anarci_wrap(ab,
                numbering_scheme: str ='IMGT',
                anarci_path: str = None,
                hmmer_path: str = None,
                debug: bool = False):
    if hmmer_path == None:
        hmmer = get_binary_path("hmmer")
    else:
        hmmer = hmmer_path
    if anarci_path == None:
        anarci = get_binary_path("hmmer")
    else:
        anarci = anarci_path
    scheme_dict = {'IMGT':'i',
                   'kabat':'k',
                   'chothia':'c',
                   'Martin':'m',
                   'Aho':'a',
                   'Wolfguy':'w',
                   }
    cmd = f"{anarci} -i {ab['sequence_aa']} --scheme {scheme_dict[numbering_scheme]} --hmmerpath {hmmer}"
    anarci_cmd = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True, universal_newlines=True)
    stdout, stderr = anarci_cmd.communicate()
    if debug:
        return stdout
    else:
        raw = stdout.splitlines()[7:-1]
        numbering = {}
        try:
            for i, e in enumerate(raw):
                elems = e.split(' ')
                numbering[i] = [elems[1], elems[-1]]
        except TypeError:
            numbering = "Sorry :-/ Numbering scheme cannot be assessed"
        return numbering
