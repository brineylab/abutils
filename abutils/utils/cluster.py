#!/usr/bin/env python
# filename: cluster.py


#
# Copyright (c) 2015 Bryan Briney
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


from __future__ import absolute_import, division, print_function, unicode_literals

from datetime import datetime
import multiprocessing as mp
import os
import platform
import random
import shutil
import sqlite3
import string
import subprocess as sp
import sys
import tempfile
import time

from Bio import SeqIO, AlignIO
from Bio.Align import AlignInfo

from .alignment import mafft
from .database import KeyValueStore
from .decorators import lazy_property
from .pipeline import make_dir
from ..core.sequence import Sequence, read_fasta


class Cluster():
    '''
    Docstring for Cluster.
    '''
    def __init__(self, name, sequences, centroid=None):
        self.name = name
        self.sequences = sequences
        self.centroid = centroid
        
    @property
    def size(self):
        return len(self.sequences)
    
    @property
    def seq_ids(self):
        return [s.id for s in self.sequences]
    
    @lazy_property
    def consensus(self):
        return self._make_consensus()

    def _make_consensus(self):
        if len(self.sequences) == 1:
            return self.sequences[0]
        aln = mafft(self.sequences)
        if aln is None:
            print("ERROR: Failed to generate an alignmnet for a consensus sequence.")
            return None
        summary_align = AlignInfo.SummaryInfo(aln)
        consensus = summary_align.gap_consensus(threshold=0.51, ambiguous='n')
        consensus_string = str(consensus).replace('-', '')
        consensus_seq = Sequence(consensus_string.upper())
        return consensus_seq


class Clusters():
    '''
    Docstring for Clusters.
    '''
    def __init__(self, clusters=None):
        self._clusters = clusters if clusters is not None else []
        
    def __iter__(self):
        for cluster in self.clusters:
            yield cluster
            
    def __getitem__(self, key):
        return self.clusters[key]
    
    def __len__(self):
        return len(self._clusters)
    
    @property
    def clusters(self):
        return sorted(self._clusters, key=lambda x: x.size, reverse=True)
    
    @property
    def largest_cluster(self):
        return self.clusters[0]
            
    @property
    def count(self):
        return len(self.clusters)
    
    def add(self, cluster):
        self._clusters.append(cluster)



def cluster(sequences, threshold=0.975, temp_dir='/tmp', in_memory_db=True, debug=False):
    # preflight
    if isinstance(sequences, Sequence):
        sequences = [sequences, ]
    elif isinstance(sequences, (list, tuple)):
        pass
    elif os.path.isfile(sequences):
        sequences = read_fasta(sequences)
    ts = datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f")
    temp_dir = os.path.join(temp_dir, 'clustering_{}'.format(ts))
    make_dir(temp_dir)
    
    # set up a sequence db
    if in_memory_db:
        db = KeyValueStore(in_memory=True)
    else:
        db = KeyValueStore(name='mmseqs', direc=temp_dir)
    for s in sequences:
        db[s.id] = s
    db.index()
    
    # make the input FASTA file
    fasta_file = os.path.join(temp_dir, 'input.fasta')
    with open(fasta_file, 'w') as f:
        fstring = [s.fasta for s in sequences]
        f.write('\n'.join(fstring))
    
    # cluster
    clu_dict = mmseqs_cluster(fasta_file, min_seq_id=threshold, temp_dir=temp_dir, debug=debug)
    
    # parse clustering output
    clus = []
    for name, seqs in clu_dict.items():
        seq_ids = seqs['seq_ids']
        centroid_id = seqs['centroid']
        seqs = db.find_many(seq_ids)
        centroid = db.find_one(centroid_id)
        clus.append(Cluster(name, seqs, centroid=centroid))
    clusters = Clusters(clus)
    
    # clean up
    if not debug:
        db.destroy()
        shutil.rmtree(temp_dir)
    
    return clusters




def mmseqs_cluster(fasta_file, min_seq_id=0.975, temp_dir='/tmp', debug=False):
    # preflight
    fasta_file = os.path.abspath(fasta_file)
    temp_dir = os.path.join(temp_dir, 'mmseqs')
    make_dir(temp_dir)
    mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    mmseqs_bin = os.path.join(mod_dir, 'bin/mmseqs_{}'.format(platform.system().lower()))

    # create the sequence DB
    db_file = os.path.join(temp_dir, 'DB')
    db_cmd = '{} createdb {} {}'.format(mmseqs_bin, fasta_file, db_file)
    p = sp.Popen(db_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    stdout, stderr = p.communicate()
    if debug:
        if sys.version_info[0] > 2:
            stdout = stdout.decode('utf-8')
            stderr = stderr.decode('utf-8')
        print('STDOUT:', stdout)
        print('')
        print('STDERR:', stderr)

    # do the clustering
    clu_file = os.path.join(temp_dir, 'CLU')
    cluster_cmd = '{} cluster {} {} {} --min-seq-id {}'.format(mmseqs_bin, db_file, clu_file, temp_dir, min_seq_id)
    p = sp.Popen(cluster_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    stdout, stderr = p.communicate()
    if debug:
        if sys.version_info[0] > 2:
            stdout = stdout.decode('utf-8')
            stderr = stderr.decode('utf-8')
        print('STDOUT:', stdout)
        print('')
        print('STDERR:', stderr)

    # generate TSV-formatted output
    tsv_file = os.path.join(temp_dir, 'CLU.tsv')
    tsv_cmd = '{0} createtsv {1} {1} {2} {3}'.format(mmseqs_bin, db_file, clu_file, tsv_file)
    p = sp.Popen(tsv_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    stdout, stderr = p.communicate()
    if debug:
        if sys.version_info[0] > 2:
            stdout = stdout.decode('utf-8')
            stderr = stderr.decode('utf-8')
        print('STDOUT:', stdout)
        print('')
        print('STDERR:', stderr)

    # parse TSV output
    clu_names = {}
    clu_seqs = {}
    with open(tsv_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                continue
            c, s = line.split()
            if c not in clu_names:
                cname = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
                clu_names[c] = cname
            clu_name = clu_names[c]
            if clu_name not in clu_seqs:
                clu_seqs[clu_name] = {'centroid': c, 'seq_ids': []}
            clu_seqs[clu_name]['seq_ids'].append(s)

    # clean up
    if not debug:
        shutil.rmtree(temp_dir)
    return clu_seqs




# __all__ = ['CDHITResult', 'Cluster', 'cluster', 'cdhit', 'parse_clusters']


# class CDHITResult(object):
#     '''
#     Docstring for CDHITResult.
#     '''
#     def __init__(self, clusters, seq_db=None, db_path=None, seq_dict=None):
#         self.clusters = sorted(clusters, key = lambda x: x.size, reverse=True)
#         self.seq_db = seq_db
#         self.db_path = db_path
#         self.seq_dict = seq_dict
    
#     def __len__(self):
#         return len(self.clusters)

#     def __iter__(self):
#         for cluster in self.clusters:
#             yield cluster


#     @property
#     def largest_cluster(self):
#         return self.clusters[0]

    
#     def delete(self):
#         self.seq_db.connection.close()
#         os.unlink(self.db_path)


# class Cluster(object):
#     """
#     Data and methods for a cluster of sequences.

#     All public attributes are evaluated lazily, so attributes that
#     require significant processing time are only computed when needed.
#     In addition, attributes are cached after calculation, so if you
#     change the Cluster object after accessing attributes, the
#     attributes will not automatically update. Setters are provided for all attributes,
#     however, so you can update them manually if necessary::

#         seqs = [Sequence1, Sequence2, ... SequenceN]
#         clust = cluster(seqs)

#         # calculate the consensus
#         consensus = clust.consensus

#         # add sequences to the Cluster
#         more_sequences = [SequenceA, SequenceB, SequenceC]
#         clust.sequences += more_sequences

#         # need to recompute the consensus manually
#         clust.consensus = clust._make_consensus()


#     Attributes:

#         ids (list): A list of all sequence IDs in the Cluster

#         size (int): Number of sequences in the Cluster

#         sequences (list): A list of all sequences in the Cluster,
#             as AbTools ``Sequence`` objects.

#         consensus (Sequence): Consensus sequence, calculated by
#             aligning all sequences with MAFFT and computing the
#             ``Bio.Align.AlignInfo.SummaryInfo.gap_consensus()``

#         centroid (Sequence): Centroid sequence, as calculated by
#             CD-HIT.
#     """
#     def __init__(self, raw_cluster, seq_db=None, db_path=None, seq_dict=None):
#         super(Cluster, self).__init__()
#         self._raw_cluster = raw_cluster
#         self._seq_db = seq_db
#         self._seq_db_path = db_path
#         self._seq_dict = seq_dict


#     @lazy_property
#     def ids(self):
#         return self._get_ids()

#     @lazy_property
#     def size(self):
#         return len(self.ids)

#     @lazy_property
#     def sequences(self):
#         if all([self._seq_db is None, self._seq_dict is None]):
#             err = "ERROR: In order to access a Cluster's sequences, you must provide "
#             err += 'either a SQLite database connection object or a dictionary of sequences at instantiation.'
#             raise RuntimeError(err)
#         return self._get_sequences()

#     @lazy_property
#     def consensus(self):
#         if all([self._seq_db is None, self._seq_dict is None]):
#             err = "ERROR: In order to compute a Cluster's consensus, you must provide "
#             err += 'either a SQLite database connection object or a dictionary of sequences at instantiation.'
#             raise RuntimeError(err)
#         return self._make_consensus()

#     @lazy_property
#     def centroid(self):
#         if all([self._seq_db is None, self._seq_dict is None]):
#             err = "ERROR: In order to retrieve a Cluster's centroid, you must provide "
#             err += 'either a SQLite database connection object or a dictionary of sequences at instantiation.'
#             raise RuntimeError(err)
#         return self._get_centroid()


#     def cleanup(self):
#         try:
#             os.unlink(self._raw_cluster)
#         except:
#             pass
#         self.terminate_db()


#     def terminate_db(self):
#         if self._seq_db is not None:
#             self._seq_db.close()
#         if self._seq_db_path is not None:
#             if os.path.isfile(self._seq_db_path):
#                 os.unlink(self._seq_db_path)


#     def _get_ids(self):
#         ids = []
#         for c in self._raw_cluster[1:]:
#             if c:
#                 ids.append(c.split()[2][1:-3])
#         return ids

#     def _get_sequences(self):
#         if self._seq_db is not None:
#             seqs = []
#             for chunk in self._chunker(self.ids):
#                 sql_cmd = '''SELECT seqs.id, seqs.sequence
#                              FROM seqs
#                              WHERE seqs.id IN ({})'''.format(','.join('?' * len(chunk)))
#                 seq_chunk = self._seq_db.execute(sql_cmd, chunk)
#                 seqs.extend(seq_chunk)
#             return [Sequence(s) for s in seqs]
#         else:
#             return [self._seq_dict[s] for s in self.ids]

#     def _get_centroid(self):
#         for line in self._raw_cluster:
#            if '*' in line:
#                 centroid_id = line.split()[2][1:-3]
#                 break
#         centroid = [s for s in self.sequences if s.id == centroid_id][0]
#         return centroid

#     def _make_consensus(self):
#         if len(self.sequences) == 1:
#             return self.sequences[0]
#         _aln = mafft(self.sequences, as_file=True)
#         aln = AlignIO.read(open(_aln, 'r'), 'fasta')
#         summary_align = AlignInfo.SummaryInfo(aln)
#         consensus = summary_align.gap_consensus(threshold=0.51, ambiguous='n')
#         consensus_string = str(consensus).replace('-', '')
#         consensus_seq = Sequence(consensus_string.upper())
#         os.unlink(_aln)
#         return consensus_seq

#     @staticmethod
#     def _chunker(l, size=900):
#         return (l[pos:pos + size] for pos in range(0, len(l), size))


# def cluster(seqs, threshold=0.975, out_file=None, temp_dir=None, make_db=True, method='cdhit',
#             quiet=False, threads=0, return_just_seq_ids=False, max_memory=800, retries=5, debug=False):
#     '''
#     Perform sequence clustering with CD-HIT.

#     Args:

#         seqs (list or filepath): An iterable of sequences, in any format that ``abutils.utils.sequence.Sequence``
#             can handle. Alternatively, ``seqs`` can be the path to a FASTA-formatted file of sequences. If a path
#             is provided, ``make_db`` is set to ``False`` and the input file will not be deleted following clustering.

#         threshold (float): Clustering identity threshold. Default is ``0.975``.

#         out_file (str): Path to the clustering output file. Default is to use
#             ``tempfile.NamedTemporaryFile`` to generate an output file name.

#         temp_dir (str): Path to the temporary directory. If not provided, ``'/tmp'`` is used.

#         make_db (bool): Whether to build a SQlite database of sequence information. Required
#             if you want to calculate consensus/centroid sequences for the resulting
#             clusters or if you need to access the clustered sequences (not just sequence IDs)
#             Default is ``True``.

#         method (str): Doesn't do anything at this point, stubbed to allow expansion of sequencing tools
#             beyond CD-HIT (UCLUST, etc).
        
#         quiet (bool): If ``True``, surpresses printing of output/progress info. Default is ``False``.

#         threads (int): Number of threads (CPU cores) to be used for clustering. Default is ``0``,
#             which results in all available cores being used.
        
#         return_just_seq_ids (bool): If ``True``, will return a 2D list of sequence IDs
#             (a list containing a list of sequence IDs for each cluster) rather than returning a
#             list of ``Cluster`` objects.

#         max_memory (int): Max memory (in MB) for CD-HIT. Will be passed directly to CD-HIT through
#             the ``-M`` runtime option. Default is ``800``.
        
#         retries (int): Occasionally (and for no discernable reason), the CD-HIT output files are empty. 
#             If this occurs, ``cdhit`` will attempt to retry running CD-HIT until the output files are 
#             not empty. Default is 5.

#         debug (bool): If ``True``, print standard output and standard error from CD-HIT. Default is ``False``.

#     Returns:

#         list: A ``CDHITResult`` object or, if ``return_just_seq_ids`` is ``True``, a 2D list of sequence IDs.
#     '''
#     if make_db:
#         ofile, cfile, seq_db, db_path = cdhit(seqs, out_file=out_file, temp_dir=temp_dir,
#                                               threshold=threshold, make_db=True, quiet=quiet,
#                                               threads=threads, max_memory=max_memory, debug=debug)            
#         return parse_clusters(ofile, cfile, seq_db=seq_db, db_path=db_path, return_just_seq_ids=return_just_seq_ids)
#     else:
#         if isinstance(seqs, (list, tuple)):
#             seqs = [Sequence(s) for s in seqs]
#             seq_dict = {s.id: s for s in seqs}
#         elif os.path.exists(seqs):
#             seq_dict = None
#         ofile, cfile, = cdhit(seqs, out_file=out_file, temp_dir=temp_dir, threads=threads,
#                               threshold=threshold, make_db=False, quiet=quiet,
#                               max_memory=max_memory, debug=debug)
#         return parse_clusters(ofile, cfile, seq_dict=seq_dict, return_just_seq_ids=return_just_seq_ids)


# def cdhit(seqs, out_file=None, temp_dir=None, threshold=0.975, make_db=True,
#           quiet=False, threads=0, max_memory=800, retries=5, debug=False):
#     '''
#     Run CD-HIT.

#     Args:

#         seqs (list or filepath): An iterable of sequences, in any format that `abutils.utils.sequence.Sequence()`
#             can handle. Alternatively, ``seqs`` can be the path to a FASTA-formatted file of sequences. If a file path
#             is provided, ``make_db`` is set to ``False`` and the input file will not be deleted following clustering.

#         threshold (float): Clustering identity threshold. Default is `0.975`.

#         out_file (str): Path to the clustering output file. Default is to use
#             `tempfile.NamedTemporaryFile` to generate an output file name.

#         temp_dir (str): Path to the temporary directory. If not provided, the default 
#             temporary directory is used (typically ``/tmp`` on MacOS and Linux).

#         make_db (bool): Whether to build a SQlite database of sequence information. Required
#             if you want to calculate consensus/centroid sequences for the resulting
#             clusters or if you need to access the clustered sequences (not just sequence IDs). 
#             Default is `True`.
        
#         quiet (bool): If `True`, surpresses printing of output/progress info. Default is `False`.

#         threads (int): Number of threads (CPU cores) to be used for clustering. Default is `0`,
#             which results in all available cores being used.

#         max_memory (int): Max memory (in MB) for CD-HIT. Will be passed directly to CD-HIT through
#             the `-M` runtime option. Default is `800`.

#         retries (int): Occasionally (and for no discernable reason), the CD-HIT output files are empty. 
#             If this occurs, ``cdhit`` will attempt to retry running CD-HIT until the output files are 
#             not empty. Default is 5.

#         debug (bool): If `True`, print standard output and standard error from CD-HIT. Default is `False`.

#     Returns:

#         If `make_db` is `True`, returns the CD-HIT output file path, the CD-HIT cluster file path,
#             a `sqlite3` database connection object, and the database path. If `make_db` is `False`, only the
#             CD-HIT output file path and CD-HIT cluster file path are returned.
#     '''
#     start_time = time.time()
#     delete_input = False if debug else True
#     if not quiet:
#         print('CD-HIT: clustering {} seqeunces'.format(len(seqs)))
#     if out_file is None:
#         out_file = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False)
#         out_file.close()
#         ofile = out_file.name
#     else:
#         ofile = os.path.expanduser(out_file)
#     cfile = ofile + '.clstr'
#     with open(ofile, 'w') as f: f.write('')
#     with open(cfile, 'w') as f: f.write('')
#     if isinstance(seqs, (list, tuple)):
#         seqs = [Sequence(s) for s in seqs]
#         ifile = _make_cdhit_input(seqs, temp_dir)
#     elif os.path.exists(seqs):
#         ifile = seqs
#         delete_input = False
#     word_size = _get_cdhit_wordsize(threshold)
#     cdhit_cmd = 'cdhit -i {} -o {} -c {} -n {} -d 0 -T {} -M {}'.format(ifile,
#                                                                         ofile,
#                                                                         threshold,
#                                                                         word_size,
#                                                                         threads,
#                                                                         max_memory)
#     while not all([os.path.getsize(cfile), os.path.getsize(cfile)]):
#         cluster = sp.Popen(cdhit_cmd,
#                         shell=True,
#                         stdout=sp.PIPE,
#                         stderr=sp.PIPE)
#         stdout, stderr = cluster.communicate()
#         if not retries:
#             break
#         retries -= 1
#     end_time = time.time()
#     if debug:
#         print(stdout)
#         print(stderr)
#     if delete_input:
#         os.unlink(ifile)
#     if not quiet:
#         print('CD-HIT: clustering took {:.2f} seconds'.format(end_time - start_time))
#     if make_db:
#         if not quiet:
#             print('CD-HIT: building a SQLite3 database')
#         seq_db, db_path = _build_seq_db(seqs, direc=temp_dir)
#         return ofile, cfile, seq_db, db_path
#     return ofile, cfile


# def parse_clusters(out_file, clust_file, seq_db=None, db_path=None, seq_dict=None, return_just_seq_ids=False):
#     '''
#     Parses CD-HIT output.

#     Args:

#         out_file (str): Path to the CD-HIT output file. Required.

#         clust_file (str): Path to the CD-HIT cluster file. Required.

#         seq_db (sqlite.Connection): SQLite3 `Connection` object. Default is ``None``. If not provided and
#             ``return_just_seq_ids`` is False, the returned ``Cluster`` objects will not contain any sequence
#             information beyond the sequence ID.

#         db_path (str): Path to a SQLite3 database file. Default is ``None``. Must be provided if
#             ``seq_db`` is also provided.

#         seq_dict (dict): A ``dict`` mapping sequence IDs to ``abutils.core.sequence.Sequence`` objects. Default
#             is ``None``. Typically used when a relatively small number of sequences are being clustered and
#             creating a ``sqlite3`` database would be overkill.

#         return_just_seq_ids (bool): If ``True``, will return a 2D list of sequence IDs
#             (a list containing a list of sequence IDs for each cluster) rather than returning a
#             list of ``Cluster`` objects.

#     Returns:

#         A ``CDHITResult`` object or, if ``return_just_seq_ids`` is ``True``, a 2D list of sequence IDs.
#     '''
#     raw_clusters = [c.split('\n') for c in open(clust_file, 'r').read().split('\n>')]
#     if return_just_seq_ids:
#         ids = []
#         for rc in raw_clusters:
#             _ids = []
#             for c in rc[1:]:
#                 if c:
#                     _ids.append(c.split()[2][1:-3])
#             ids.append(_ids)
#         os.unlink(out_file)
#         os.unlink(clust_file)
#         return ids
#     os.unlink(out_file)
#     os.unlink(clust_file)
#     clusters = [Cluster(rc, seq_db, db_path, seq_dict) for rc in raw_clusters]
#     return CDHITResult(clusters, seq_db=seq_db, db_path=db_path, seq_dict=seq_dict)


# def _make_cdhit_input(seqs, temp_dir):
#     ifile = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False)
#     ifile.close()
#     fastas = [s.fasta for s in seqs]
#     with open(ifile.name, 'w') as f:
#         f.write('\n'.join(fastas))
#     return ifile.name


# def _get_cdhit_wordsize(threshold):
#     if threshold > 0.7:
#         return 5
#     elif threshold > 0.6:
#         return 4
#     elif threshold > 0.5:
#         return 3
#     return 2


# def _build_seq_db(seqs, direc=None):
#     # '''
#     # Builds a SQLite3 database of sequences.

#     # Inputs are a list of Sequence objects and an optional directory to store the database.
#     # If ::direc:: is not provided, '/tmp' will be used.

#     # Returns a SQLite3 connection object and the database path.
#     # '''
#     direc = direc if direc is not None else '/tmp'
#     db_name = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
#     db_path = os.path.join(direc, db_name)
#     conn = sqlite3.connect(db_path)
#     c = conn.cursor()
#     create_cmd = '''CREATE TABLE seqs (id text, sequence text)'''
#     insert_cmd = 'INSERT INTO seqs VALUES (?,?)'
#     c.execute('DROP TABLE IF EXISTS seqs')
#     c.execute(create_cmd)
#     c.executemany(insert_cmd, [(str(s.id), str(s.sequence)) for s in seqs])
#     c.execute('CREATE INDEX seq_index ON seqs (id)')
#     return c, db_path
