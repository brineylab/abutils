#!/usr/bin/env python
# filename: seqio.py


#
# Copyright (c) 2018 Bryan Briney
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


from abc import ABCMeta, abstractmethod
import json
import os
import sys

from Bio import SeqIO

from . import mongodb
from .decorators import lazy_property
from .pipeline import list_files
from ..core.sequence import Sequence


if sys.version_info[0] > 2:
    STR_TYPES = [
        str,
    ]
else:
    STR_TYPES = [str, unicode]


# def read_input(input, data_type,
#                collection=None, mongo_ip='localhost', mongo_port=27017, mongo_user=None, mongo_password=None,
#                query=None, projection=None, verbose=False, **kwargs):
#     '''
#     Returns an Input class based on the data information provided.

#     Args:

#         data_type (str): One of the following: `'fasta'`, `'json'`, or `'mongodb'`.

#         input (str): Path to an input file for FASTA or JSON data types, or the database name for MongoDB data.

#         collection (str): Name of a MongoDB collection. Required for the MongoDB data type.

#         mongo_ip (str): IP address of the MongoDB server. Default is `'localhost'` if not provided.

#         mongo_port (int): Port of the MongoDB server. Default is `27017` if not provided.

#         query (dict): Query to limit the results returned from a MongoDB database.

#         projection (dict): Projection to specify fields to be retuired from a MongoDB database.
#     '''
#     if data_type.lower() == 'mongodb':
#         return MongoDBInput(database=input, collection=collection, ip=mongo_ip, port=mongo_port,
#                             user=mongo_user, password=mongo_password, query=query, projection=projection)
#     elif data_type.lower() == 'json':
#         return JSONInput(input)
#     elif data_type.lower() == 'fasta':
#         return FASTAInput(input)
#     else:
#         err = '\n\nERROR: data_type must be one of the following:\n'
#         err += 'json, mongodb, or fasta\n\n'
#         print(err)
#         sys.exit(1)


def read_fasta(fasta_file, verbose=False):
    return FASTAInput(fasta_file, verbose=verbose)


def from_fasta(fasta_file, verbose=False):
    return FASTAInput(fasta_file, verbose=verbose)


def from_json(json_file, seq_field="vdj_nt", verbose=False):
    return JSONInput(json_file, seq_field=seq_field, verbose=verbose)


def from_mongodb(
    db,
    collection=None,
    ip="localhost",
    port=27017,
    user=None,
    password=None,
    query=None,
    projection=None,
    seq_field="vdj_nt",
    verbose=False,
):
    return MongoDBInput(
        database=db,
        collection=collection,
        ip=ip,
        port=port,
        user=user,
        password=password,
        query=query,
        projection=projection,
        seq_field=seq_field,
        verbose=verbose,
    )


class BaseInput:
    """
    Base class for parsing inputs (JSON, MongoDB, FASTA files, etc)
    """

    __metaclass__ = ABCMeta

    def __init__(self):
        pass

    @property
    @abstractmethod
    def data_type(self):
        "Returns the data type"
        pass

    @property
    @abstractmethod
    def as_list(self):
        "Returns the input as a list of Sequence objects"
        pass

    @property
    @abstractmethod
    def as_generator(self):
        " Returns the input as a genarator of Sequence objects"
        pass


class FASTAInput(BaseInput):
    """
    Representation of FASTA input data.
    """

    def __init__(self, data, verbose=False):
        self.input = data
        self.verbose = verbose

    @property
    def data_type(self):
        return "fasta"

    @property
    def files(self):
        if type(self.input) in STR_TYPES:
            if os.path.isdir(self.input):
                return list_files(self.input, "json")
            else:
                return [
                    self.input,
                ]
        else:
            return self.input

    @lazy_property
    def as_list(self):
        sequences = []
        for input_file in self.files:
            if self.verbose:
                print(input_file)
            with open(input_file, "r") as f:
                for seq in SeqIO.parse(f, "fasta"):
                    sequences.append(Sequence(str(seq.seq), id=seq.id))
        return sequences

    @property
    def as_generator(self):
        for input_file in self.files:
            if self.verbose:
                print(input_file)
            with open(input_file, "r") as f:
                for seq in SeqIO.parse(f, "fasta"):
                    yield Sequence(str(seq.seq), id=seq.id)


class JSONInput(BaseInput):
    """
    Representation of JSON input data
    """

    def __init__(self, data, seq_field="vdj_nt", verbose=False):
        self.input = data
        self.seq_field = seq_field
        self.verbose = verbose

    @property
    def data_type(self):
        return "json"

    @property
    def files(self):
        if type(self.input) in STR_TYPES:
            if os.path.isdir(self.input):
                return list_files(self.input, "json")
            else:
                return [
                    self.input,
                ]
        else:
            return self.input

    @lazy_property
    def as_list(self):
        sequences = []
        for input_file in self.files:
            if self.verbose:
                print(input_file)
            with open(input_file, "r") as f:
                for line in f:
                    j = json.loads(line.strip().lstrip("]").rstrip("]").rstrip(","))
                    sequences.append(Sequence(j, seq_key=self.seq_field))
        return sequences

    @property
    def as_generator(self):
        for input_file in self.files:
            if self.verbose:
                print(input_file)
            with open(input_file, "r") as f:
                for line in f:
                    j = json.loads(
                        line.strip().lstrip("[").rstrip("]").rstrip().rstrip(",")
                    )
                    yield Sequence(j, seq_key=self.seq_field)


class MongoDBInput(BaseInput):
    """
    Representation of MongoDB input data
    """

    def __init__(
        self,
        database,
        collection,
        ip,
        port,
        user,
        password,
        query,
        projection,
        seq_field="vdj_nt",
        verbose=False,
    ):
        self.db_name = database
        self.raw_collections = collection
        self.ip = ip
        self.port = port
        self.user = user
        self.password = password
        self.query = query
        self.projection = projection
        self.seq_field = seq_field
        self.verbose = verbose

    @property
    def data_type(self):
        return "mongodb"

    @property
    def db(self):
        return mongodb.get_db(
            self.db_name,
            ip=self.ip,
            port=self.port,
            user=self.user,
            password=self.password,
        )

    @property
    def collections(self):
        if type(self.raw_collections) in STR_TYPES:
            return [
                self.raw_collections,
            ]
        elif self.raw_collections is None:
            return mongodb.get_collections(self.db)
        else:
            return self.raw_collections

    @lazy_property
    def as_list(self):
        sequences = []
        for collection in self.collections:
            if self.verbose:
                print(collection)
            res = self.db[collection].find(self.query, self.projection)
            for r in res:
                if self.seq_field not in r:
                    continue
                sequences.append(Sequence(r, seq_key=self.seq_field))
        return sequences

    @property
    def as_generator(self):
        for collection in self.collections:
            if self.verbose:
                print(collection)
            res = self.db[collection].find(self.query, self.projection)
            for r in res:
                if self.seq_field not in r:
                    continue
                yield Sequence(r, seq_key=self.seq_field)

    def _process_collections(self, collection):
        if type(collection) in STR_TYPES:
            return [
                collection,
            ]
        elif collection is None:
            return mongodb.get_collections(self.db)
        else:
            return collection

