#!/usr/local/bin/python
# filename: databases.py

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


from abc import ABCMeta, abstractmethod
import os
import sqlite3
import sys
import time

if sys.version_info[0] > 2:
    STR_TYPES = [
        str,
    ]
    import pickle
else:
    STR_TYPES = [str, unicode]
    import cPickle as pickle


class SQLiteDatabase:
    """
    base SQLiteDatabase class
    """

    __metaclass__ = ABCMeta

    def __init__(self, name=None, direc=None, in_memory=False, table_name=None):
        super(SQLiteDatabase, self).__init__()
        self.name = name
        self.initialized = False
        if all([name is not None, direc is not None]):
            self.dir = os.path.abspath(direc)
            self.path = os.path.join(self.dir, self.name)
            if os.path.isfile(self.path):
                self.initialized = True
        elif in_memory:
            self.dir = None
            self.path = ":memory:"
        self.table_name = "seqs" if table_name is None else table_name
        self._connection = None
        self._cursor = None
        self._create_table_cmd = None
        if not self.initialized:
            self.create_table()

    @property
    @abstractmethod
    def structure(self):
        pass

    @property
    def columns(self):
        return [s[0] for s in self.structure]

    @property
    def insert_cmd(self):
        return "INSERT INTO {} VALUES ({})".format(
            self.table_name, ",".join(["?"] * len(self.structure))
        )

    @property
    def connection(self):
        if self._connection is None:
            self._connection = sqlite3.connect(self.path)
        return self._connection

    @connection.setter
    def connection(self, connection):
        self._connection = connection

    @property
    def cursor(self):
        if self._cursor is None:
            self._cursor = self.connection.cursor()
        return self._cursor

    @cursor.setter
    def cursor(self, cursor):
        self._cursor = cursor

    @property
    def count(self):
        result = self.cursor.execute("SELECT COUNT(*) FROM {}".format(self.table_name))
        return result.fetchone()[0]

    @property
    def create_table_cmd(self):
        if self._create_table_cmd is None:
            field_string = ", ".join([" ".join(s) for s in self.structure])
            self._create_table_cmd = "CREATE TABLE {} ({})".format(
                self.table_name, field_string
            )
        return self._create_table_cmd

    # @abstractmethod
    # def insert_one(self, data, value=None):
    #     '''
    #     Inserts a single database entry
    #     '''
    #     pass

    # @abstractmethod
    # def insert_many(self, data):
    #     '''
    #     Inserts multiple database entries.
    #     '''
    #     pass

    # @abstractmethod
    # def find_one(self):
    #     '''
    #     Returns a single match.
    #     '''
    #     pass

    # @abstractmethod
    # def find_many(self):
    #     '''
    #     Returns multiple matches.
    #     '''
    #     pass

    # @abstractmethod
    # def find_all(self):
    #     '''
    #     Returns all values in a SQLite database.
    #     '''
    #     pass

    # @abstractmethod
    # def delete(self):
    #     '''
    #     Deletes matching database entries
    #     '''
    #     pass

    def commit(self):
        self.connection.commit()

    def close(self):
        self.connection.close()
        del self._cursor
        del self._connection
        self._cursor = None
        self._connection = None

    def create_table(self):
        self.cursor.execute("DROP TABLE IF EXISTS {}".format(self.table_name))
        self.cursor.execute(self.create_table_cmd)
        self.initialized = True

    def destroy(self):
        self.close()
        if self.path != ":memory:":
            os.unlink(self.path)

    # def insert(self, *args, **kwargs):
    #     '''
    #     Alias for ``insert_many``
    #     '''
    #     return self.insert_many(*args, **kwargs)

    # def find(self, *args, **kwargs):
    #     '''
    #     Alias for ``find_many``
    #     '''
    #     return self.find_many(*args, **kwargs)

    def index(self, fields):
        """
        Indexes the database

        Inputs:
          field - the field or list of field on which to create the index.
        """
        if type(fields) in STR_TYPES:
            fields = [fields]
        index_name = "__".join(fields) + "__index"
        self.cursor.execute(
            "CREATE INDEX {} ON {} ({})".format(
                index_name, self.table_name, ", ".join(fields)
            )
        )

    @staticmethod
    def chunker(l, n=900):
        """
        Yields successive n-sized chunks from l.
        """
        for i in range(0, len(l), n):
            yield l[i : i + n]


class KeyValueStore(SQLiteDatabase):
    """
    A simple Database with a single key/value per entry. Values are pickled upon insert and
    unpickled upon retrieval, so any pickleable object can be stored. Also provides dictionary-style
    insert and retrieval of values.
    """

    def __init__(self, name=None, direc=None, in_memory=False, table_name=None):
        super(KeyValueStore, self).__init__(
            name=name, direc=direc, in_memory=in_memory, table_name=table_name
        )

    def __getitem__(self, key):
        return self.find_one(key)

    def __setitem__(self, key, value):
        return self.insert_one(key, value)

    @property
    def structure(self):
        return [("key", "text"), ("value", "text")]

    @property
    def insert_cmd(self):
        return "INSERT INTO {} VALUES ({})".format(
            self.table_name, ",".join(["?"] * len(self.structure))
        )

    def keys(self):
        results = self.cursor.execute(
            """SELECT {0}.key
               FROM {0}""".format(
                self.table_name
            )
        )
        return [r[0] for r in results]

    def values(self):
        results = self.cursor.execute(
            """SELECT {0}.value
               FROM {0}""".format(
                self.table_name
            )
        )
        return [pickle.loads(r[0]) for r in results]

    def items(self):
        results = self.cursor.execute(
            """SELECT {0}.key, {0}.value
               FROM {0}""".format(
                self.table_name
            )
        )
        return [(r[0], pickle.loads(r[1])) for r in results]

    def insert_one(self, data, value=None):
        """
        Inserts a single key/value pair

        Inputs:
          data - data to be inserted (as a list or tuple), or optionally, a key name
          value - optional, only used if a key is passed to ``data``
        """
        if all([value is not None, type(data) in STR_TYPES]):
            data = (data, pickle.dumps(value, protocol=0))
        else:
            data = (data[0], pickle.dumps(data[1], protocol=0))
        with self.connection as conn:
            conn.execute(self.insert_cmd, data)

    def insert_many(self, data):
        """
        Inserts multiple key/value pairs.

        Inputs:
          data - a list/generator of iterables (lists or tuples) with each iterable containing
                a single key/value pair.
        """
        # for kv in data:
        #     if len(kv) != len(self.structure):
        #         err = 'Mismatch between one of the supplied values:\n{}\n'.format(kv)
        #         err += 'and the structure of the table:\n{}'.format(self.structure)
        #         raise RuntimeError(err)
        data = ((d[0], pickle.dumps(d[1], protocol=0)) for d in data)
        with self.connection as conn:
            conn.executemany(self.insert_cmd, data)

    def find_one(self, key):
        """
        Searches a SQLite database.

        Inputs:
          key - a single key (as a string)

        Returns: a single unpickled value
        """
        query_str = """SELECT {0}.key, {0}.value
                       FROM {0}
                       WHERE {0}.key LIKE ?""".format(
            self.table_name
        )
        self.cursor.execute(query_str, (key,))
        try:
            return pickle.loads(self.cursor.fetchone()[1])
        except TypeError:
            return None

    def find_many(self, keys):
        """
        Searches a SQLite database.

        Inputs:
          keys - a single key (string) or iterable (list/tuple) containing one or more keys

        Returns: a list of unpickled values
        """
        if type(keys) in STR_TYPES:
            keys = [
                keys,
            ]
        results = []
        for chunk in self.chunker(keys):
            result_chunk = self.cursor.execute(
                """SELECT {0}.key, {0}.value
                   FROM {0}
                   WHERE {0}.key IN ({1})""".format(
                    self.table_name, ",".join("?" * len(chunk))
                ),
                chunk,
            )
            results.extend(result_chunk)
        return [pickle.loads(r[1]) for r in results]

    def find_all(self):
        """
        Returns all values in a SQLite database.

        Inputs:
          keys - a single key (string) or iterable (list/tuple) containing one or more keys

        Returns: a list of unpickled values
        """
        results = self.cursor.execute(
            """SELECT {0}.value
               FROM {0}""".format(
                self.table_name
            )
        )
        return [pickle.loads(r[0]) for r in results]

    def delete(self, keys):
        if type(keys) in STR_TYPES:
            keys = [
                keys,
            ]
        with self.connection as conn:
            conn.executemany(
                """DELETE
                   FROM {0}
                   WHERE {0}.key == ?""".format(
                    self.table_name
                ),
                keys,
            )

    def index(self, fields="key"):
        """
        Indexes the database

        Inputs:
          field - the field or list of field on which to create the index.
        """
        if type(fields) in STR_TYPES:
            fields = [fields]
        index_name = "__".join(fields) + "__index"
        self.cursor.execute(
            "CREATE INDEX {} ON {} ({})".format(
                index_name, self.table_name, ", ".join(fields)
            )
        )

