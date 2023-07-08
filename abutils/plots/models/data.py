#!/usr/bin/env python
# filename: data.py


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


import sys
from typing import Callable, Iterable, Optional, Union

import numpy as np
import pandas as pd

from scipy import stats

from natsort import natsorted


class PlotData:
    """
    Provides methods for processing input data.
    """

    def __init__(
        self,
        data: Union[pd.DataFrame, Iterable, None] = None,
        x: Union[Iterable, str, None] = None,
        y: Union[Iterable, str, None] = None,
        hue: Union[Iterable, str, None] = None,
        column_labels: Optional[Iterable] = None,
        row_labels: Optional[Iterable] = None,
    ):
        """
        Initializes a PlotData object.

        Parameters
        ----------
        data : Union[pd.DataFrame, Iterable, None], optional
            Input data. Can be a ``DataFrame``, an iterable, or ``None``. If ``None``, the x, y, and hue
            data must be provided as iterables to `x`, `y`, and `hue`, respectively.

        x : Union[Iterable, str, None], optional
            The x-axis data. Can be an iterable, a string, or ``None``. If an iterable, the length should
            match the length of all other provided inputs (`data`, `y`, and `hue`). If a string, it should
            be the name of a column in `data`.

        y : Union[Iterable, str, None], optional
            The y-axis data. Can be an iterable, a string, or ``None``. If an iterable, the length should
            match the length of all other provided inputs (`data`, `x`, and `hue`). If a string, it should
            be the name of a column in `data`.

        hue : Union[Iterable, str, None], optional
            The hue data. Can be an iterable, a string, or ``None``. If an iterable, the length should
            match the length of all other provided inputs (`data`, `x`, and `y`). If a string, it should
            be the name of a column in `data`.

        column_labels : Optional[Iterable], optional
            Labels for the columns of `data`. If provided, the length should match the number of columns
            in `data`.

        row_labels : Optional[Iterable], optional
            Labels for the rows of `data`. If provided, the length should match the number of rows in
            `data`.

        """
        self.raw_data = data
        self.raw_x = x
        self.raw_y = y
        self.raw_hue = hue
        self.raw_column_labels = column_labels
        self.raw_row_labels = row_labels

        self.x_name = None
        self.y_name = None
        self.hue_name = None

        self.raw_df = self.process_input(
            data=data,
            x=x,
            y=y,
            hue=hue,
            column_labels=column_labels,
            row_labels=row_labels,
        )
        self.df = self.raw_df.copy()

    @property
    def transform_funcs(self):
        funcs = {
            "log": np.log,
            "log2": np.log2,
            "log10": np.log10,
            "log1p": np.log1p,
            "exp": np.exp,
            "linear": lambda x: x,
        }
        return funcs

    @property
    def agg_funcs(self):
        funcs = {
            "mean": "mean",
            "average": "mean",
            "median": "median",
            "mode": stats.mode,
            "sum": "sum",
            "count": "count",
            "size": "size",
            "nunique": "nunique",
            "max": "max",
            "min": "min",
            "idxmax": "idxmax",
            "idxmin": "idxmin",
            "first": "first",
            "last": "last",
            "cov": "cov",
            "sem": "sem",
            "std": "std",
            "var": "var",
        }
        return funcs

    def process_input(
        self,
        data: Optional[pd.DataFrame] = None,
        x: Union[str, Iterable, None] = None,
        y: Union[str, Iterable, None] = None,
        hue: Union[str, Iterable, None] = None,
        column_labels: Optional[Iterable] = None,
        row_labels: Optional[Iterable] = None,
    ) -> pd.DataFrame:
        """ """
        # input can be provided as iterables passed to x, y and/or hue
        # the column names for x, y and hue are set to 'x', 'y', and 'hue', respectively
        if data is None:
            _data = {}
            _data["x"] = x
            self.x_name = "x"
            _data["y"] = y
            self.y_name = "y"
            if hue is not None:
                _data["hue"] = hue
                self.hue_name = "hue"
            df = pd.DataFrame(_data)

        # or input can be passed to data as a DataFrame or iterable
        else:
            if isinstance(data, pd.DataFrame):
                df = data.copy()
            else:
                df = pd.DataFrame(data)

            # if any of the other inputs are iterables, add them to the DataFrame
            # also, set the x, y, and hue names
            if x is not None:
                if not isinstance(x, str) and len(x) == df.shape[0]:
                    df["x"] = x
                    self.x_name = "x"
                elif isinstance(x, str):
                    self.x_name = x
            if y is not None:
                if not isinstance(y, str) and len(y) == df.shape[0]:
                    df["y"] = y
                    self.y_name = "y"
                elif isinstance(y, str):
                    self.y_name = y
            if hue is not None:
                if not isinstance(hue, str) and len(hue) == df.shape[0]:
                    df["hue"] = hue
                    self.hue_name = "hue"
                elif isinstance(hue, str):
                    self.hue_name = hue

        # set column/row labels
        if column_labels is not None and len(column_labels) == df.shape[1]:
            df.columns = column_labels
        if row_labels is not None and len(row_labels) == df.shape[0]:
            df.index = row_labels
        return df

    def norm(
        self,
        subset: Union[Iterable, str, None] = None,
        axis: Union[int, str] = "columns",
        as_percent: bool = False,
        use_raw: bool = False,
        inplace: bool = True,
    ) -> Optional[pd.DataFrame]:
        """ """
        df = self.raw_df.copy() if use_raw else self.df.copy()
        if isinstance(subset, str):
            subset = [subset]
        multiplier = 100 if as_percent else 1

        # normalize each column separately
        if axis in [0, "columns", "column"]:
            subset = subset if subset is not None else df.columns
            subset = [s for s in subset if s in df.columns]
            for c in subset:
                df[c] = df[c] / df[c].sum() * multiplier

        # normalize each row separately
        if axis in [1, "rows", "row"]:
            df = df.T
            subset = subset if subset is not None else df.columns
            subset = [s for s in subset if s in df.columns]
            for c in subset:
                df[c] = df[c] = df[c] / df[c].sum() * multiplier
            df = df.T

        # update or return
        if inplace:
            self.df = df
        else:
            return df

    def scale(
        self,
        subset: Union[Iterable, str, None] = None,
        axis: Union[int, str, None] = None,
        use_raw: bool = False,
        inplace: bool = True,
    ) -> Optional[pd.DataFrame]:
        """ """
        df = self.raw_df.copy() if use_raw else self.df.copy()
        if isinstance(subset, str):
            subset = [subset]

        # scale across the entire entire dataframe
        if axis is None:
            df = (df - df.min().min()) / (df.max().max() - df.min().min())

        # scale each column separately
        if axis in [0, "columns", "column"]:
            subset = subset if subset is not None else df.columns
            subset = [s for s in subset if s in df.columns]
            for c in subset:
                df[c] = (df[c] - df[c].min()) / (df[c].max() - df[c].min())

        # scale each row separately
        if axis in [1, "rows", "row"]:
            df = df.T
            subset = subset if subset is not None else df.columns
            subset = [s for s in subset if s in df.columns]
            for c in subset:
                df[c] = (df[c] - df[c].min()) / (df[c].max() - df[c].min())
            df = df.T

        # update or return
        if inplace:
            self.df = df
        else:
            return df

    def transform(
        self,
        subset: Union[Iterable, str, None] = None,
        func: Union[Callable, str, None] = None,
        use_raw: bool = False,
        inplace: bool = True,
    ) -> Optional[pd.DataFrame]:
        """ """
        if isinstance(func, str):
            func = self.transform_funcs.get(func.lower(), None)
            if func is None:
                func_block = "\n  - ".join(self.transform_funcs.keys())
                err = f"\nERROR: invalid transform function.\n"
                err += f"Options are: {func_block}\n"
                err + +f"You provided: '{func}'\n"
                print(err)
                sys.exit()
        elif not callable(func):
            return
        df = self.raw_df.copy() if use_raw else self.df.copy()

        # subset of the dataframe to transform
        if isinstance(subset, str):
            subset = [subset]
        if subset is None:
            subset = df.columns
        subset = [s for s in subset if s in df.columns]

        # transform
        for c in subset:
            df[c] = func(df[c])

        # update or return
        if inplace:
            self.df = df
        else:
            return df

    def reorder(
        self,
        by: Union[Iterable, str, None] = None,
        axis: Union[int, str] = "rows",
        row_order: Optional[Iterable] = None,
        column_order: Optional[Iterable] = None,
        ascending: bool = True,
        use_raw: bool = False,
        inplace: bool = True,
    ) -> Optional[pd.DataFrame]:
        """ """
        df = self.raw_df.copy() if use_raw else self.df.copy()

        # reorder based on row/column values
        if by is not None:
            if axis in [0, "columns", "column"]:
                df = df.sort_values(by=by, axis=0, ascending=ascending)
            if axis in [1, "rows", "row"]:
                df = df.sort_values(by=by, axis=1, ascending=ascending)

        # reorder based on column names
        if column_order is not None:
            df = df[column_order]

        # reorder based on row names (index)
        if row_order is not None:
            df = df.loc[row_order]

        # update or return
        if inplace:
            self.df = df
        else:
            return df

    def to_squareform(
        self,
        rows: Optional[str] = None,
        columns: Optional[str] = None,
        values: Optional[str] = None,
        agg: Union[Iterable, Callable, str] = "size",
        use_raw: bool = False,
        reset_index: bool = False,
        inplace: bool = True,
    ) -> Optional[pd.DataFrame]:
        """ """
        df = self.raw_df.copy() if use_raw else self.df.copy()
        # get aggregation function(s)
        if isinstance(agg, str) or callable(agg):
            agg = [agg]
        agg_funcs = []
        for a in agg:
            if isinstance(a, str):
                a = self.agg_funcs.get(a.lower(), None)
                if a is None:
                    agg_block = "\n  - ".join(self.agg_funcs.keys())
                    err = f"\nERROR: invalid aggregation function.\n"
                    err += f"Options are: {agg_block}\n"
                    err + +f"You provided: '{a}'\n"
                    print(err)
                    sys.exit()
            elif not callable(a):
                agg_block = "\n  - ".join(self.agg_funcs.keys())
                err = f"\nERROR: invalid aggregation function.\n"
                err += f"<agg> must be a callable function or the name of a built-in aggregation function.\n"
                err += f"Built-in options are: {agg_block}\n"
                err + +f"You provided: '{a}'\n"
                print(err)
                sys.exit()
            agg_funcs.append(a)
        agg_funcs = [a for a in agg_funcs if a is not None]
        if len(agg_funcs) == 1:
            agg_funcs = agg_funcs[0]

        # group by columns/rows
        if columns is None:
            columns = self.x_name
        if rows is None:
            rows = self.y_name
        g = df.groupby([rows, columns])

        # aggregate
        agg_df = g.agg(agg_funcs)
        if values is not None:
            agg_df = agg_df[values]

        # pivot to squareform
        sq_df = agg_df.unstack()
        if reset_index:
            sq_df.reset_index(inplace=True)

        # update or return
        if inplace:
            self.df = sq_df
        else:
            return sq_df
