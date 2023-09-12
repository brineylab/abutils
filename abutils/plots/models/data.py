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

from ...core.sequence import Sequence


class PlotData:
    """
    Provides methods for processing input data.
    """

    def __init__(
        self,
        x: Union[Iterable, str, None] = None,
        y: Union[Iterable, str, None] = None,
        data: Union[pd.DataFrame, Iterable, dict, None] = None,
        sequences: Optional[Iterable[Sequence]] = None,
        hue: Union[Iterable, str, None] = None,
        column_labels: Optional[Iterable] = None,
        row_labels: Optional[Iterable] = None,
    ):
        """
        Initializes a PlotData object.

        Parameters
        ----------
        data : Union[pd.DataFrame, Iterable, None], optional
            Input data. Can be a ``DataFrame``, a ``dict``, or an iterable. If both `data` and `sequences`
            are provided, only `sequences` is used. If neither `data` nor `sequences` are provided, the
            x, y, and hue data must be provided as iterables to `x`, `y`, and `hue`, respectively.

        sequences : Optional[Iterable[Sequence]], optional
            Input data. Must be an iterable of ``abutils.Sequence`` objects. If both `data` and `sequences`
            are provided, only `sequences` is used. If neither `data` nor `sequences` are provided, the
            x, y, and hue data must be provided as iterables to `x`, `y`, and `hue`, respectively.

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
        self.raw_sequences = sequences
        self.raw_x = x
        self.raw_y = y
        self.raw_hue = hue
        self.raw_column_labels = column_labels
        self.raw_row_labels = row_labels

        self.x_name = None
        self.y_name = None
        self.hue_name = None

        self.raw_df = self.process_input(
            x=x,
            y=y,
            data=data,
            sequences=sequences,
            hue=hue,
            column_labels=column_labels,
            row_labels=row_labels,
        )
        self.df = self.raw_df.copy()

    @property
    def transform_funcs(self) -> dict:
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
    def agg_funcs(self) -> dict:
        funcs = {
            "mean": "mean",
            "average": "mean",
            "median": "median",
            "mode": pd.Series.mode,
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
        x: Union[str, Iterable, None] = None,
        y: Union[str, Iterable, None] = None,
        data: Union[pd.DataFrame, Iterable, dict, None] = None,
        sequences: Optional[Iterable[Sequence]] = None,
        hue: Union[str, Iterable, None] = None,
        column_labels: Optional[Iterable] = None,
        row_labels: Optional[Iterable] = None,
    ) -> pd.DataFrame:
        """
        Reads input datasets and returns a ``DataFrame`` containing x, y, and hue data.
        """
        # input can be provided as iterables passed to x, y and/or hue
        # the column names for x, y and hue are set to 'x', 'y', and 'hue', respectively
        if all([data is None, sequences is None]):
            _data = {}
            _data["x"] = x
            self.x_name = "x"
            _data["y"] = y
            self.y_name = "y"
            if hue is not None:
                _data["hue"] = hue
                self.hue_name = "hue"
            df = pd.DataFrame(_data)

        # input can also be passed to data as a DataFrame or iterable
        # or to sequences as an iterable of abutils.Sequence objects
        else:
            if sequences is not None:
                df = [s.annotations for s in sequences]
            elif isinstance(data, pd.DataFrame):
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

        # drop any columns that aren't x, y, or hue
        df = df[[self.x_name, self.y_name, self.hue_name]]

        # set column/row labels
        if column_labels is not None and len(column_labels) == df.shape[1]:
            df.columns = column_labels
        if row_labels is not None and len(row_labels) == df.shape[0]:
            df.index = row_labels
        return df

    def dropna(
        self,
        axis: Union[int, str] = "rows",
        how: str = "any",
        thresh: Optional[int] = None,
        subset: Optional[Union[str, Iterable]] = None,
        use_raw: bool = False,
        inplace: bool = True,
    ) -> Optional[pd.DataFrame]:
        """
        Drop missing or NA values from the data. Default settings will drop all rows with any
        missing values. See the `pandas.DataFrame.dropna`_ documentation for more information.

        .. _pandas.DataFrame.dropna:
        https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.dropna.html

        """
        df = self.raw_df.copy() if use_raw else self.df.copy()
        df.dropna(axis=axis, how=how, thresh=thresh, subset=subset, inplace=True)
        if inplace:
            self.df = df
        else:
            return df

    def fillna(
        self,
        value: Optional[Union[int, float, str]] = None,
        axis: Union[int, str] = "columns",
        subset: Optional[Union[str, Iterable]] = None,
        method: Optional[str] = None,
        limit: Optional[int] = None,
        downcast: Optional[dict] = None,
        use_raw: bool = False,
        inplace: bool = True,
    ) -> Optional[pd.DataFrame]:
        """
        Fill missing or NA values in the data. The default settings will fill all missing values
        with ``0``. See the `pandas.DataFrame.fillna`_ documentation for more information.

        Parameters
        ----------
        value : Optional[Union[int, float, str]]
            Value to use to fill holes. If ``None``, the ``method`` parameter must be specified.

        axis : Union[int, str]
            Only used in combination with `subset`. If ``0``, ``"column"``, or ``"columns"``,
            columns in `subset` will be filled. If ``1``, ``"row"``, or ``"rows"``, rows in
            `subset` will be filled. Default is ``"columns"``, so columns in `subset` will be
            filled.

        subset : Optional[Union[str, Iterable]]
            Subset of columns or rows to fill. If ``None``, all columns or rows will be filled.

        method : Optional[str]
            Method to use to fill holes. See the `pandas.DataFrame.fillna`_ documentation for
            available methods.

        limit : Optional[int]
            Maximum number of consecutive NaN values to forward/backward fill.

        downcast : Optional[dict]
            A dict of item->dtype of what to downcast if possible, or the string "infer" which
            will try to downcast to an appropriate equal type (e.g. float64 to int64 if possible).
            See the `pandas.DataFrame.fillna`_ documentation for more information.

        use_raw : bool
            If ``True``, fill missing values in the raw data. If ``False``, fill missing values in
            the processed data. Default is ``False``

        inplace : bool
            If ``True``, fill missing values in place. If ``False``, return a copy of the DataFrame
            with missing values filled. Default is ``True``.

        Returns
        -------
        Optional[pd.DataFrame]
            If ``inplace`` is ``False``, returns a copy of the DataFrame with missing values filled.
            Otherwise, returns ``None``.


        .. _pandas.DataFrame.fillna:
        https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.fillna.html

        """
        df = self.raw_df.copy() if use_raw else self.df.copy()
        if method is None:
            value = 0 if value is None else value

        # fill columns
        if axis in [0, "column", "columns"]:
            subset = subset if subset is not None else df.columns
            subset = [s for s in subset if s in df.columns]
            for s in subset:
                df[s].fillna(
                    value=value,
                    method=method,
                    limit=limit,
                    downcast=downcast,
                    inplace=True,
                )

        # fill rows
        if axis in [1, "row", "rows"]:
            df = df.T
            subset = subset if subset is not None else df.columns
            subset = [s for s in subset if s in df.columns]
            for s in subset:
                df[s].fillna(
                    value=value,
                    method=method,
                    limit=limit,
                    downcast=downcast,
                    inplace=True,
                )
            df = df.T

        # update or return
        if inplace:
            self.df = df
        else:
            return df

    def norm(
        self,
        subset: Union[Iterable, str, None] = None,
        axis: Union[int, str] = "columns",
        as_percent: bool = False,
        use_raw: bool = False,
        inplace: bool = True,
    ) -> Optional[pd.DataFrame]:
        """
        Normalize the data. By default, each column is normalized to sum to 1.

        Parameters
        ----------
        subset : Union[Iterable, str, None]
            The columns or rows to normalize. If ``None``, all columns or rows will be normalized.
            Default is ``None``.

        axis : Union[int, str]
            The axis to normalize. If ``0``, ``"column"``, or ``"columns"``, each column will
            be  normalized. If ``1``, ``"row"``, or ``"rows"``, each row will be normalized.
            Default is ``"columns"``.

        as_percent : bool
            If ``True``, the data will be multiplied by ``100`` after normalization. Default is
            ``False``.

        use_raw : bool
            If ``True``, the raw data will be normalized. Otherwise, the processed data will be
            normalized. Default is ``False``.

        inplace : bool
            If ``True``, the normalized data will replace ``self.df``. Otherwise, the normalized
            ``DataFrame`` will be returned. Default is ``True``.

        Returns
        -------
        Optional[pd.DataFrame]
            If ``inplace`` is ``False``, the normalized ``DataFrame`` will be returned.

        """
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
        """
        Scale the data to be between 0 and 1. By default, the entire dataset is scaled
        according to the global minimum and maximum values in the dataset. Alternatively, each
        column or row can be scaled separately.

        Parameters
        ----------
        subset : Union[Iterable, str, None]
            The columns or rows to scale. Can be either a single row/column name, or an iterable
            of row/column names. If ``None``, all columns or rows will be scaled. Default is
            ``None``.

        axis : Union[int, str, None]
            The axis to scale. If ``0``, ``"column"``, or ``"columns"``, each column will be
            scaled separately. If ``1``, ``"row"``, or ``"rows"``, each row will be scaled separately.
            If ``None``, the entire dataset will be scaled. Default is ``None``.

        use_raw : bool
            If ``True``, the raw data will be scaled. Otherwise, the processed data will be scaled.
            Default is ``False``.

        inplace : bool
            If ``True``, the scaled data will replace ``self.df``. Otherwise, the scaled
            ``DataFrame`` will be returned. Default is ``True``.

        Returns
        -------
        Optional[pd.DataFrame]
            If ``inplace`` is ``False``, the scaled ``DataFrame`` will be returned.

        """
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
        func: Union[Callable, str, None] = None,
        subset: Union[Iterable, str, None] = None,
        use_raw: bool = False,
        inplace: bool = True,
    ) -> Optional[pd.DataFrame]:
        """
        Transform the data (log, expoonential, etc). By default, the function is applied
        to all columns the dataset. Alternatively, a subset of columns can be transformed.

        Parameters
        ----------
        func : Union[Callable, str, None]
            The function to apply to the data. If ``None``, the data will not be transformed.
            If a string, the function must be one of the following:
                * ``"log"``: natural log
                * ``"log10"``: log base 10
                * ``"log2"``: log base 2
                * ``"log1p"``: natural log of (value + 1)
                * ``"exp"``: exponential
                * ``"linear"``: no transformation

        subset : Union[Iterable, str, None]
            The columns to transform. Can be either a single column name, or an iterable
            of column names. If ``None``, all columns will be transformed. Default is
            ``None``.

        use_raw : bool
            If ``True``, the raw data will be transformed. Otherwise, the processed data will be
            transformed. Default is ``False``.

        inplace : bool
            If ``True``, the transformed data will replace ``self.df``. Otherwise, the transformed
            ``DataFrame`` will be returned. Default is ``True``.

        Returns
        -------
        Optional[pd.DataFrame]
            If ``inplace`` is ``False``, the transformed ``DataFrame`` will be returned.

        """
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

    def clip(
        self,
        lower: Union[int, float, None] = None,
        upper: Union[int, float, None] = None,
        axis: Union[int, str] = "columns",
        subset: Union[Iterable, str, None] = None,
        use_raw: bool = False,
        inplace: bool = True,
    ) -> Optional[pd.DataFrame]:
        """
        Clip the data. Values less than `lower` will be replaced with `lower` and values
        greater than `upper` will be replaced with `upper`. By default, the entire dataset
        will be clipped. Alternatively, a subset of columns can be clipped.

        Parameters
        ----------
        lower : Union[int, float, None]
            The lower bound. If ``None``, no lower bound will be applied. Default is ``None``.

        upper : Union[int, float, None]
            The upper bound. If ``None``, no upper bound will be applied. Default is ``None``.

        axis : Union[int, str, None]
            The axis to clip. Only used if `subset` is also provided. If ``0``, ``"column"``,
            or ``"columns"``, only column names in `subset` will be clipped. If ``1``,
            ``"row"``, or ``"rows"``, only row names in `subset` will be clipped. Default is
            ``"columns"``.

        subset : Union[Iterable, str, None]
            The columns or rows to clip. Can be either a single column/row name, or an iterable
            of column/row names. If ``None``, all columns/rows will be clipped. Default is
            ``None``.

        use_raw : bool
            If ``True``, the raw data will be clipped. Otherwise, the processed data will be
            clipped. Default is ``False``.

        inplace : bool
            If ``True``, the transformed data will replace ``self.df``. Otherwise, the transformed
            ``DataFrame`` will be returned. Default is ``True``.

        Returns
        -------
        Optional[pd.DataFrame]
            If ``inplace`` is ``False``, the transformed ``DataFrame`` will be returned.

        """
        df = self.raw_df.copy() if use_raw else self.df.copy()
        if isinstance(subset, str):
            subset = [subset]

        # clip by column
        if axis in [0, "columns", "column"]:
            subset = subset if subset is not None else df.columns
            subset = [s for s in subset if s in df.columns]
            for c in subset:
                df[c] = df[c].clip(lower=lower, upper=upper)

        # clip by row
        if axis in [1, "rows", "row"]:
            df = df.T
            subset = subset if subset is not None else df.columns
            subset = [s for s in subset if s in df.columns]
            for c in subset:
                df[c] = df[c].clip(lower=lower, upper=upper)
            df = df.T

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

    def groupby(
        self,
        by: Union[Iterable, str],
        axis: Union[int, str] = "rows",
        values: Optional[Union[Iterable, str]] = None,
        agg: Union[Iterable, Callable, str, dict] = "size",
        use_raw: bool = False,
        inplace: bool = True,
    ) -> Optional[pd.DataFrame]:
        """ 
        Group the data by row or column values.

        Parameters
        ----------
        by : Union[Iterable, str]
            The column or row values to group by. Can be either a single column/row name, or an
            iterable of column/row names.

        axis : Union[int, str]
            The axis to group by. If ``0``, ``"column"``, or ``"columns"``, the data will be
            grouped by column values. If ``1``, ``"row"``, or ``"rows"``, the data will be
            grouped by row values. Default is ``"rows"``.

        values : Optional[Union[Iterable, str]]
            The columns or rows to aggregate. Can be either a single column/row name, or an
            iterable of column/row names. If ``None``, all columns/rows will be aggregated.
            Default is ``None``.

        agg : Union[Iterable, Callable, str, dict]
            The aggregation function(s) to apply to the grouped data. Can be a single function,
            an iterable of functions, or a dictionary mapping column/row names to functions.
            Default is ``"size"``.

        use_raw : bool
            If ``True``, the raw data will be grouped. Otherwise, the processed data will be
            grouped. Default is ``False``.

        inplace : bool
            If ``True``, the transformed data will replace ``self.df``. Otherwise, the transformed
            ``DataFrame`` will be returned. Default is ``True``.

        Returns
        -------
        Optional[pd.DataFrame]
            If ``inplace`` is ``False``, the transformed ``DataFrame`` will be returned.

        """
        df = self.raw_df.copy() if use_raw else self.df.copy()
        # do the groupby
        if axis in [1, "rows", "row"]:
            g = df.groupby(by=by, axis=1)
        else:
            g = df.groupby(by=by, axis=0)
        # filter values
        if values is not None:
            g = g[values]
        # aggregate
        agg_df = g.agg(agg)
        if inplace:
            self.df = agg_df
        else:
            return agg_df

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
        """ 
        Pivots the data into a squareform matrix.

        Parameters
        ----------
        rows : Optional[str]
            The column name to use for the rows. If ``None``, ``self.x_name`` will be used. 
            Default is ``None``.
        
        columns : Optional[str]
            The column name to use for the columns. If ``None``, ``self.y_name`` will be used. 
            Default is ``None``.

        values : Optional[str]
            The column name to use for the aggregated values. If ``None``, all columns will 
            be used. Default is ``None``.

        agg : Union[Iterable, Callable, str]
            The aggregation function(s) to use. If a single function is provided, it will be
            applied to all values. If multiple functions are provided, they will be applied
            to the corresponding values. If a string is provided, it must be a valid aggregation
            function. Default is ``"size"``.

        use_raw : bool
            If ``True``, the raw data will be used. Otherwise, the processed data will be used.
            Default is ``False``.

        reset_index : bool
            If ``True``, the index will be reset. Default is ``False``.

        inplace : bool
            If ``True``, the transformed data will replace ``self.df``. Otherwise, the transformed
            ``DataFrame`` will be returned. Default is ``True``.

        Returns
        -------
        Optional[pd.DataFrame]
            If ``inplace`` is ``False``, the transformed ``DataFrame`` will be returned.

        """
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
                err += f"<agg> must be a function or the name of a built-in aggregation function.\n"
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
        if values is not None:
            agg_df = g[values].agg(agg_funcs)
        else:
            agg_df = g.agg(agg_funcs)

        # pivot to squareform
        sq_df = agg_df.unstack()
        if reset_index:
            sq_df.reset_index(inplace=True)

        # update or return
        if inplace:
            self.df = sq_df
        else:
            return sq_df
