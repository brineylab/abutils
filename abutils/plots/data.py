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


from collections import Counter
from copy import deepcopy
from typing import Callable, Iterable, Optional, Union

import pandas as pd

from ..core.sequence import Sequence


def process_input_data(
    x: Union[str, Iterable, None] = None,
    y: Union[str, Iterable, None] = None,
    hue: Union[str, Iterable, None] = None,
    sequences: Optional[Iterable[Sequence]] = None,
    data: Union[pd.DataFrame, dict, None] = None,
):
    """
    Processes input data for plotting functions.


    """
    # sequences are provided
    _data = {}
    if sequences is not None:
        # X values
        if isinstance(x, str):
            _data[x] = [s[x] for s in sequences]
        elif isinstance(x, Iterable):
            if len(x) == len(sequences):
                _data["x"] = x
                x = "x"
            else:
                raise ValueError(
                    "If `x` is an iterable, it must be the same length as `sequences`."
                )
        # Y values
        if isinstance(y, str):
            _data[y] = [s[y] for s in sequences]
        elif isinstance(y, Iterable):
            if len(y) == len(sequences):
                _data["y"] = y
                y = "y"
            else:
                raise ValueError(
                    "If `y` is an iterable, it must be the same length as `sequences`."
                )
        # hue values
        if isinstance(hue, str):
            _data[hue] = [s[hue] for s in sequences]
        elif isinstance(hue, Iterable):
            if len(hue) == len(sequences):
                _data["hue"] = hue
                hue = "hue"
            else:
                raise ValueError(
                    "If `hue` is an iterable, it must be the same length as `sequences`."
                )
        data = pd.DataFrame(_data, index=[s.id for s in sequences])
    # data is provided
    elif data is not None:
        # X values
        if isinstance(x, str):
            _data[x] = data[x]
        elif isinstance(x, Iterable):
            _data["x"] = x
            x = "x"
        # Y values
        if isinstance(y, str):
            _data[y] = data[y]
        elif isinstance(y, Iterable):
            _data["y"] = y
            y = "y"
        # hue values
        if isinstance(hue, str):
            _data[hue] = data[hue]
        elif isinstance(hue, Iterable):
            _data["hue"] = hue
            hue = "hue"
        data = pd.DataFrame(_data)
    # neither sequences nor data are provided
    else:
        if x is not None:
            _data["x"] = x
            x = "x"
        if y is not None:
            _data["y"] = y
            y = "y"
        if hue is not None:
            _data["hue"] = hue
            hue = "hue"
        data = pd.DataFrame(_data)

    return data, x, y, hue


class InputData:
    """ """

    def __init__(
        self,
        data: Union[pd.DataFrame, dict, None] = None,
        sequences: Optional[Iterable[Sequence]] = None,
        x: Union[str, Iterable, None] = None,
        y: Union[str, Iterable, None] = None,
        values: Union[str, Iterable, dict, None] = None,
        categories: Union[str, Iterable, None] = None,
        hue: Union[str, Iterable, None] = None,
        agg_method: Union[str, Callable] = "count",
    ):
        """
        Data can be either a dict or a DataFrame

        Sequences must be an iterable of Sequence objects.

        X and Y values can be provided as a string (column name), an iterable, or
        ``None``. If ``None``, the data will be aggregated by the ``agg_method``.

        """
        self.input_data = data
        self.data = deepcopy(data) if data is not None else {}
        self.sequences = sequences
        self.input_x = x
        self.input_y = y
        self.input_values = values
        self.input_categories = categories
        self.input_hue = hue
        self.x_name = None
        self.y_name = None
        self.values_name = None
        self.categories_name = None
        self.hue_name = None
        self.agg_method = agg_method

    @property
    def agg(self):
        if callable(self.agg_method):
            return self.agg_method
        elif self.agg_method == "count":
            return self._count

    # @property
    # def needs_aggragation(self):
    #     """
    #     Determines whether aggregation is needed.

    #     Note::
    #     this should only be called after `load`, since certain cases (such as
    #     when a ``dict` is provided as `values`) will cause this function to
    #     incorrectly return ``True`` when aggregation is not needed.
    #     """
    #     if isinstance(self.x, (list, tuple)) and self.y is None:
    #         return True
    #     if isinstance(self.values, (list, tuple)) and self.categories is None:
    #         return True
    #     return False

    def load(self):
        # reformat inputs if `x` or `categories` is a ``dict``
        if isinstance(self.input_x, dict):
            self.x_vals = list(self.input_x.keys())
            self.y_vals = list(self.input_x.values())
        else:
            self.x_vals = self.input_x
            self.y_vals = self.input_y
        if isinstance(self.input_categories, dict):
            self.categories_vals = list(self.input_categories.keys())
            self.values_vals = list(self.input_categories.values())
        else:
            self.categories_vals = self.input_categories
            self.values_vals = self.input_values
        # read all of the input values (from either sequences or data)
        d = {}
        val_names = ["x", "y", "categories", "values", "hue"]
        vals = [
            self.x_vals,
            self.y_vals,
            self.categories_vals,
            self.values_vals,
            self.hue_vals,
        ]
        for val, name in zip(vals, val_names):
            if val is not None:
                val_name, parsed_vals = self.parse_vals(val, name=name)
                d[val_name] = parsed_vals
                setattr(self, f"{name}_name", val_name)
        df = pd.DataFrame(d)
        # split into hue groups
        if self.input_hue is not None:
            hue_groups = []
            for h in df[self.hue_name].unique():
                hue_groups.append(df[df[self.hue_name] == h])
        else:
            hue_groups = [df]
        # process hue groups separately
        hue_dfs = []
        for group in hue_groups:
            hue_data = {}
            if any([self.x_name is not None, self.y_name is not None]):
                hue_x, hue_y = self.parse_xy(data=group, x=self.x, y=self.y)
                hue_data[self.x_name] = hue_x
                hue_data[self.y_name] = hue_y
                if self.input_hue is not None:
                    hue_val = group[self.hue_name].unique()[0]
                    hue_data[self.hue_name] = [hue_val] * len(hue_x)
            else:
                hue_categories, hue_values = self.parse_catval(
                    data=group,
                    categories=self.input_categories,
                    values=self.input_values,
                )
                hue_data[self.categories_name] = hue_categories
                hue_data[self.values_name] = hue_values
                if self.input_hue is not None:
                    hue_val = group[self.hue_name].unique()[0]
                    hue_data[self.hue_name] = [hue_val] * len(hue_x)
            hue_dfs.append(pd.DataFrame(hue_data))
        self.data = pd.concat(hue_dfs)

    def parse_xy(self, data, x, y):
        if y is None:
            if isinstance(x, (list, tuple)):
                x, y = self.agg(x)
            elif isinstance(x, str):
                if x in data:
                    x, y = self.agg(data[x])
                elif self.sequences is not None:
                    vals = [s[x] for s in self.sequences]
                    x, y = self.agg(vals)
                else:
                    raise ValueError(
                        "If `x` is a ``str``, it must be a column name in `data` or a "
                        "property of `sequences`."
                    )
            else:
                raise ValueError("`x` must be a ``str`` or ``Iterable``.")
        else:
            x = self.parse_data_item(x)
            y = self.parse_data_item(y)
        return x, y

    def parse_catval(self, data, categories, values):
        if values is None:
            if isinstance(categories, dict):
                values = list(categories.values())
                categories = list(categories.keys())
            elif isinstance(categories, (list, tuple)):
                categories, values = self.agg(categories)
            elif isinstance(categories, str):
                if categories in data:
                    categories, values = self.agg(data[categories])
                elif self.sequences is not None:
                    vals = [s[categories] for s in self.sequences]
                    categories, values = self.agg(vals)
                else:
                    raise ValueError(
                        "If `categories` is a ``str``, it must be a column name in `data` or a "
                        "property of `sequences`."
                    )
            else:
                raise ValueError(
                    "`categories` must be a ``dict``, ``str``, or ``Iterable``."
                )
        else:
            categories = self.parse_data_item(categories)
            values = self.parse_data_item(values)
        return categories, values

    def parse_vals(self, item, name=None):
        if isinstance(item, str):
            if item in self.data:
                name = item
                vals = self.data[item]
            elif self.sequences is not None:
                name = item
                return [s[item] for s in self.sequences]
            else:
                raise ValueError(
                    f"If `{name}` is a ``str``, it must be a column name in `data` or a "
                    "property of `sequences`."
                )
        return name, item

    def _count(self, val_list):
        counts = Counter(val_list).most_common()
        return list(counts.keys()), list(counts.values())
