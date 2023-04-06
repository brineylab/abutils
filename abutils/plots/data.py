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


from copy import deepcopy
from typing import Iterable, Optional, Union

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
