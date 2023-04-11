#!/usr/bin/env python
# filename: utils.py


#
# Copyright (c) 2022 Bryan Briney
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


def get_inset_axes_bounds(loc, bbox_to_anchor, width, height):
    if bbox_to_anchor is None:
        loc_dict = {
            "upper left": [0, 1 - height],
            "upper center": [0.5 - width / 2, 1 - height],
            "upper right": [1 - width, 1 - height],
            "center left": [0, 0.5 - height / 2],
            "center": [0.5 - width / 2, 0.5 - height / 2],
            "center right": [1 - width, 0.5 - height / 2],
            "lower left": [0, 0],
            "lower center": [0.5 - width / 2, 0],
            "lower right": [1 - width, 0],
        }
        x0, y0 = loc_dict.get(loc, [0, 0])
    else:
        x, y = bbox_to_anchor[:2]
        loc_dict = {
            "upper left": [x, y - height],
            "upper center": [x - width / 2, y - height],
            "upper right": [x - width, y - height],
            "center left": [x, y - height / 2],
            "center": [x - width / 2, y - height / 2],
            "center right": [x - width, y - height / 2],
            "lower left": [x, y],
            "lower center": [x - width / 2, y],
            "lower right": [x - width, y],
        }
        x0, y0 = loc_dict.get(loc, [0, 0])
    return [x0, y0, width, height]
