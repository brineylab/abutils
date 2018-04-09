#!/usr/bin/env python
# filename: utilities.py


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

import operator
import sys

if sys.version_info[0] > 2:
    STR_TYPES = [str, ]
    from functools import reduce
else:
    STR_TYPES = [str, unicode]




#================================
#       USEFUL SNIPPITS
#================================


def nested_dict_lookup(d, key_list):
    '''
    Retrieves a nested value from a dictionary given a list of keys::

        mydict = {'key1a': {'key1b': 'val1b}, 'key2a': {'key2b': {'key2c': 'val2c}}}

        keylist1 = ['key1a', 'key1b']
        nested_dict_lookup(mydict, keylist1)
        >>> val1b

        keylist2 = ['key2a', 'key2b', 'key2c']
        nested_dict_lookup(mydict, keylist2)
        >>> val2c

    Args:

        d (dict): Nested dictionary

        keylist (list(str)): list of strings that correspond to dict keys (in nested order)


    Returns:

        Value from the nested dictionary. If one or more of the nested keys is not present
        in the nested dictionary, `None` will be returned.
    '''
    try:
        return reduce(operator.getitem, key_list, d)
    except KeyError:
        return None


