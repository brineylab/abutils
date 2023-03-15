#!/usr/bin/env python
# filename: progbar.py


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

import sys
from datetime import datetime


def progress_bar(
    finished,
    total,
    start_time=None,
    extra_info=None,
    autocomplete=True,
    complete=False,
    completion_string="\n",
):
    """
    Prints an ASCII progress bar.

    Each call to ``progress_bar`` will update the progress bar. An example
    of tracking the progress of a list of items would look like::

        job_list = [job1, job2, job3, ... jobN]
        total_jobs = len(job_list)

        #initialize the progress bar
        progress_bar(0, total_jobs)

        # do the jobs
        for i, job in enumerate(job_list):
            do_job(job)
            progress_bar(i + 1, total_jobs)

    Args:

        finished (int): Number of finished jobs.

        total (int): Total number of jobs.

        start_time (datetime): Start time, as a ``datetime.datetime`` object.
            Only required if you want to display execution time alongside
            the progress bar. If not provided, execution time is not shown.
        
        extra_info (str): A string containing extra information to be displayed
            at the end of the progbar string. Examples include the number of failed
            jobs, the name of the job batch currently being processed, etc.

        complete (bool): If `True`, will append `completion_string` to the end
            of the progbar string.

        completion_string (str): Will be appended to the progbar string if
            `complete` is `True`. Default is `'\n'`.

    """
    pct = int(100.0 * finished / total)
    ticks = int(pct / 2)
    spaces = int(50 - ticks)
    if start_time is not None:
        elapsed = (datetime.now() - start_time).seconds
        minutes = int(elapsed / 60)
        seconds = int(elapsed % 60)
        minute_str = "0" * (2 - len(str(minutes))) + str(minutes)
        second_str = "0" * (2 - len(str(seconds))) + str(seconds)
        prog_bar = "\r({}/{}) |{}{}|  {}%  ({}:{})  ".format(
            finished, total, "|" * ticks, " " * spaces, pct, minute_str, second_str
        )
    else:
        prog_bar = "\r({}/{}) |{}{}|  {}%  ".format(
            finished, total, "|" * ticks, " " * spaces, pct
        )
    if extra_info is not None:
        prog_bar += str(extra_info)
    if autocomplete and finished == total:
        prog_bar += completion_string
    if complete:
        prog_bar += completion_string
    sys.stdout.write(prog_bar)
    sys.stdout.flush()
