#!/usr/bin/env python
# filename: alignment.py


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


import multiprocessing as mp
import time
from contextlib import contextmanager
from typing import Iterable, Optional

from tqdm.auto import tqdm

from . import progbar

__all__ = ["monitor_mp_jobs", "monitor_celery_jobs"]


# def monitor_mp_jobs(
#     results: Iterable,
#     start_time: time.time = None,
#     completion_string: str = "\n",
#     print_progress: bool = True,
# ):
#     """
#     Monitors the progress of a set of multiprocessing jobs.

#     Parameters
#     ----------
#     results : Iterable[mp.pool.AsyncResult]
#         A list of AsyncResult objects from multiprocessing jobs.

#     start_time : time.time
#         The time at which the jobs started. If not provided, elapsed time will not be shown.

#     completion_string : str
#         A string to print when the jobs are complete. Default is a newline.

#     print_progress : bool
#         Whether or not to print a progress bar. Default is ``True``.

#     """
#     finished = 0
#     jobs = len(results)
#     while finished < jobs:
#         time.sleep(1)
#         ready = [ar for ar in results if ar.ready()]
#         finished = len(ready)
#         if print_progress:
#             progbar.progress_bar(finished, jobs, start_time=start_time)
#     if print_progress:
#         progbar.progress_bar(
#             finished,
#             jobs,
#             start_time=start_time,
#             complete=True,
#             completion_string=completion_string,
#         )


def monitor_mp_jobs(
    results: Iterable,
    completion_string: Optional[str] = None,
    print_progress: bool = True,
    start_time: Optional[time.time] = None,
):
    """
    Monitors the progress of a set of multiprocessing jobs.

    Parameters
    ----------
    results : Iterable[mp.pool.AsyncResult]
        A list of AsyncResult objects from multiprocessing jobs.

    completion_string : str, optional
        A string to print after the progress bar when the jobs are complete.

    print_progress : bool
        Whether or not to print a progress bar. Default is ``True``.

    start_time : time.time, optional
        DEPRECATED and is ignored. ``tqdm`` automatically handles the start time,
        but the argument is maintained for backwards compatibility.

    """
    jobs = len(results)
    finished = 0
    if print_progress:
        # initialize the TQDM progress bar
        pbar = tqdm(total=jobs)
        prev = 0
    # monitor job progress
    while finished < jobs:
        time.sleep(1)
        ready = [ar for ar in results if ar.ready()]
        finished = len(ready)
        if print_progress:
            update = finished - prev
            pbar.update(update)
            prev = finished
    if print_progress:
        pbar.close()
    if completion_string is not None:
        print(completion_string)


def monitor_celery_jobs(
    results: Iterable,
    start_time: time.time = None,
    completion_string: str = "\n",
    print_progress: bool = True,
):
    """
    Monitors the progress of a set of Celery jobs.

    Parameters
    ----------
    results : Iterable[celery.result.AsyncResult]
        A list of AsyncResult objects from Celery jobs.

    start_time : time.time
        The time at which the jobs started. If not provided, elapsed time will not be shown.

    completion_string : str
        A string to print when the jobs are complete. Default is a newline.

    print_progress : bool
        Whether or not to print a progress bar. Default is ``True``.

    """
    finished = 0
    jobs = len(results)
    while finished < jobs:
        time.sleep(1)
        succeeded = [ar for ar in results if ar.successful()]
        failed = [ar for ar in results if ar.failed()]
        finished = len(succeeded) + len(failed)
        if print_progress:
            progbar.progress_bar(finished, jobs, start_time=start_time)
    if print_progress:
        progbar.progress_bar(
            finished,
            jobs,
            start_time=start_time,
            complete=True,
            completion_string=completion_string,
        )


@contextmanager
def poolcontext(*args, **kwargs):
    pool = mp.Pool(*args, **kwargs)
    yield pool
    pool.terminate()
