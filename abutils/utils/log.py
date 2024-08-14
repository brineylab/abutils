#!/usr/bin/env python
# filename: log.py


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


import logging
import os

__all__ = [
    "SingleLineHandler",
    "NotebookLogger",
    "setup_logging",
    "null_logger",
    "get_logger",
]


class SingleLineHandler(logging.StreamHandler):
    """
    A handler that writes to the stream without a newline at the end. Can be used
    interchangeably with the built-in ``logging.StreamHandler``.
    """

    def emit(self, record):
        try:
            msg = self.format(record)
            stream = self.stream
            # Prevent newline character at the end
            stream.write(msg)
            self.flush()
        except Exception:
            self.handleError(record)


class NotebookLogger:
    """
    A logger that prints instead of logging, which works better in Jupyter notebooks.
    """

    def __init__(self, verbose: bool = True, end: str = "\n"):
        self.verbose = verbose
        self.end = end

    def info(self, msg, end: str = "", *args, **kwargs):
        if self.verbose:
            end = end if end else self.end
            print(msg, end=end, *args, **kwargs)

    def debug(self, msg, end: str = "", *args, **kwargs):
        if self.verbose:
            end = end if end else self.end
            print(msg, end=end, *args, **kwargs)


def setup_logging(
    logfile: str,
    print_log_location: bool = True,
    add_stream_handler: bool = True,
    single_line_handler: bool = False,
    debug: bool = False,
):
    """
    Set up logging using the built-in ``logging`` package.

    An optional stream handler can be added to the ``logger``, so that logs at or above
    ``logging.INFO`` level are printed to screen as well as written to the log file.

    Parameters
    ----------

    logfile : str
        Path to the log file. If the parent directory
        does not exist, it will be created. Required.

    print_log_location : bool, default=True
        If ``True``, the log path will be written to the log upon initialization.

    add_stream_handler : bool, default=True
        If ``True``, a stream handler will be added to the log.

    debug : bool, default=False
        If true, the log level will be set to ``logging.DEBUG``.
        If ``False``, the log level will be set to ``logging.INFO``.

    """
    from ..io import make_dir

    log_dir = os.path.dirname(logfile)
    make_dir(log_dir)
    fmt = "[%(levelname)s] %(name)s %(asctime)s %(message)s"
    if debug:
        # log everything at DEBUG level and above
        logging.basicConfig(
            filename=logfile, filemode="w", format=fmt, level=logging.DEBUG
        )
    else:
        # log everything at INFO level and above
        logging.basicConfig(
            filename=logfile, filemode="w", format=fmt, level=logging.INFO
        )
    logger = logging.getLogger("log")
    if add_stream_handler:
        logger = _add_stream_handler(
            logger=logger,
            handler_class=SingleLineHandler
            if single_line_handler
            else logging.StreamHandler,
        )
    if print_log_location:
        logger.info("LOG LOCATION: {}".format(logfile))


def null_logger():
    """
    Get a logger that does nothing.
    """
    return logging.getLogger("null")


def get_logger(
    name: str = None, add_stream_handler: bool = True, single_line_handler: bool = False
) -> logging.Logger:
    """
    Get a logging handle.

    As with ``setup_logging``, an optional stream handler can be added to the logger.

    Parameters
    ----------

    name : str
        Name of the log handle.

    add_stream_handler : bool, default=True
        If ``True``, a stream handler will be added to the log.

    Returns
    -------

    logger : logging.Logger

    """
    logger = logging.getLogger(name)
    # only add the stream handler if there are no other handlers
    if len(logger.handlers) == 0 and add_stream_handler:
        logger = _add_stream_handler(
            logger=logger,
            handler_class=SingleLineHandler
            if single_line_handler
            else logging.StreamHandler,
        )
    return logger


def _add_stream_handler(
    logger: logging.Logger, handler_class: type = logging.StreamHandler
) -> logging.Logger:
    """
    Add a stream handler to the logger.

    Parameters
    ----------

    logger : logging.Logger
        The logger to add the stream handler to.

    handler_class : type, default=logging.StreamHandler
        The class of the handler to add to the logger.

    Returns
    -------

    logger : logging.Logger
        The logger with the stream handler added.
    """
    formatter = logging.Formatter("%(message)s")
    handler = handler_class()
    handler.setFormatter(formatter)
    handler.setLevel(logging.INFO)
    logger.addHandler(handler)
    return logger
