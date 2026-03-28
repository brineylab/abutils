# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT


import logging
import os
from typing import Optional

from ..io import make_dir

_all__ = [
    "LoggingMixin",
    "SingleLineHandler",
    "NotebookLogger",
    "SimpleLogger",
    "setup_logging",
    "null_logger",
    "get_logger",
]


class LoggingMixin:
    """
    Mixin for logging messages and exceptions.

    Methods:
    --------
    log(self, *args, separator: str = " ") -> None:
        Records a log message, with each argument joined by the separator.

    exception(self, *args, separator: str = "\n") -> None:
        Records an exception, with each argument joined by the separator.

    format_log(self, separator: str = "\n") -> str:
        Formats the log, including exceptions, as a string to be written to file or stdout.

    format_exceptions(self, separator: str = "\n") -> str:
        Formats the exceptions as a string to be written to file or stdout.

    """

    def __init__(self):
        super().__init__()
        self.logs = []
        self.exceptions = []

    def log(self, *args, separator: str = " ") -> None:
        """
        Records a log message

        Parameters:
        ----------
        *args : str
            The log message to record.

        separator : str, default: " "
            The separator to join the arguments.

        """
        log_str = separator.join([str(a) for a in args])
        self.logs.append(log_str)

    def exception(self, *args, separator: str = "\n") -> None:
        """
        Records an exception

        Parameters:
        ----------
        *args : str
            The exception to record.

        separator : str, default: "\n"
            The separator to join the arguments.

        """
        exception_str = separator.join([str(a) for a in args])
        self.exceptions.append(exception_str)

    def format_log(self, separator="\n") -> str:
        """
        Formats the log, including exceptions.

        Log formatting will only be performed on sequences that had an
        error during annotation, unless AbStar is run in debug
        mode. In debug mode, all sequences will be logged.

        Parameters:
        ----------
        separator : str, default: "\n"
            The separator to join the log messages.

        Returns:
        --------
        str: Formatted log string.

        """
        output = ""
        output += separator.join(self.logs)
        if self.exceptions:
            output += "\n\n"
            output += self._format_exceptions()
        output += "\n\n"
        return output

    def _format_exceptions(self) -> str:
        """
        Formats the exceptions as a string to be written to file or stdout.

        Returns:
        --------
        str: Formatted exceptions string.

        """
        exception_str = "\n\nEXCEPTIONS\n"
        exception_str += "----------\n\n"
        exception_str += "\n\n".join([e for e in self.exceptions])
        return exception_str


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


class SimpleLogger:
    """
    A simple logger that writes log and exception info to a file.

    The log file can be specified at initialization, or a log file can be provided to
    the ``write()`` method.
    """

    def __init__(self, log_file: Optional[str] = None):
        self.logs = []
        self.exceptions = []
        self.log_file = log_file
        self.formatted = ""

    def log(self, *args, separator: str = " ") -> None:
        """
        Records a log message

        Parameters:
        ----------
        *args : str
            The log message to record.

        separator : str, default: " "
            The separator to join the arguments.

        """
        log_str = separator.join([str(a) for a in args])
        self.logs.append(log_str)

    def exception(self, *args, separator: str = "\n") -> None:
        """
        Records an exception

        Parameters:
        ----------
        *args : str
            The exception to record.

        separator : str, default: "\n"
            The separator to join the arguments.

        """
        exception_str = separator.join([str(a) for a in args])
        self.exceptions.append(exception_str)

    def write(
        self, log_file: Optional[str] = None, sep: str = "\n", append: bool = False
    ):
        """
        Writes the log and exceptions to a file.

        Parameters:
        -----------
        log_file : str
            The path to the log file.

        append : bool, default: False
            If ``True``, the log will be appended to the file.

        sep : str, default: "\n"
            The separator to join the log messages.

        """
        if log_file is None and self.log_file is None:
            raise ValueError(
                "log_file must be provided, either at initialization or to the write() method"
            )
        log_file = log_file or self.log_file
        # logs
        self.checkpoint()
        # make sure the logfile's directory exists
        make_dir(os.path.dirname(log_file))
        # write to file
        mode = "a" if append else "w"
        with open(log_file, mode) as f:
            if append:
                f.write("\n")
            f.write(self.formatted)

    def checkpoint(self):
        """
        Checkpoints the log and exceptions by formatting them and adding to the formatted log string.

        Useful when iterating over a large number of items, and you want to save the log and exceptions
        for each item such that the log and exceptions are adjacent to each other in the final log file.
        """
        if self.logs:
            self.formatted += "\n".join(self.logs)
            self.logs = []
        if self.exceptions:
            self.formatted += "\n\nEXCEPTIONS\n"
            self.formatted += "\n".join(self.exceptions)
            self.exceptions = []


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
