#!/usr/bin/env python
# filename: pipeline.py


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


import glob
import os
import re
from typing import Iterable, Optional, Union

from natsort import natsorted

# =======================
#   General file I/O
# =======================


def make_dir(directory: str) -> None:
    """
    Makes a directory, if it doesn't already exist.

    Parameters
    ----------
    directory : str
        Path to a directory.

    """
    if not os.path.exists(directory):
        os.makedirs(os.path.abspath(directory))


def list_files(
    directory: str,
    extension: Union[str, Iterable, None] = None,
    recursive: bool = False,
    match: Optional[str] = None,
    ignore_dot_files: bool = True,
) -> Iterable[str]:
    """
    Lists files in a given directory.

    Parameters
    ----------
    directory : str
        Path to a directory. If a file path is passed instead, the returned list of files will contain
        only that file path.

    extension : str
        If supplied, only files that end with the specificied extension(s) will be returned. Can be either
        a string or a list of strings. Extension evaluation is case-insensitive and can match complex
        extensions (e.g. '.fastq.gz'). Default is ``None``, which returns all files in the directory,
        regardless of extension.

    recursive : bool, default=False
        If ``True``, the directory will be searched recursively, and all files in all subdirectories will be returned.

    match : str, optional
        If supplied, only files that match the specified pattern will be returned. Regular expressions are supported.

    ignore_dot_files : bool, default=True
        If ``True``, dot files (hidden files) will be ignored.

    Returns
    -------
    Iterable[str]

    """
    directory = os.path.abspath(directory)
    if os.path.exists(directory):
        if os.path.isdir(directory):
            if recursive:
                files = []
                for root, dirs, _files in os.walk(directory):
                    for f in _files:
                        file_path = os.path.join(root, f)
                        files.append(file_path)
            else:
                files = natsorted(glob.glob(directory + "/*"))
        else:
            files = [directory]
    else:
        raise ValueError(f"Directory {directory} does not exist.")
    if extension is not None:
        if isinstance(extension, str):
            extension = [
                extension,
            ]
        files = [
            f
            for f in files
            if any(
                [
                    any([f.lower().endswith(e.lower()) for e in extension]),
                    any([f.endswith(e.upper()) for e in extension]),
                    any([f.endswith(e.lower()) for e in extension]),
                ]
            )
        ]
    if ignore_dot_files:
        files = [f for f in files if not os.path.basename(f).startswith(".")]
    if match is not None:
        files = [f for f in files if re.match(match, os.path.basename(f)) is not None]
    return files


def rename_file(file: str, new_name: str) -> None:
    """
    Renames a file.

    Parameters
    ----------
    file : str
        Path to the file to be renamed.

    new_name : str
        New name for the file.

    """
    os.rename(file, new_name)


def delete_files(files: Union[str, Iterable]) -> None:
    """
    Deletes files.

    Parameters
    ----------
    files : Union[str, Iterable]
        Path to a file or an iterable of file paths.

    """
    if isinstance(files, str):
        files = [
            files,
        ]
    for f in files:
        if os.path.exists(f):
            os.remove(f)


# from __future__ import absolute_import, division, print_function, unicode_literals

# import glob
# import os
# import sys
# from typing import Iterable, Optional, Union

# from . import log


# # for backward compatibility
# def list_files(*args, **kwargs):
#     from ..io import list_files

#     return list_files(*args, **kwargs)


# def make_dir(*args, **kwargs):
#     from ..io import make_dir

#     return make_dir(*args, **kwargs)


# def initialize(log_file, project_dir=None, debug=False):
#     """
#     Initializes an AbTools pipeline.

#     Initialization includes printing the AbTools splash, setting up logging,
#     creating the project directory, and logging both the project directory
#     and the log location.

#     Parameters
#     ----------
#     log_file : str
#         Path to the log file. Required.

#     project_dir : str
#         Path to the project directory. If not provided,
#         the project directory won't be created and the location won't be logged.

#     debug : bool
#         If ``True``, the logging level will be set to ``logging.DEBUG``.

#     Returns
#     -------
#     logger
#         A logger instance.

#     """
#     print_splash()
#     log.setup_logging(log_file, print_log_location=False, debug=debug)
#     logger = log.get_logger("pipeline")
#     if project_dir is not None:
#         make_dir(os.path.normpath(project_dir))
#         logger.info("PROJECT DIRECTORY: {}".format(project_dir))
#         logger.info("")
#     logger.info("LOG LOCATION: {}".format(log_file))
#     print("")
#     return logger


# # def make_dir(directory: str) -> None:
# #     """
# #     Makes a directory, if it doesn't already exist.

# #     Parameters
# #     ----------
# #     directory : str
# #         Path to a directory.

# #     """
# #     if not os.path.exists(directory):
# #         os.makedirs(os.path.abspath(directory))


# # def list_files(
# #     directory: str, extension: Union[str, Iterable, None] = None
# # ) -> Iterable[str]:
# #     """
# #     Lists files in a given directory.

# #     Parameters
# #     ----------
# #     directory : str
# #         Path to a directory.

# #     extension : str
# #         If supplied, only files that end with the specificied extension(s) will be returned. Can be either
# #         a string or a list of strings. Extension evaluation is case-insensitive and can match complex
# #         extensions (e.g. '.fastq.gz'). Default is ``None``, which returns all files in the directory,
# #         regardless of extension.

# #     Returns
# #     -------
# #     Iterable[str]

# #     """
# #     if os.path.isdir(directory):
# #         expanded_dir = os.path.expanduser(directory)
# #         files = sorted(glob.glob(expanded_dir + "/*"))
# #     else:
# #         files = [
# #             directory,
# #         ]
# #     if extension is not None:
# #         if isinstance(extension, str):
# #             extension = [
# #                 extension,
# #             ]
# #         files = [
# #             f
# #             for f in files
# #             if any(
# #                 [
# #                     any([f.lower().endswith(e.lower()) for e in extension]),
# #                     any([f.endswith(e.upper()) for e in extension]),
# #                     any([f.endswith(e.lower()) for e in extension]),
# #                 ]
# #             )
# #         ]
# #     return files


# def print_splash():
#     splash = """
#      _    _     _   _ _   _ _       ____  _            _ _
#     / \  | |__ | | | | |_(_) |___  |  _ \(_)_ __   ___| (_)_ __   ___
#    / _ \ | '_ \| | | | __| | / __| | |_) | | '_ \ / _ \ | | '_ \ / _ \\
#   / ___ \| |_) | |_| | |_| | \__ \ |  __/| | |_) |  __/ | | | | |  __/
#  /_/   \_\_.__/ \___/ \__|_|_|___/ |_|   |_| .__/ \___|_|_|_| |_|\___|
#                                            |_|
#     """
#     print("")
#     print(splash)
#     print(
#         "(c) 2023 Bryan Briney\nDistributed under the MIT License (http://opensource.org/licenses/MIT)"
#     )
#     print("")
