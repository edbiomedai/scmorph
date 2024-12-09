import logging
import sys
from typing import TextIO

LOG_LEVEL = logging.INFO
LOG_HANDLE = "scmorph"


def set_logger(stream: TextIO = sys.stdout) -> None:
    """
    Set logger handle

    Parameters
    ----------
    stream
        Output stream
    """
    root = logging.getLogger(LOG_HANDLE)
    root.setLevel(LOG_LEVEL)

    handler = logging.StreamHandler(stream)
    handler.setLevel(LOG_LEVEL)
    formatter = logging.Formatter("%(levelname)s: %(message)s")
    handler.setFormatter(formatter)
    root.addHandler(handler)


def get_logger() -> logging.Logger:
    """
    Get package logger handle

    Returns
    -------
    logging.Logger
        Handle to logger
    """
    return logging.getLogger(LOG_HANDLE)


set_logger()
