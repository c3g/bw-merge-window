import functools
import logging

__all__ = ["get_logger"]


@functools.cache
def get_logger() -> logging.Logger:
    logger = logging.getLogger("bw-merge-window")
    logger.setLevel(logging.INFO)
    return logger
