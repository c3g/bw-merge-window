from importlib import metadata

__all__ = [
    "__version__",
]

__version__ = metadata.version(__name__)
