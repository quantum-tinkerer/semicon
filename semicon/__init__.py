from . import models
from . import parameters
from . import peierls

from ._version import __version__


def test(verbose=True):
    from pytest import main
    import os.path

    return main(
        [os.path.dirname(os.path.abspath(__file__)), "-s"] + (["-v"] if verbose else [])
    )


test.__test__ = False

__all__ = ["parameters", "models", "peierls", "__version__"]
