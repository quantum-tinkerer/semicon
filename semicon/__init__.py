__all__ = []

from ._version import __version__
__all__.append('__version__')

for module in ['parameters', 'models', 'peierls']:
    exec('from . import {0}'.format(module))
    __all__.append(module)


def test(verbose=True):
    from pytest import main
    import os.path

    return main([os.path.dirname(os.path.abspath(__file__)),
                 "-s"] + (['-v'] if verbose else []))

test.__test__ = False
