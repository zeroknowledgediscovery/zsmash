# details here: https://stackoverflow.com/a/52729181/8728279
from importlib.abc import MetaPathFinder
from importlib.machinery import ExtensionFileLoader
import sys

# Chooses the right init function
class CythonPackageMetaPathFinder(MetaPathFinder):
    def __init__(self, name_filter):
        super(CythonPackageMetaPathFinder, self).__init__()
        self.name_filter =  name_filter

    def find_module(self, fullname, path):
        if fullname.startswith(self.name_filter):
            # use this extension-file but PyInit-function of another module:
            return ExtensionFileLoader(fullname,__file__)


# injecting custom finder/loaders into sys.meta_path:
def bootstrap_cython_submodules():
    sys.meta_path.append(CythonPackageMetaPathFinder('zsmash.'))
