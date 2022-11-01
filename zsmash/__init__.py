# details here: https://stackoverflow.com/a/52714500/8728279
# bootstrap is the only module which
# can be loaded with default Python-machinery
# because the resulting extension is called `bootstrap`:
from . import bootstrap

# injecting our finders into sys.meta_path
# after that all other submodules can be loaded
bootstrap.bootstrap_cython_submodules()