import pandas as pd
import shutil
import uuid
import os
import atexit
import traceback
import glob
import inspect

class CythonWrapper:
    def __init__(self):
        self.temp_dir = "zed_temp"

    def _add_flag(self, flag: str) -> None:
        self.argv.append(flag)
        self.argc = len(self.argv)

    def _add_param(self, flag: str, value: str) -> None:
        self.argv.append(flag)
        self.argv.append(value)
        self.argc = len(self.argv)


def os_remove(filename):
    try:
        os.remove(filename)
    except OSError:
        pass


def run_once_per_process(f):
    def wrapper(*args, **kwargs):
        if os.getpgrp() != wrapper.has_run:
            wrapper.has_run = os.getpgrp()
            return f(*args, **kwargs)

    wrapper.has_run = 0
    return wrapper


def _clean_up_temp_folder(path, signal=None, frame=None):
    if signal is not None:
        traceback.print_stack(frame)
    try:
        [
            os_remove(x) if os.path.isfile(x) else shutil.rmtree(x)
            for x in glob.glob(path + "*")
        ]
        os.rmdir(os.path.dirname(path))
    except OSError:
        pass

    # print(frame.print_stack())
    # aulthandler.dump_traceback()
    # print(traceback)


def getValidKwargs(func, argsDict):

    args = set(inspect.getfullargspec(func).args)
    kwargs = set(inspect.getfullargspec(func).kwonlyargs)
    validKwargs = args.union(kwargs)

    validKwargs.discard('self')
    return dict((key, value) for key, value in argsDict.items()
                if key in validKwargs)


def callwithValidKwargs(func, argsDict):
    return func(**getValidKwargs(func, argsDict))


@run_once_per_process
def add_signal(path):
    atexit.register(_clean_up_temp_folder, path)
    # signal.signal(signal.SIGINT, partial(_clean_up_temp_folder, path))
    # signal.signal(signal.SIGTERM, partial(_clean_up_temp_folder, path))
    # faulthandler.register(signal.SIGINT)


def RANDOM_NAME(clean=True, path="zed_temp"):
    if not os.path.isdir(path):
        os.mkdir(path)
    path = os.path.join(path, "clean_")
    if clean:
        # path = path + str(os.getpgrp()) + "_"
        add_signal(path)
    random_name = str(uuid.uuid4())
    full = path + random_name
    while os.path.isfile(full):
        # print("name_double")
        random_name = str(uuid.uuid4())
        full = path + random_name

    return full
