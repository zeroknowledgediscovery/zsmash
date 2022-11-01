import os
import warnings
import ctypes
import pandas as pd
import numpy as np
from cython import boundscheck, nonecheck, wraparound, cdivision
from cython cimport sizeof
from zsmash.utils import RANDOM_NAME, CythonWrapper
from zsmash.smash cimport _smash, _lsmash
from libc.stdlib cimport malloc, free
from libc.string cimport strcpy


class Smash(CythonWrapper):
    def __init__(
        self,
        helpmsg: bool = False,
        version: bool = False,
        data: [[]] = None,
        data_file: str = None,
        config_file: str = None,
        data_dir: str = None,
        data_type: str = None,
        seq_len: ctypes.c_uint = None,
        num_each: ctypes.c_uint = None,
        partition: list = None,
        timer: bool = None,
        outfile: str = None,
    ) -> None:

        super().__init__()

        self.argc = 1
        self.argv = ['smash']

        # Program information:
        if helpmsg:
            self._add_flag('-h')
        if version:
            self._add_flag('-V')

        # Usage:
        self.data = None
        if data is not None:
            if type(data) is pd.DataFrame or type(data) is pd.Series:
                self.data = data.values.tolist()
            else:
                self.data = data
        if data_file is not None:
            self._add_param('-f', data_file)
        if config_file is not None:
            self._add_param('-c', config_file)
        if data_dir is not None:
            self._add_param('-D', data_dir)
        if data_type is not None:
            if data_type.lower() != 'continuous' and data_type.lower() != 'symbolic':
                data_type = 'symbolic'
                warnings.warn("Invalid data type passed. Must be \"continuous\" or \"symbolic\". Assuming symbolic.", Warning)
            elif data_type.lower() == 'continuous':
                pass
            self._add_param('-T', data_type)
        if seq_len is not None:
            self._add_param('-L', str(seq_len))
        if num_each is not None:
            self._add_param('-n', str(num_each))
        if partition is not None:
            self._add_param('-p', ' '.join([str(n) for n in partition]))
        if timer is not None:
            self._add_param('-t', str(timer))
        if outfile is not None:
            self._add_param('-o', outfile)

        # Results
        self.dist_matrix = None

    @boundscheck(False)
    @nonecheck(False)
    @wraparound(False)
    @cdivision(True)
    def run(self):
        # convert self.argv to a char** to be passed to wrapped c++ function
        cdef char** argv = <char **> malloc(self.argc * sizeof(ctypes.c_char_p))
        for i, arg in enumerate(self.argv):
            argv[i] = <char *> malloc(len(arg) * sizeof(char))
            strcpy(argv[i], bytes(arg, 'utf-8'))

        # _smash() defined in src/_smash.cc
        self.dist_matrix = _smash(self.argc, argv, self.data)

        for i in range(self.argc):
            free(argv[i])
        free(argv)

        return self.dist_matrix


class Lsmash(CythonWrapper):
    def __init__(
        self,
        helpmsg: bool = False,
        version: bool = False,
        data: [[]] = None,
        seqfile: str = None,
        data_dir: str = None,
        data_type: str = None,
        seq_len: ctypes.c_uint = None,
        partition: list = None,
        derivative: bool = None,
        pfsafiles: list = None,
        timer: bool = None,
        sae: bool = None,
        repeat: ctypes.c_uint = None,
        outfile: str = None,
        random_mc: ctypes.c_uint = None,
        print_mc: bool = None,
    ) -> None:

        super().__init__()

        self.argc = 1
        self.argv = ['lsmash']

        # Program information:
        if helpmsg:
            self._add_flag('-h')
        if version:
            self._add_flag('-V')

        # Usage:
        self.data = None
        if data is not None:
            if type(data) is pd.DataFrame or type(data) is pd.Series:
                self.data = data.values.tolist()
            else:
                self.data = data
        if seqfile is not None:
            self._add_param('-f', seqfile)
        if data_dir is not None:
            self._add_param('-D', data_dir)
        if data_type is not None:
            self._add_param('-T', data_type)
        if seq_len is not None:
            self._add_param('-x', str(seq_len))
        if partition is not None:
            self._add_param('-P', ' '.join([str(n) for n in partition]))
        if derivative is not None:
            self._add_param('-u', str(derivative))
        if pfsafiles is not None:
            self._add_param('-F', ' '.join([str(n) for n in pfsafiles]))
        if timer is not None:
            self._add_param('-t', str(timer))
        if sae is not None:
            self._add_param('-S', str(sae))
        if repeat is not None:
            self._add_param('-n', str(repeat))
        if outfile is not None:
            self._add_param('-o', outfile)
        if random_mc is not None:
            self._add_param('-R', str(random_mc))
        if print_mc is not None:
            self._add_param('-m', str(print_mc))

        # Results
        self.dist_matrix = None

    @boundscheck(False)
    @nonecheck(False)
    @wraparound(False)
    @cdivision(True)
    def run(self):
        # convert self.argv to a char** to be passed to wrapped c++ function
        cdef char** argv = <char **> malloc(self.argc * sizeof(ctypes.c_char_p))
        for i, arg in enumerate(self.argv):
            argv[i] = <char *> malloc(len(arg) * sizeof(char))
            strcpy(argv[i], bytes(arg, 'utf-8'))

        # _lsmash() defined in _lsmash.cc
        self.dist_matrix = _lsmash(self.argc, argv, self.data)

        # testing parallelization without GIL
        # cdef unsigned int argc = self.argc
        # cdef matrix_dbl dist_matrix
        # openmp.omp_set_dynamic(1)
        # with nogil, parallel():
        #     dist_matrix = _lsmash(argc, argv)
        # self.dist_matrix = dist_matrix

        for i in range(self.argc):
            free(argv[i])
        free(argv)

        return self.dist_matrix

