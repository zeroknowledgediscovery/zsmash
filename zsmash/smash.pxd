from libcpp.vector cimport vector
from libcpp.map cimport map

ctypedef map[unsigned int, map[unsigned int, double]] matrix_dbl

cdef extern from "src/_smash.cc":
    matrix_dbl _smash(int argc, char *argv[], vector[vector[unsigned int]] &data_in)

cdef extern from "src/_lsmash.cc":
    matrix_dbl _lsmash(int argc, char *argv[], vector[vector[unsigned int]] &data_in)
