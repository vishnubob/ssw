import os
import itertools
from ctypes import *

def load_ssw_library():
    mod_path = os.path.split(__file__)[0]
    search_path = os.path.join(mod_path, "..")
    libname = "_libssw"
    for (root, dirs, files) in os.walk(search_path):
        for fn in files:
            if fn.lower().startswith(libname) and fn.lower().endswith(".so"):
                libpath = os.path.join(root, fn)
                return cdll.LoadLibrary(libpath)
    return cdll.LoadLibrary(libname)

# load our library
libssw = load_ssw_library()

# Types
class AlignmentResult(Structure):
    _fields_ = [
        ("score", c_uint16),
        ("score2", c_uint16),
        ("ref_begin", c_int32),
        ("ref_end", c_int32),
        ('query_begin', c_int32),
        ('query_end', c_int32),
        ('ref_end2', c_int32),
        ('cigar', POINTER(c_uint32)),
        ('cigarLen', c_int32),
    ]

ssw_profile_p = c_void_p
matrix_type = c_int8
symbol_type = c_int8

# ssw_init method
ssw_profile_init = libssw.ssw_init
ssw_profile_init.argtypes = [POINTER(c_int8), c_int32, POINTER(c_int8), c_int32, c_int8]
ssw_profile_init.restype = ssw_profile_p

# init_destroy function
ssw_profile_del = libssw.init_destroy
ssw_profile_del.argtypes =  [ssw_profile_p]
ssw_profile_del.restype = None

# ssw_align function
ssw_align_init = libssw.ssw_align
ssw_align_init.argtypes = [ssw_profile_p, POINTER(c_int8), c_int32, c_uint8, c_uint8, c_uint8, c_uint16, c_int32, c_int32]
ssw_align_init.restype = POINTER(AlignmentResult)

# align_destroy function
ssw_align_del = libssw.align_destroy
ssw_align_del.argtypes = [POINTER(AlignmentResult)]
ssw_align_del.restype = None

# cigar_int_to_len function
cigar_int_to_len = libssw.cigar_int_to_len
cigar_int_to_len.argtypes = [c_int32]
cigar_int_to_len.restype = c_int32

# cigar_int_to_op function
cigar_int_to_op = libssw.cigar_int_to_op
cigar_int_to_op.argtypes = [c_int32]
cigar_int_to_op.restype = c_char

# flags
FLAG_BEST_POS = 0
FLAG_FILTER_SCORE = 1
FLAG_FILTER_DISTANCE = 2
FLAG_BUILD_CIGAR = 3
