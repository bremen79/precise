from ctypes import c_uint, c_int, c_uint32, c_uint64, c_float, c_double, c_char_p, c_void_p, CFUNCTYPE, c_bool, Structure, POINTER, byref
import numpy as np
import numpy.ctypeslib as npct

arraylf = npct.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
array2f = npct.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS')

lib = npct.load_library("lib_opt.so", ".")

lib.find_mean_max_log_wealth_constrained.restype = c_float
lib.find_mean_max_log_wealth_constrained.argtypes= [arraylf, c_uint32, c_float, c_float, POINTER(c_float)]


def find_mean_max_log_wealth_constrained(
    X: np.ndarray,
    bmin: float,
    bmax: float,
    m_try: float
):
    x_out = c_float()
    
    X_ = X - m_try
    log_wealth = lib.find_mean_max_log_wealth_constrained(
        X_.T.astype(np.float32, order='C'),
        c_uint32(X.size),
        c_float(bmin),
        c_float(bmax),
        byref(x_out),
    )
    return (log_wealth, x_out)
