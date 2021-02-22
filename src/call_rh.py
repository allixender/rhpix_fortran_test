from ctypes import CDLL, POINTER, c_int, c_double, byref, pointer
import numpy as np

librhpix = CDLL("../build2/librhpix.so")
lam=-185.0
c_lam=c_double(-185.0)
radians=0

librhpix.c_wrap_longitude(byref(c_lam), byref(c_int(radians)))
print(c_lam.value)

phi=-135.0
librhpix.c_wrap_latitude.restype = c_double
res2 = librhpix.c_wrap_latitude(byref(c_double(phi)), byref(c_int(radians)))
print(res2)
