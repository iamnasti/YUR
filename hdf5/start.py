from make_bn800 import *
from mg import *

def calculate_1():
    make_model_continuous()

def calculate_2():
    nameoflib = 'mgxs.h5'
    make_libgroup_macro(nameoflib)
    make_model_mg(nameoflib)

def calculate_3():
    nameoflib = 'libxs.h5'
    make_libgroup_micro(nameoflib)
    make_model_mg_by_nuclide(nameoflib)


calculate_1()