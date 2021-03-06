import os
import openmc
import time
from Table_Mend import *
from sys import argv
import numpy as np
import h5py
import xml.dom.minidom

PATH_ASSEMBLY = "/Reactor#000000001/Reactor/Assembly00000"
PATH_HEADER = "/Reactor#000000001/Reactor/HEADER"
PATH_DATA = "/Assembly/Rod0/RodTVEL/NC/MatrixCF/Data"
PATH_HEADO = "/Assembly/HEAD0"
PATH_MATRIX = "/Assembly/Rod0/RodTVEL/CMC/SectionPEL/CF000000001/MatrixCF/Data"
PATH_NAME_EL = "/Reactor#000000001/Reactor/Assembly000000001/Assembly/Rod0/RodTVEL/CMC/SectionPEL/NuclideIndex"
PATH_NATR = '/Assembly/PEL000000005/SectionPEL/CF000000002/MatrixCF/Data'
PATH_TYPE_FUEL = '/Assembly/HeigthFuel/MatrixCF/Data'
PATH_TYPE_ASSEMBLY = '/Assembly/HEADER'

# different number of assembly
def deter_const(path):
    with h5py.File(path, 'r') as f:
        return f[PATH_HEADER].attrs['Amount Of Assemblies'][0]


def for_all_dir(number, path):
    arr = []
    for i in range(1, deter_const(path) + 1):
        if number == 1:
            data = PATH_ASSEMBLY + '{:04}'.format(i) + PATH_DATA
            arr.append(data)
        elif number == 2:
            el = PATH_ASSEMBLY + '{:04}'.format(i) + PATH_MATRIX
            arr.append(el)
        elif number == 3:
            el = PATH_ASSEMBLY + '{:04}'.format(i) + PATH_NATR
            arr.append(el)
        elif number == 4:
            el = PATH_ASSEMBLY + '{:04}'.format(i) + PATH_TYPE_FUEL
            arr.append(el)
        elif number == 5:
            el = PATH_ASSEMBLY + '{:04}'.format(i) + PATH_TYPE_ASSEMBLY
            arr.append(el)
        else:
            heado = PATH_ASSEMBLY + '{:04}'.format(i) + PATH_HEADO
            arr.append(heado)
    return arr


def xml_f(path):
    xmlf = xml.dom.minidom.parse(path)
    namenuclides = xmlf.getElementsByTagName("namenuclides")[0]
    return ''.join(namenuclides.firstChild.data).split()


def name_el(path):
    arr = []
    with h5py.File(path, 'r') as f:
        for num, val in enumerate(f[PATH_NAME_EL]):
            # except B because we have nuclides
            if str(val) == '50100':
                arr.append("B10")
            elif str(val) == '50110':
                arr.append("B11")
            elif val != 0:
                if len(str(val)) == 6:
                    arr.append(determ_el(str(val)[0:5]))
                else:
                    arr.append(determ_el(str(val)[0:4]))
    return arr


def run(path):

    path_assembly = for_all_dir(1, path)
    path_el_matrix = for_all_dir(2, path)
    path_natr = for_all_dir(3,path)
    path_type_fuel = for_all_dir(4, path)
    path_type_assembly = for_all_dir(5, path)
    path_heado = for_all_dir(10, path)

    dzz = []; hss = []; rods = []
    with h5py.File(path, 'r') as f:

        diction_1 = dict()
        diction_2 = dict()
        diction_3 = dict()
        densna = dict()
        diction_5 = dict()
        diction_6 = dict()

        print(len(path_heado))
        for i in range(1, len(path_heado)):

            date_new_1 = np.transpose(f[path_assembly[i - 1]])
            date_new_2 = np.transpose(f[path_el_matrix[i - 1]])
            date_new_3 = np.array(f[path_natr[i-1]])
            date_new_3 = np.array(f[path_natr[i - 1]])
            date_new_4 = np.array(f[path_type_fuel[i - 1]])

            densna[i] = date_new_3[0]
            diction_5[i] = date_new_4[1]
            diction_6[i] = f[path_type_assembly[i - 1]].attrs['Type Of Assembly'][0]
            g = f[PATH_ASSEMBLY + '{:04}'.format(i)]

            if f[path_heado[i - 1]].attrs['is Assembly'][0] == 1:
                dzz = g['Assembly']['HeigthFuel']['MatrixCF']['Data'][0]
                rods.append(False)
                for num, val in enumerate(date_new_1):
                    t = (i, num + 1)
                    diction_1[t] = val
                for num, val in enumerate(date_new_2):
                    t = (i, num + 1)
                    # cut our array because we have only 25 elements
                    diction_2[t] = val[0:25]
            else:
                hss = g['Assembly']['Zones']['MatrixCF']['Data'][0]
                rods.append(True)
                for num, val in enumerate(date_new_1):
                    t = (i, num + 1)
                for num, val in enumerate(date_new_2):
                    t = (i, num)
                    diction_3[t] = val[0:25]
    # check(diction_1, "check_1.txt")
    # check(diction_2, "check_2.txt")
    # check(diction_5, "check_4.txt")
    # check(diction_6, "check_5.txt")
    return diction_1, diction_2, diction_3, dzz, hss, rods, densna, diction_5, diction_6


def check(diction, name):
    with open(name, 'w') as out:
        if type(diction) == dict:
            for key, val in diction.items():
                out.write('{}:{}\n'.format(key, val))
        else:
            for num, val in enumerate(diction):
                out.write('{}:{}\n'.format(num, val))
# my path
# nuclide = xml_f("/home/barakuda/PycharmProjects/hdf5/reactions.xml")
# element = name_el("/home/barakuda/yur/archive/800/InOut6800.dat.hdf5")

