import os
from sys import argv
import numpy as np
import h5py
import xml.dom.minidom

PATH_ASSEMBLY = "/Reactor#000000001/Reactor/Assembly00000"
PATH_HEADER = "/Reactor#000000001/Reactor/HEADER"
PATH_DATA = "/Assembly/Rod0/RodTVEL/NC/MatrixCF/Data"
PATH_HEADO = "/Assembly/HEAD0"
PATH_MATRIX = "/Assembly/Rod0/RodTVEL/CMC/SectionPEL/CF000000001/MatrixCF/Data"


# т.к может быть разное к-во сборок
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
        else:
            heado = PATH_ASSEMBLY + '{:04}'.format(i) + PATH_HEADO
            arr.append(heado)
    return arr


def xml_f(path):
    xmlf = xml.dom.minidom.parse(path)
    namenuclides = xmlf.getElementsByTagName("namenuclides")[0]
    return ''.join(namenuclides.firstChild.data).split()


def run(path):
    # path_assembly - массив с путями до всех массивов с нуклидами
    path_assembly = for_all_dir(1, path)
    # path_el_matrix - массив с путями до массивов элементов
    path_el_matrix = for_all_dir(2, path)
    # path_heado - массив с индикаторами о сузах
    path_heado = for_all_dir(3, path)

    with h5py.File(path, 'r') as f:
        diction_1 = dict()
        diction_2 = dict()
        diction_3 = dict()
        for i in range(1, len(path_heado) + 1):
            date_new_1 = np.transpose(f[path_assembly[i - 1]])
            date_new_2 = np.transpose(f[path_el_matrix[i - 1]])
            if f[path_heado[i - 1]].attrs['is Assembly'][0] == 1:
                # 1 словарь
                for num, val in enumerate(date_new_1):
                    t = (i, num + 1)
                    diction_1[t] = val
                # 2 словарь
                for num, val in enumerate(date_new_2):
                    t = (i, num + 1)
                    diction_2[t] = val
                    # 3 словарь
                    diction_3[t] = []
            else:
                # 1 словарь
                for num, val in enumerate(date_new_1):
                    t = (i, num + 1)
                    diction_1[t] = []
                # 2 словарь
                for num, val in enumerate(date_new_2):
                    t = (i, num)
                    diction_2[t] = []
                    # 3 словарь
                    diction_3[t] = val
    check(diction_1, "check_1.txt")
    check(diction_2, "check_2.txt")
    check(diction_3, "check_3.txt")


# для проверки
def check(diction, name):
    with open(name, 'w') as out:
        if type(diction) == dict:
            for key, val in diction.items():
                out.write('{}:{}\n'.format(key, val))
        else:
            for num, val in enumerate(diction):
                out.write('{}:{}\n'.format(num, val))


if __name__ == "__main__":
    n_arg = argv[1]
    m_arg = argv[2]
    nuclide = xml_f(m_arg)
    check(nuclide, "check_nuclide.txt")
    run(n_arg)

