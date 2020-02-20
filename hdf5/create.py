import os
from sys import argv
import numpy as np
import h5py

PATH_ASSEMBLY = "/Reactor#000000001/Reactor/Assembly00000"
PATH_DATA = "/Assembly/Rod0/RodTVEL/NC/MatrixCF/Data"
PATH_HEADO = "/Assembly/HEAD0"

def for_all_dir(number):
    arr = []
    for i in range(1,967):
        if number == 1:
            data = PATH_ASSEMBLY + str('{:04}'.format(i)) + PATH_DATA
            arr.append(data)
        else:
            heado = PATH_ASSEMBLY + str('{:04}'.format(i)) + PATH_HEADO
            arr.append(heado)
    return arr


def run(PATH):
    # path_assembly - массив с путями до всех массивов с нуклидами
    path_assembly = for_all_dir(1)
    # path_heado - массив с индикаторами о сузах
    path_heado = for_all_dir(2)
    with h5py.File(PATH, 'r') as f:
        diction = dict()
        for i in range(1, len(path_heado)+1):
            if f[path_heado[i-1]].attrs['is Assembly'][0] == 1:
                date_new = np.transpose(f[path_assembly[i-1]])
                for j in range(len(date_new)):
                    t = (i,j)
                    diction[t] = date_new[j]
            else:
                date_new = np.transpose(f[path_assembly[i-1]])
                for j in range(1,len(date_new)+1):
                    t = (i,j)
                    diction[t] = []
    check(diction)

#для проверки
def check(diction):
    with open('check.txt', 'w') as out:
        for key, val in diction.items():
            out.write('{}:{}\n'.format(key, val))


if __name__ == "__main__":
    narg = argv[1]
    run(narg)








