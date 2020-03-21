import os
import sys
import re
import numpy as np
import openmc
import openmc.mgxs
from make_bn800 import *
# Start global variables
N_LEGEND = 2
N_DELAY = 1
# PATH
PATH = "/home/barakuda/Рабочий стол/hdf5_openmc/XMAS172jeff2p2woZr71"
# PATH
ENERGIES = [19640330,
17332530,
14918250,
13840310,
11618340,
10000000,
8187308,
6703200,
6065307,
5488116,
4493290,
3678794,
3011942,
2465970,
2231302,
2018965,
1652989,
1353353,
1224564,
1108032,
1002588,
907179.5,
820850,
608100.625,
550232.188,
497870.688,
450492,
407622,
301973.812,
273237.188,
247235.297,
183156.406,
122773.398,
111090,
82297.4688,
67379.4688,
55165.6406,
40867.7109,
36978.6406,
29283,
27394.4492,
24787.5195,
16615.5703,
15034.3896,
11137.75,
9118.82031,
7465.85791,
5530.84424,
5004.51416,
3526.62207,
3354.62598,
2248.6731,
2034.68396,
1507.33105,
1433.81702,
1234.09802,
1010.39398,
914.24231,
748.518311,
677.287415,
453.999298,
371.703186,
304.324799,
203.994995,
148.625397,
136.742004,
91.660881,
75.6735687,
67.9040527,
55.5951309,
51.5780182,
48.2515984,
45.5174408,
40.1689987,
37.2665291,
33.72015,
30.5112591,
27.6077309,
24.9804993,
22.6032906,
19.4548397,
15.9282703,
13.70959,
11.2244596,
9.90555382,
9.18981361,
8.31528664,
7.523983,
6.1601162,
5.34642982,
5.04347706,
4.12925005,
4,
3.38074994,
3.29999995,
2.76792002,
2.72000003,
2.5999999,
2.54999995,
2.3599999,
2.13000011,
2.0999999,
2.01999998,
1.92999995,
1.84000003,
1.755,
1.66999996,
1.59000003,
1.5,
1.47500002,
1.44000006,
1.37,
1.33749998,
1.29999995,
1.23500001,
1.16999996,
1.14999998,
1.12300003,
1.11000001,
1.097,
1.07099998,
1.04499996,
1.03499997,
1.01999998,
0.995999992,
0.986000001,
0.972000003,
0.949999988,
0.930000007,
0.910000026,
0.860000014,
0.850000024,
0.790000021,
0.779999971,
0.704999983,
0.625,
0.540000021,
0.5,
0.485000014,
0.432999998,
0.400000006,
0.391000003,
0.349999994,
0.319999993,
0.314500004,
0.300000012,
0.280000001,
0.247999996,
0.219999999,
0.188999996,
0.180000007,
0.159999996,
0.140000001,
0.134000003,
0.115000002,
0.100000001,
9.50E-02,
8.00E-02,
7.70E-02,
6.70E-02,
5.80E-02,
5.00E-02,
4.20E-02,
3.50E-02,
3.00E-02,
2.50E-02,
2.00E-02,
1.50E-02,
1.00E-02,
6.90E-03,
5.00E-03,
3.00E-03,
1.00E-05]
tempdata = [300]
# tempdata = [300, 600, 900, 1200, 1500, 1800, 2100]
libgroup = {t : {} for t in tempdata}
ENERGIES = ENERGIES[::-1]
# Instantiate the energy group data
groups = openmc.mgxs.EnergyGroups(np.array(ENERGIES))
# end global variables
class groupsection:
    def __init__(self, name, egroup, nlegendr, numdel, temp):
        self.csname = name
        self.energy = egroup
        self.ng = len(self.energy) - 1
        self.order = nlegendr
        self.numdel = numdel
        self.temperature = temp
    def read_section(self, path):
        f = open(path, 'r')
        self._is_finished = False
        self.nusfdel = []
        self.chidel = []
        self.scatxs = {i : [] for i in range(self.order)}
        self.xsmatrix = []
        _matxs = {i: [] for i in range(self.order)}
        _curlegendr = 0
        for line in f:
            if (re.search("\sSF\s", line)):
                self.sf = self._read_cs_data(f)
            if (re.search("STOT", line)):
                self.stot = self._read_cs_data(f)
            if (re.search("SABS", line)):
                self.sabs = self._read_cs_data(f)
            if (re.search("SCAPT", line)):
                self.scapt = self._read_cs_data(f)
            if (re.search("CHI0", line)):
                self.schi = self._read_cs_data(f)
            if (re.search("NUSF0", line)):
                self.nusf = self._read_cs_data(f)
            if (re.search("CHIDEL", line)):
                self.chidel.append(self._read_cs_data(f))
            if (re.search("NUSFDEL", line)):
                self.nusfdel.append(self._read_cs_data(f))
            if (re.search("SIGS", line)):
                _curlegendr = int(line.split()[-1][-1])
                self.scatxs[_curlegendr].append(self._read_scat_data(f))
            if (re.search("\sSCAT\d+", line)):
                _curlegendr = int(line.split()[-1][-1])
            if (re.search("SCATTERING FROM GROUP", line)):
                _matxs[_curlegendr].append(self._read_scat_data(f))
        for k, v in _matxs.items():
            if (len(v) > 0):
                self.xsmatrix.append(np.array(v).reshape(self.ng, self.ng))
        self.xsmatrix = np.array(self.xsmatrix)
        if (len(self.xsmatrix.shape) > 2):
            self.xsmatrix = np.rollaxis(self.xsmatrix, 0 , 3)
        f.close()
    def _read_cs_data(self, f):
        arr = []
        for e in self.energy[1:]:
            line = next(f)
            arr.append(float(line.strip()))
        return np.array(arr)
    def _read_scat_data(self, f):
        arr = []
        while (len(arr) < self.ng):
            line = next(f)
            arr.extend([float(s) for s in line.split()])
        return np.array(arr)
### Go into Function!!!
def temp_libgroup(tempdata):
    for t in tempdata:
        tree = [tt for tt in os.walk(os.path.join(PATH, str(t)))]
        #print(len( tree))
        for fn in tree[0][-1]:
            libgroup[t][fn] = groupsection(fn, ENERGIES, N_LEGEND, N_DELAY, t)
            libgroup[t][fn].read_section(os.path.join(PATH, str(t), fn))
            #print("Temp {} : {} ".format(t, fn))
    return libgroup
# return libgroup
### Go into Function!!!
"""
conc = {"fuel" : [("U238" , 1.8744e-2), ("O16" , 0.039235), ("U235" , 8.737e-4)],
        "clad" : [("Zr0" , 0.0423)],
        "water" : [("H1_H2O" , 0.06694), ("O16" , 0.03347), ("B10" , 6.6262e-6), ("B11" , 2.6839E-5)]}

conc = {"fuel" : [("U238" , 1.8744e-2), ("O16" , 0.039235), ("U235" , 8.737e-4), ("Xe135" , 9.4581E-9), ("Sm149" , 7.3667E-8)],
        "clad" : [("Zr0" , 0.0423)],
        "water" : [("H1_H2O" , 0.04783), ("O16" , 0.02391), ("B10" , 4.7344e-6), ("B11" , 1.9177E-5)]}
mattemp = {"fuel" : 1027, "clad" : 579, "water" : 579}
"""


def temp_interolate(temparray, T):
    if (T <= temparray[0]):
        return [(temparray[0], 1)]
    elif (T >= temparray[-1]):
        return [(temparray[-1], 1)]
    else:
        for j,_t in enumerate(temparray):
            if (_t > T):
                i = j-1
                t = temparray[i]
                break
        return [(t, (temparray[i+1] -T)/(temparray[i+1] - temparray[i])), (temparray[i+1], (T - t)/(temparray[i+1] - temparray[i]))]



# TEMPERATURE INDEPENDENSE CASE
"""
for name, val in conc.items():
    openmclib[name] = openmc.XSdata(name, groups)
    openmclib[name].order = 1
    stot = np.zeros(len(ENERGIES) - 1, dtype=np.double)
    sabs = np.zeros(len(ENERGIES) - 1, dtype=np.double)
    scapt = np.zeros(len(ENERGIES) - 1, dtype=np.double)
    sf = np.zeros(len(ENERGIES) - 1, dtype=np.double)
    nusf = np.zeros(len(ENERGIES) - 1, dtype=np.double)
    chi = np.zeros(len(ENERGIES) - 1, dtype=np.double)
    scatter = np.zeros((len(ENERGIES) - 1, len(ENERGIES) - 1, 2), dtype = np.double)
    concentration = 0.0
    for v in val:
        for el in temp_interolate(tempdata, mattemp[name]):
            tt = el[0]; wt = el[1]
            stot += libgroup[tt][v[0]].stot * v[1] * wt
            sabs += libgroup[tt][v[0]].sabs * v[1] * wt
            scapt += libgroup[tt][v[0]].scapt * v[1] * wt
            sf += libgroup[tt][v[0]].sf * v[1] * wt
            nusf += libgroup[tt][v[0]].nusf * v[1] * wt
            scatter += libgroup[tt][v[0]].xsmatrix * v[1] * wt
        if (libgroup[300][v[0]].sf.sum() > 0):
            concentration += v[1]
            chi += libgroup[300][v[0]].schi * v[1]
    if (concentration > 0):
        chi = chi/concentration
    openmclib[name].set_total(stot, temperature=294.)
    openmclib[name].set_absorption(sabs, temperature=294.)
    openmclib[name].set_scatter_matrix(scatter, temperature=294.)
    openmclib[name].set_fission(sf, temperature=294.)
    openmclib[name].set_nu_fission(nusf, temperature=294.)
    openmclib[name].set_chi(chi, temperature=294.)

# TEMPERATURE DEPENDET CASE
for name, val in conc.items():
    openmclib[name] = openmc.XSdata(name, groups, temperatures=tempdata)
    openmclib[name].order = 1
    for tt in tempdata:
        stot = np.zeros(len(ENERGIES) - 1, dtype=np.double)
        sabs = np.zeros(len(ENERGIES) - 1, dtype=np.double)
        scapt = np.zeros(len(ENERGIES) - 1, dtype=np.double)
        sf = np.zeros(len(ENERGIES) - 1, dtype=np.double)
        nusf = np.zeros(len(ENERGIES) - 1, dtype=np.double)
        chi = np.zeros(len(ENERGIES) - 1, dtype=np.double)
        scatter = np.zeros((len(ENERGIES) - 1, len(ENERGIES) - 1, 2), dtype=np.double)
        concentration = 0.0
        for v in val:
            stot += libgroup[tt][v[0]].stot * v[1]
            sabs += libgroup[tt][v[0]].sabs * v[1]
            scapt += libgroup[tt][v[0]].scapt * v[1]
            sf += libgroup[tt][v[0]].sf * v[1]
            nusf += libgroup[tt][v[0]].nusf * v[1]
            scatter += libgroup[tt][v[0]].xsmatrix * v[1]
            if (libgroup[tt][v[0]].sf.sum() > 0):
                concentration += v[1]
                chi += libgroup[tt][v[0]].schi * v[1]
        if (concentration > 0):
            chi = chi / concentration
        openmclib[name].set_total(stot, temperature=tt)
        openmclib[name].set_absorption(sabs, temperature=tt)
        openmclib[name].set_scatter_matrix(scatter, temperature=tt)
        openmclib[name].set_fission(sf, temperature=tt)
        openmclib[name].set_nu_fission(nusf, temperature=tt)
        openmclib[name].set_chi(chi, temperature=tt)
"""
def prepare_temperature_independed_mg(libgroup, conc, namenuclide, tempdata,
                                                                     groups, mattemp):
    """
    Prepare multigroup data for calculation based with temperature independent
    constant
    Paramertres:
    -----------
    libgroup : dictionary
    {temperature : dictonary { name of nuclide : element groupsection class }};
    :param conc: dict
    - dictionary with name of material : np.array - R*8 concentration of nuclide;
    :param namenuclide:
    - a list of nuclide names;
    :param temperature:
    - dictionary with name of material : temperature value;
    groups - global variable groups
    mattemp - dictionary with name of material : temperature
    :return:
    - dictionary with name of material : openmc.MGXS class element
    """
    # TEMPERATURE INDEPENDET CASE
    openmclib = {}
    for name, val in conc.items():
        openmclib[name] = openmc.XSdata(name, groups, temperatures=[tempdata[0]])
        openmclib[name].order = 1
        stot = np.zeros(len(ENERGIES) - 1, dtype=np.double)
        sabs = np.zeros(len(ENERGIES) - 1, dtype=np.double)
        scapt = np.zeros(len(ENERGIES) - 1, dtype=np.double)
        sf = np.zeros(len(ENERGIES) - 1, dtype=np.double)
        nusf = np.zeros(len(ENERGIES) - 1, dtype=np.double)
        chi = np.zeros(len(ENERGIES) - 1, dtype=np.double)
        scatter = np.zeros((len(ENERGIES) - 1, len(ENERGIES) - 1, 2), dtype=np.double)
        concentration = 0.0
        for n, v in zip(namenuclide, val):
            for el in temp_interolate(tempdata, mattemp[name]):
                tt = el[0]
                wt = el[1]
                if (n in libgroup[tt].keys()):
                    print("Nuclide ",n, tt)
                    stot += libgroup[tt][n].stot * v * wt
                    sabs += libgroup[tt][n].sabs * v * wt
                    scapt += libgroup[tt][n].scapt * v * wt
                    sf += libgroup[tt][n].sf * v * wt
                    nusf += libgroup[tt][n].nusf * v * wt
                    scatter += libgroup[tt][n].xsmatrix * v * wt
            if (n in libgroup[tempdata[0]].keys()):
                if (libgroup[tempdata[0]][n].sf.sum() > 0):
                    concentration += v
                    chi += libgroup[tempdata[0]][n].schi * v
        if (concentration > 0):
            chi = chi / concentration
        openmclib[name].set_total(stot, temperature=tempdata[0])
        openmclib[name].set_absorption(sabs, temperature=tempdata[0])
        openmclib[name].set_scatter_matrix(scatter, temperature=tempdata[0])
        openmclib[name].set_fission(sf, temperature=tempdata[0])
        openmclib[name].set_nu_fission(nusf, temperature=tempdata[0])
        openmclib[name].set_chi(chi, temperature=tempdata[0])
    return openmclib
#
def prepare_temperature_depended_mg(libgroup, conc, namenuclide, temperature):
    """
    Prepare multigroup data for calculation based with temperature dependent
    constant
    Paramertres:
    -----------
    libgroup : dictionary
    {temperature : dictonary { name of nuclide : element groupsection class }};
    :param conc: dict
    - dictionary with name of material : np.array - R*8 concentration of nuclide;
    :param namenuclide:
    - a list of nuclide names;
    :param temperature:
    - dictionary with name of material : temperature value;
    :return:
    - dictionary with name of material : openmc.MGXS class element
    """
    # TEMPERATURE DEPENDET CASE
    openmclib = {}
    for name, val in conc.items():
        openmclib[name] = openmc.XSdata(name, groups, temperatures=tempdata)
        openmclib[name].order = 1
        for tt in tempdata:
            stot = np.zeros(len(ENERGIES) - 1, dtype=np.double)
            sabs = np.zeros(len(ENERGIES) - 1, dtype=np.double)
            scapt = np.zeros(len(ENERGIES) - 1, dtype=np.double)
            sf = np.zeros(len(ENERGIES) - 1, dtype=np.double)
            nusf = np.zeros(len(ENERGIES) - 1, dtype=np.double)
            chi = np.zeros(len(ENERGIES) - 1, dtype=np.double)
            scatter = np.zeros((len(ENERGIES) - 1, len(ENERGIES) - 1, 2), dtype=np.double)
            concentration = 0.0
            for n, v in zip(namenuclide, val):
                stot += libgroup[tt][n].stot * v
                sabs += libgroup[tt][n].sabs * v
                scapt += libgroup[tt][n].scapt * v
                sf += libgroup[tt][n].sf * v
                nusf += libgroup[tt][n].nusf * v
                scatter += libgroup[tt][n].xsmatrix * v
                if (libgroup[tt][n].sf.sum() > 0):
                    concentration += v
                    chi += libgroup[tt][n].schi * v
            if (concentration > 0):
                chi = chi / concentration
            openmclib[name].set_total(stot, temperature=tt)
            openmclib[name].set_absorption(sabs, temperature=tt)
            openmclib[name].set_scatter_matrix(scatter, temperature=tt)
            openmclib[name].set_fission(sf, temperature=tt)
            openmclib[name].set_nu_fission(nusf, temperature=tt)
            openmclib[name].set_chi(chi, temperature=tt)
    return openmclib


res,resl, dicel, el = make_model()
openmclib = {}
libgroup = temp_libgroup(tempdata)
for key,value in res.items():
    conc = {}
    temp = {}
    key_1 = ''+str(key)+''
    res[key],resl=collapse_into_nuclide(res[key],resl, dicel[key], el)
    conc[key_1] = res[key]
    temp[key_1] = 600.0  # for all zones a temperature the same : 600.0
    openmclib[key] = prepare_temperature_independed_mg(libgroup, conc, resl, tempdata, groups, temp)

######
# key = (546,7)
# key_1 = ''+str(key)+''
# res[key],resl=collapse_into_nuclide(res[key],resl, dicel[key], el) # from Table_Mend.py
# conc = {}; conc[key_1] = res[key]
# temp = {}; temp[key_1] = 600.0 # for all zones a temperature the same : 600.0
# openmclib[key] = prepare_temperature_independed_mg(libgroup, conc, resl, tempdata, groups, temp)
######

mg_cross_sections_file = openmc.MGXSLibrary(groups)
print("next step")
mg_cross_sections_file.add_xsdatas([openmclib[o] for o in openmclib])
print("last step")
mg_cross_sections_file.export_to_hdf5()
