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
# Instantiate the energy group data
groups = openmc.mgxs.EnergyGroups(np.array(ENERGIES))



def temp_interolate(temparray, T):
    if (T <= temparray[0]):
        return [(temparray[0], 1)]
    elif (T >= temparray[-1]):
        return [(temparray[-1], 1), (temparray[-1], 0)]
    else:
        for j,_t in enumerate(temparray):
            if (_t > T):
                i = j-1
                t = temparray[i]
                break
        return [(i, (temparray[i+1] -T)/(temparray[i+1] - temparray[i])), (i + 1, (T - t)/(temparray[i+1] - temparray[i]))]
        #return [(t, (temparray[i+1] -T)/(temparray[i+1] - temparray[i])), (temparray[i+1], (T - t)/(temparray[i+1] - temparray[i]))]

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

openmclib = {}

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
def prepare_mg(libgroup, conc, namenuclide, tempdata, groups, mattemp):
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
    :return:
    - dictionary with name of material : openmc.MGXS class element
    """
    # TEMPERATURE INDEPENDET CASE
    openmclib = {}
    t0 = time.time()
    nsize = len([k for k in conc])
    values = np.array([v for v in conc.values()])
    values = values.reshape(nsize, len(values[0]))
    indices = nsize * [(tempdata[0], 1.0)]
    indarray=np.zeros((nsize, 2), dtype = np.int)
    wgtarray=np.zeros((nsize, 2), dtype = np.double)
    stot = np.zeros((nsize, len(ENERGIES) - 1), dtype=np.double)
    sabs = np.zeros((nsize,len(ENERGIES) - 1), dtype=np.double)
    scapt = np.zeros((nsize,len(ENERGIES) - 1), dtype=np.double)
    sf = np.zeros((nsize,len(ENERGIES) - 1), dtype=np.double)
    nusf = np.zeros((nsize,len(ENERGIES) - 1), dtype=np.double)
    chi = np.zeros((nsize,len(ENERGIES) - 1), dtype=np.double)
    scatter = np.zeros((nsize,len(ENERGIES) - 1, len(ENERGIES) - 1, 2),
                        dtype=np.double)
    for i, name in enumerate(mattemp.keys()):
        openmclib[name] = openmc.XSdata(name, groups, temperatures=[tempdata[0]])
        openmclib[name].order = 1
        indices[i] = temp_interolate(tempdata, mattemp[name])
        indarray[i, 0]= indices[i][0][0];indarray[i, 1]= indices[i][1][0]
        wgtarray[i, 0]= indices[i][0][1];wgtarray[i, 1]= indices[i][1][1]
    nuclind = np.zeros((len(namenuclide), len(tempdata)), dtype=np.int)
    for i, n in enumerate(namenuclide):
        for j, tt in enumerate(tempdata):
            if (n in libgroup[tt].keys()):
                nuclind[i][j] = i + 1
    t1 = time.time()
    for i in range(nsize):
        for ind in range(len(namenuclide)):
            for j in [0, 1]:
                if (nuclind[ind][indarray[i, j]] > 0):
                    stot[i, :] += libgroup[tempdata[indarray[i, j]]][namenuclide[ind]].stot * values[i, ind] * wgtarray[i, j]
                    sabs[i, :] += libgroup[tempdata[indarray[i, j]]][namenuclide[ind]].sabs * values[i, ind] * wgtarray[i, j]
                    scapt[i, :] += libgroup[tempdata[indarray[i, j]]][namenuclide[ind]].scapt * values[i, ind] * wgtarray[i, j]
                    sf[i, :] += libgroup[tempdata[indarray[i, j]]][namenuclide[ind]].sf * values[i, ind] * wgtarray[i, j]
                    nusf[i, :] += libgroup[tempdata[indarray[i, j]]][namenuclide[ind]].nusf * values[i, ind] * wgtarray[i, j]
                    scatter[i, :, :] += libgroup[tempdata[indarray[i, j]]][namenuclide[ind]].xsmatrix * values[i, ind] * wgtarray[i, j]
            concentration = 0.0
            if (namenuclide[ind] in libgroup[tempdata[0]].keys()):
                if (libgroup[tempdata[0]][namenuclide[ind]].sf.sum() > 0):
                    concentration += values[i, ind]
                    chi[i, :] += libgroup[tempdata[0]][namenuclide[ind]].schi * values[i, ind]
    for i, name in enumerate(conc.keys()):
        openmclib[name]._total[0]=stot[i, :]
        openmclib[name]._absorption[0]=sabs[i, :]
        openmclib[name]._scatter_matrix[0]=scatter[i, :, :]
        openmclib[name]._fission[0]=sf[i, :]
        if (sum(sf[i, :]) > 0):
            openmclib[name]._fissionable = True
        openmclib[name]._nu_fission[0]=nusf[i, :]
        openmclib[name]._chi[0]=chi[i, :]
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



def make_libgroup_micro(nameoflib):
    libgroup = temp_libgroup(tempdata)
    openmclib = {}
    for tt in tempdata:
        for n in libgroup[tt]:
            scatter = np.zeros((len(ENERGIES) - 1, len(ENERGIES) - 1, 2), dtype=np.double)
            if (n in openmclib):
                openmclib[n].set_total(libgroup[tt][n].stot, temperature=tt)
                openmclib[n].set_absorption(libgroup[tt][n].sabs, temperature=tt)
                scatter += libgroup[tt][n].xsmatrix*1.0
                openmclib[n].set_scatter_matrix(scatter, temperature=tt)
                openmclib[n].set_fission(libgroup[tt][n].sf, temperature=tt)
                openmclib[n].set_nu_fission(libgroup[tt][n].nusf, temperature=tt)
                openmclib[n].set_chi(libgroup[tt][n].schi, temperature=tt)
            else:
                openmclib[n] = openmc.XSdata(n, groups, temperatures=tempdata)
                if (libgroup[tt][n].xsmatrix.shape[-1] < 2):
                    openmclib[n].order = 0
                else:
                    openmclib[n].order = 1
                openmclib[n].order = 1
                openmclib[n].set_total(libgroup[tt][n].stot, temperature=tt)
                openmclib[n].set_absorption(libgroup[tt][n].sabs, temperature=tt)
                scatter += libgroup[tt][n].xsmatrix*1.0
                openmclib[n].set_scatter_matrix(scatter, temperature=tt)
                openmclib[n].set_fission(libgroup[tt][n].sf, temperature=tt)
                openmclib[n].set_nu_fission(libgroup[tt][n].nusf, temperature=tt)
                openmclib[n].set_chi(libgroup[tt][n].schi, temperature=tt)
    mg_cross_sections_file = openmc.MGXSLibrary(groups)
    mg_cross_sections_file.add_xsdatas([openmclib[o] for o in openmclib])
    mg_cross_sections_file.export_to_hdf5(nameoflib)



###start
from make_bn800 import *
def make_libgroup_macro(nameoflib):
   res,resl, dicel, el, rodel, dolna = make_model()
   libgroup = temp_libgroup(tempdata)
   nuclval = len(resl)
   indexNa = el.index("Na")
   for key, value in dicel.items():
       value[indexNa] = densna(600.0) * 0.6022 / 23 * dolna[key[0]][key[1]]
   for key, value in rodel.items():
       value[indexNa] = densna(600.0) * 0.6022 / 23 * dolna[key[0]][key[1]]
   nuclist, element_from, element_ind, element_val = get_unite_list(resl, el)
   concentration = np.zeros(len(nuclist))
   conc = {}
   temp = {}
   for key,value in res.items():
       t0 = time.time()
       concentration[:nuclval] = value
       for f, i, v in zip(element_from, element_ind, element_val):
           concentration[i] += dicel[key][f] * v
       key_1 = ''+str(key)+''
       conc[key_1] = concentration
       temp[key_1] = 600.0  # for all zones a temperature the same : 600.0
       concentration = 0.0*concentration[:]
       print("Estimated time is {} min".format((time.time() - t0)))
   rodnuclist, element_from, element_ind, element_val = get_unite_list([], el)
   concentration = np.zeros(len(rodnuclist))
   concrod = {}
   temprod = {}
   for key,value in rodel.items():
       t0 = time.time()
       for f, i, v in zip(element_from, element_ind, element_val):
           concentration[i] += rodel[key][f] * v
       key_1 = ''+str(key)+''
       concrod[key_1] = concentration
       temprod[key_1] = 600.0  # for all zones a temperature the same : 600.0
       concentration = 0.0*concentration[:]
       print("Estimated time is {} min".format((time.time() - t0)))
   openmclib = prepare_mg(libgroup, conc, resl, tempdata,
                             groups, temp)
   openmclib.update(prepare_mg(libgroup, concrod, rodnuclist, tempdata,
                                 groups, temprod))
   mg_cross_sections_file = openmc.MGXSLibrary(groups)
   mg_cross_sections_file.add_xsdatas([openmclib[o] for o in openmclib])
   mg_cross_sections_file.export_to_hdf5(nameoflib)
