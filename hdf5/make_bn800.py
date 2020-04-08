import openmc
import numpy as np
from math import *
from create import *
###############################################################################
#                      Function
###############################################################################

##
NASS   = 817
NRING  = 17
HPITCH = 10.04
#
#
PATH_REACTION_XML = '/home/barakuda/Рабочий стол/hdf5_openmc/hdf5/reactions.xml'
PATH_INOUT_DATA = '/home/barakuda/yur/archive/800/InOut6800.dat.hdf5'
def make_csnucls():
    import lxml.etree as ET
    root = ET.parse(str(os.environ['OPENMC_CROSS_SECTIONS']))
    elements = []
    for p in root.iter():
        if 'materials' in p.keys():
           elements.append(p.get('materials'))
    return elements

def get_nuclids_from_hdf5(csname):
    import h5py
    f = h5py.File(csname)
    return [k for k in f]

def densna(TS):
    TS1 = TS + 273.15
    SRO = (0.89660679 + 0.5161343E-03*TS1 - 1.8297218E-06*TS1*TS1 +
           2.2016247E-09*TS1**3 - 1.3975634E-12*TS1**4 +
           0.44866894E-15*TS1**5 - 0.057963628E-18*TS1**6)
    return SRO
#
def translate_actinide(inp_dict, names, mode_nuclide):
    out_dict = {}
    _nuclide_elemements = {}
    sections = make_csnucls()
    if (not mode_nuclide):
        for n in names:
            if (not n[-2:].isdigit()):
                _nuclide_elemements[n] = openmc.Element(n).expand(1.0, 'ao')
    for k, v in inp_dict.items():
        print("proceed material {}".format(k))
        mat = openmc.Material(material_id = 100*k[0] + k[1] + 1,name=str(k))
        mat.set_density('sum')
        for name, value in zip(names, v):
            if (value > 0.0):
                if (mode_nuclide):
                    #mat.add_nuclide(name, value)
                    if (name in sections): mat._nuclides.append((name, value, 'ao'))
                else:
                    if (name[-2:].isdigit()):
                        mat._nuclides.append((name, value, 'ao'))
                    else:
                        for vv in _nuclide_elemements[name]:
                            mat._nuclides.append((vv[0], vv[1] * value, 'ao'))
                #if (name[-2:].isdigit()):
                #    mat.add_nuclide(name, value)
                #else:
                #    mat.add_element(name, value)
        out_dict[k] = mat
    return out_dict

def translate_cm(inp_dict, names, out_dict, mode_nuclide):
    _nuclide_elemements = {}
    sections = make_csnucls()
    if (not mode_nuclide):
        for n in names:
            if (not n[-2:].isdigit()):
                _nuclide_elemements[n] = openmc.Element(n).expand(1.0, 'ao')
    for k, v in inp_dict.items():
        print("proceed material in element {}".format(k))
        mat = out_dict[k]
        for name, value in zip(names, v):
            if (value > 0):
                if (mode_nuclide):
                    if (name in sections): mat.add_nuclide(name, value)
                else:
                    if (name[-2:].isdigit()):
                        mat._nuclides.append((name, value, 'ao'))
                    else:
                        for vv in _nuclide_elemements[name]:
                            mat._nuclides.append((vv[0], vv[1] * value, 'ao'))
                #if (name[-2:].isdigit()):
                #    mat.add_nuclide(name, value)
                #else:
                #    mat.add_element(name, value)
    return out_dict

def translate_mg(nkas, rods, dzz, hss, namelib=None, namen=None, ron = None,
                    namer=None, ror = None):
    """Create multigroup library
    Parameters:
    ----------
    nkas : int
       number of assemblies in a library
    rods : iterable of boolean
       whether assembly or rod
    dzz : iterable of float
       list with asseblies axial discretiziation
    hss : iterable of float
       list with asseblies axial discretiziation
    --lower Parametres for by nuclide mg mode--
    namelib : str
       name of mgxs library. By default None
    namen : list of str
       name of nuclide in mgxs library for assembly. By default None
    ron : dict
       (Nass, Nheigtlayher) : np.array of concetration. By default None
    namer : list of str for rods
       name of nuclide in mgxs library for rods. By default None
    ror : dict for rods
       (Nass, Nheigtlayher) : np.array of concetration. By default None
    Returns:
    -------
    out_dict : dict
       dictionary with key (Nass, Nheight_number) and value openmc.Material
    """
    out_dict = {}
    if namelib is not None:
        validnuclides = get_nuclids_from_hdf5(namelib)
    for i in range(nkas):
        print("current assembly is {}".format(i))
        if (rods[i]):
            for j, v in enumerate(hss):
                string = '({}, {})'.format(i + 1, j + 1)
                mat = openmc.Material(material_id = 100*(i+1) + j + 1,
                                         name = string)
                mat.set_density('sum')
                if namelib is None:
                    mat.add_macroscopic(string)
                else:
                    for name, value in zip(namer, ror[string]):
                        if name in validnuclides:
                            mat._nuclides.append((name, value, 'ao'))
                out_dict[(i + 1, j + 1)] = mat
        else:
            for j, v in enumerate(dzz):
                string = '({}, {})'.format(i + 1, j + 1)
                mat = openmc.Material(material_id = 100*(i+1) + j + 1,
                                         name = string)
                mat.set_density('sum')
                if namelib is None:
                    mat.add_macroscopic(string)
                else:
                    for name, value in zip(namen, ron[string]):
                        if name in validnuclides:
                            mat._nuclides.append((name, value, 'ao'))
                out_dict[(i + 1, j + 1)] = mat
    return out_dict

def init_openmc_spiral(openmc_object,nring = 2):
    tot_arr = [] # Resulting array of rings
    for ring in range(nring - 1,-1,-1):
        arr = []
        arr.append(openmc_object)
        for cell in range(ring * 6 - 1):
            arr.append(openmc_object)
        tot_arr.append(arr)
    return tot_arr
#
def get_hex_num_by_spiral(spiral_num,nring = 2):
    if(spiral_num > 3*(nring - 1)*nring + 1): raise Exception
    for ring in range(nring,0,-1):
        if (spiral_num == 1):
            return [nring-1,0]
        else:
            if (spiral_num > 3*(ring - 2)*(ring - 1) + 1):
                c1 = nring - ring
                c2 = 3*(ring - 1)*ring + 1 - spiral_num
                if (c2 == 6*(ring - 1) - 1):
                    return [c1,0]
                else:
                    return [c1,c2+1]
#
###3D
def make_load(name, outer, assemblies, number, hpitch, lenz, load_dict = {}):
    """Parametres:
       ----------
       name : str
           - load called
       outer : openmc.Universe
           - outer universe
       assemblies : dictionary
           - loaded type of assembly
       load_dict : dictionary
           - a pairs 'name of lattice3d' : a 'list of spiral numbers'
       Return:
       ----------
       openmc.Lattice
    """
    lattice = openmc.HexLattice(lattice_id=1001, name=name)
    lattice.center = (0, 0, 0)
    lattice.orientation = "x"
    lattice.pitch = (hpitch, lenz)
    lattice.outer = outer
    univs = []
    univs.append(init_openmc_spiral(outer, number))
    lattice.universes = univs
    for k,v in load_dict.items():
        #
        for g in v:
            zz = get_hex_num_by_spiral(g, number)
            i = zz[0]
            j = zz[1]
            lattice.universes[0][i][j] = assemblies[k]
    return lattice
#
zcore = np.ones(25)
rdzcore = [zcore.sum()]
#
def make_surfaces(zcore, hpitch):
    overall_size = sum(zcore)
    Z = [zcore[0]]
    for i, zz in enumerate(zcore[1:]):
        Z.append(zz + Z[i])
    surfaces = []
    surfaces.append(openmc.ZPlane(z0=0.0 - overall_size/2, boundary_type = 'vacuum'))
    for i,z in enumerate(Z[:-1]):
        surfaces.append(openmc.ZPlane(z0=z - overall_size/2))
    surfaces.append(openmc.ZPlane(z0=Z[-1] - overall_size/2, boundary_type = 'vacuum'))
    edge_length = (1./np.sqrt(3.0)) * hpitch
    hexsurface = openmc.model.get_hexagonal_prism(edge_length=edge_length,
                                                  origin=(0.0, 0.0))
    return hexsurface, surfaces

def _make_outer_universe():
    mat  = openmc.Material(material_id=100000)
    mat.set_density('g/cm3',0.88)
    mat.add_nuclide("Na23", 1.0)
    cell = openmc.Cell(cell_id=100000)
    univ = openmc.Universe(universe_id=100000)
    cell.fill = mat
    univ.add_cell(cell)
    outsurface = openmc.ZCylinder(R=175.25, name='Reactor preasure vessel out', boundary_type = 'vacuum')
    return univ, outsurface

def _make_outer_mg_universe(outmat):
    cell = openmc.Cell(cell_id=100000)
    univ = openmc.Universe(universe_id=100000)
    cell.fill = outmat
    univ.add_cell(cell)
    outsurface = openmc.ZCylinder(R=175.25, name='Reactor preasure vessel out', boundary_type = 'vacuum')
    return univ, outsurface

def make_abl_cells(materials, hexsurface, surfaces, iass=0, densna=[]):
    cells = {}
    univ = openmc.Universe(universe_id = iass)
    for j, s in enumerate(surfaces[:-1]):
        cells[(iass, j)] = openmc.Cell(cell_id = iass*100 + j + 1)
        cells[(iass, j)].region = hexsurface & +surfaces[j] & -surfaces[j+1]
        cells[(iass, j)].temperature = 300.0
        materials[(iass, j + 1)].add_nuclide("Na23", densna[iass + 1][j]*0.6022/23.0)
        cells[(iass, j)].fill = materials[(iass, j + 1)]
        univ.add_cell(cells[(iass, j)])
    return univ

def make_abl_mg_cells(materials, hexsurface, surfaces, iass=0, temperature=None):
    cells = {}
    univ = openmc.Universe(universe_id = iass)
    for j, s in enumerate(surfaces[:-1]):
        cells[(iass, j)] = openmc.Cell(cell_id = iass*100 + j + 1)
        cells[(iass, j)].region = hexsurface & +surfaces[j] & -surfaces[j+1]
        if temperature is None:
            cells[(iass, j)].temperature = 300.0
        else:
            cells[(iass, j)].temperature = temperature[(iass, j + 1)]
        cells[(iass, j)].fill = materials[(iass, j + 1)]
        univ.add_cell(cells[(iass, j)])
    return univ


def make_model_continuous():
    element = name_el(PATH_INOUT_DATA)
    nuclide = xml_f(PATH_REACTION_XML)
    print("Start")
    print(nuclide, element)
    diction_1, diction_2, diction_3, dzz, hss, rods, densna, type_fuel, type_assembly = run(PATH_INOUT_DATA)
    print('First point', dzz, hss,rods)
    outuniv, outsurface = _make_outer_universe()
    diction_1_out = translate_actinide(diction_1, nuclide, True)
    diction_ass   = translate_cm(diction_2, element, diction_1_out, False)
    diction_rod   = translate_actinide(diction_3, element, False)
    #
    assemblies = {}
    load_dict  = {}
    hexsurf, rodfaces = make_surfaces([sum(dzz)], HPITCH)
    hexsurf, zfaces   = make_surfaces(dzz, HPITCH)
    for i in range(NASS):
        print("PROCEED {} assembly".format(i+1))
        if (rods[i]):
            assemblies[i+1] = make_abl_cells(diction_rod, hexsurf, rodfaces, i+1, densna)
            load_dict[i+1]  = [i+1]
        else:
            assemblies[i+1] = make_abl_cells(diction_ass, hexsurf, zfaces, i+1, densna)
            load_dict[i+1]  = [i+1]
    print("Making load ..")

    ##################################################################################
    load = make_load("BN-800", outuniv, assemblies, NRING, HPITCH, sum(dzz), load_dict)
    root_cell = openmc.Cell(cell_id=0, region=-outsurface, fill=load)
    root_univ = openmc.Universe(universe_id=0, cells = [root_cell])
    # Instantiate a Geometry, register the root Universe, and export to XML
    openmc_geometry = openmc.Geometry()
    openmc_geometry.root_universe = root_univ
    #openmc_geometry.export_to_xml()
    materials_list = openmc_geometry.get_all_materials()
    materials_file = openmc.Materials([v for v in materials_list.values()])
    #materials_file.cross_sections = 'mgxs.h5'
    #materials_file.export_to_xml()
    # Instantiate a Settings object, set all runtime parameters, and export to XML
    settings_file = openmc.Settings()
    ##########################SETTINGS
    settings_file.batches = 120
    settings_file.inactive = 20
    settings_file.particles = 1000
    ###########################
    settings_file.temperature = {'range':(250,2500)}
    # Tell OpenMC this is a multi-group problem
    #settings_file.energy_mode = 'multi-group'
    #settings_file.survival_biasing = True
    #settings_file.cutoff = {'energy_neutron' : 1.0,'weight':1.e-10}
    settings_file.source = openmc.Source(space=openmc.stats.Box(
        [-160.04, -160.04, 0.0], [160.04, 160.04, 170], only_fissionable=True))
    settings_file.export_to_xml()
    plot_1 = openmc.Plot(plot_id=1)
    plot_1.filename = 'plot_1'
    plot_1.origin = [0.0, 0.0, 8.]
    plot_1.width = [520.26, 520.26]
    plot_1.pixels = [2400, 2400]
    plot_1.color = 'mat'
    plot_1.basis = 'xy'

    plot_2 = openmc.Plot(plot_id=2)
    plot_2.filename = 'plot_2'
    plot_2.origin = [0.0, 0.0, 0.0]
    plot_2.width = [520.26, 350.2]
    plot_2.pixels = [1200, 1200]
    plot_2.color = 'mat'
    plot_2.basis = 'xz'

    # Instantiate a Plots collection and export to XML
    plot_file = openmc.Plots([plot_1, plot_2])
    #plot_file.export_to_xml()

def make_datas():
    res,resl, dicel, el, rodel, dolna, dzz, hss, rods, type_fuel, type_assembly = make_model()
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
        temp[key] = 600.0  # for all zones a temperature the same : 600.0
        concentration = 0.0*concentration[:]
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
        temprod[key] = 600.0  # for all zones a temperature the same : 600.0
        concentration = 0.0*concentration[:]
    print("DATAS MADE!!!")
    return nuclist, conc, temp, rodnuclist, concrod, temprod, dzz, hss, rods, type_fuel, type_assembly

def make_model_mg(nameoflib):
    element = name_el(PATH_INOUT_DATA)
    nuclide = xml_f(PATH_REACTION_XML)
    print("Start")
    print(nuclide, element)
    diction_1, diction_2, diction_3, dzz, hss, rods, densna, type_fuel, type_assembly = run(PATH_INOUT_DATA)
    print('First point', dzz, hss,rods)
    nkas = len(rods)
    diction_mg   = translate_mg(nkas, rods, dzz, [sum(dzz)])
    outuniv, outsurface = _make_outer_mg_universe(diction_mg[(1,1)])
    #
    assemblies = {}
    load_dict  = {}
    hexsurf, rodfaces = make_surfaces([sum(dzz)], HPITCH)
    hexsurf, zfaces   = make_surfaces(dzz, HPITCH)
    for i in range(NASS):
        if (rods[i]):
            assemblies[i+1] = make_abl_mg_cells(diction_mg, hexsurf, rodfaces, i+1)
            load_dict[i+1]  = [i+1]
        else:
            assemblies[i+1] = make_abl_mg_cells(diction_mg, hexsurf, zfaces, i+1)
            load_dict[i+1]  = [i+1]
    print("Making load ..")
    ##################################################################################
    load = make_load("BN-800", outuniv, assemblies, NRING, HPITCH, sum(dzz), load_dict)
    root_cell = openmc.Cell(cell_id=0, region=-outsurface, fill=load)
    root_univ = openmc.Universe(universe_id=0, cells = [root_cell])
    # Instantiate a Geometry, register the root Universe, and export to XML
    openmc_geometry = openmc.Geometry()
    openmc_geometry.root_universe = root_univ
    openmc_geometry.export_to_xml()
    materials_list = openmc_geometry.get_all_materials()
    materials_file = openmc.Materials([v for v in materials_list.values()])
    materials_file.cross_sections = nameoflib
    materials_file.export_to_xml()
    # Instantiate a Settings object, set all runtime parameters, and export to XML
    settings_file = openmc.Settings()
    ##########################SETTINGS
    settings_file.batches = 120
    settings_file.inactive = 20
    settings_file.particles = 1000
    ###########################
    settings_file.temperature = {'range':(250,2500)}
    # Tell OpenMC this is a multi-group problem
    settings_file.energy_mode = 'multi-group'
    settings_file.source = openmc.Source(space=openmc.stats.Box(
        [-160.04, -160.04, 0.0], [160.04, 160.04, 170], only_fissionable=True))
    settings_file.export_to_xml()
    return openmc
    plot_1 = openmc.Plot(plot_id=1)
    plot_1.filename = 'plot_1'
    plot_1.origin = [0.0, 0.0, 8.]
    plot_1.width = [520.26, 520.26]
    plot_1.pixels = [2400, 2400]
    plot_1.color = 'mat'
    plot_1.basis = 'xy'

    plot_2 = openmc.Plot(plot_id=2)
    plot_2.filename = 'plot_2'
    plot_2.origin = [0.0, 0.0, 0.0]
    plot_2.width = [520.26, 350.2]
    plot_2.pixels = [1200, 1200]
    plot_2.color = 'mat'
    plot_2.basis = 'xz'

    # Instantiate a Plots collection and export to XML
    plot_file = openmc.Plots([plot_1, plot_2])
    #plot_file.export_to_xml()


def make_model_mg_by_nuclide(nameoflib):

    nuclist, conc, temp, rodnuclist, concrod, temprod, dzz, hss, rods, type_fuel, type_assembly = make_datas()
    nkas = len(rods)
    diction_mg = translate_mg(nkas, rods, dzz, [sum(dzz)], nameoflib,
                                 nuclist, conc, rodnuclist, concrod)
    outuniv, outsurface = _make_outer_mg_universe(diction_mg[(1,1)])
    assemblies = {}
    load_dict  = {}
    hexsurf, rodfaces = make_surfaces([sum(dzz)], HPITCH)
    hexsurf, zfaces   = make_surfaces(dzz, HPITCH)
    for i in range(NASS):
        if (rods[i]):
            assemblies[i+1] = make_abl_mg_cells(diction_mg, hexsurf, rodfaces, i+1, temprod)
            load_dict[i+1]  = [i+1]
        else:
            assemblies[i+1] = make_abl_mg_cells(diction_mg, hexsurf, zfaces, i+1, temp)
            load_dict[i+1]  = [i+1]
    print("Making load ..")
    ##################################################################################
    load = make_load("BN-800", outuniv, assemblies, NRING, HPITCH, sum(dzz), load_dict)
    root_cell = openmc.Cell(cell_id=0, region=-outsurface, fill=load)
    root_univ = openmc.Universe(universe_id=0, cells = [root_cell])
    # Instantiate a Geometry, register the root Universe, and export to XML
    openmc_geometry = openmc.Geometry()
    openmc_geometry.root_universe = root_univ
    openmc_geometry.export_to_xml()
    materials_list = openmc_geometry.get_all_materials()
    materials_file = openmc.Materials([v for v in materials_list.values()])
    materials_file.cross_sections = nameoflib
    materials_file.export_to_xml()
    # Instantiate a Settings object, set all runtime parameters, and export to XML
    settings_file = openmc.Settings()
    ##########################SETTINGS
    settings_file.batches = 120
    settings_file.inactive = 20
    settings_file.particles = 1000
    ###########################
    settings_file.temperature = {'range':(250,2500)}
    # Tell OpenMC this is a multi-group problem
    settings_file.energy_mode = 'multi-group'
    #settings_file.survival_biasing = True
    #settings_file.cutoff = {'energy_neutron' : 1.0,'weight':1.e-10}
    settings_file.source = openmc.Source(space=openmc.stats.Box(
        [-160.04, -160.04, 0.0], [160.04, 160.04, 170], only_fissionable=True))
    settings_file.export_to_xml()
    plot_1 = openmc.Plot(plot_id=1)
    plot_1.filename = 'plot_1'
    plot_1.origin = [0.0, 0.0, 8.]
    plot_1.width = [520.26, 520.26]
    plot_1.pixels = [2400, 2400]
    plot_1.color = 'mat'
    plot_1.basis = 'xy'

    plot_2 = openmc.Plot(plot_id=2)
    plot_2.filename = 'plot_2'
    plot_2.origin = [0.0, 0.0, 0.0]
    plot_2.width = [520.26, 350.2]
    plot_2.pixels = [1200, 1200]
    plot_2.color = 'mat'
    plot_2.basis = 'xz'

def make_model():
    element = name_el(PATH_INOUT_DATA)
    nuclide = xml_f(PATH_REACTION_XML)
    print("Start")
    print(nuclide, element)
    diction_1, diction_2, diction_3, dzz, hss, rods, densna, type_fuel, type_assembly = run(PATH_INOUT_DATA)
    return diction_1, nuclide, diction_2, element, diction_3, densna

def draw_plot():
    print("draw")