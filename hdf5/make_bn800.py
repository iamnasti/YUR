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
def make_csnucls():
    import lxml.etree as ET
    root = ET.parse(str(os.environ['OPENMC_CROSS_SECTIONS']))
    elements = []
    for p in root.iter():
        if 'materials' in p.keys():
           elements.append(p.get('materials'))
    return elements
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

def make_abl_cells(materials, hexsurface, surfaces, iass=0, densna=[]):
    cells = {}
    univ = openmc.Universe(universe_id = iass)
    for j, s in enumerate(surfaces[:-1]):
        cells[(iass, j)] = openmc.Cell(cell_id = iass*100 + j + 1)
        cells[(iass, j)].region = hexsurface & +surfaces[j] & -surfaces[j+1]
        cells[(iass, j)].temperature = 300.0
        materials[(iass, j + 1)].add_nuclide("Na23", densna[j,iass]*0.6022/23.0)
        cells[(iass, j)].fill = materials[(iass, j + 1)]
        univ.add_cell(cells[(iass, j)])
    return univ

def make_model():
    start_time = time.time()
    m_arg = '/home/barakuda/Рабочий стол/hdf5_openmc/hdf5/reactions.xml'
    n_arg = '/home/barakuda/yur/archive/800/InOut6800.dat.hdf5'
    element = name_el(n_arg)
    nuclide = xml_f(m_arg)
    print("Start")
    print(nuclide, element)
    diction_1, diction_2, diction_3, dzz, hss, rods, densna, type_fuel, type_assembly = run(n_arg)
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
    openmc_geometry.export_to_xml()    
    materials_list = openmc_geometry.get_all_materials()
    materials_file = openmc.Materials([v for v in materials_list.values()])
    materials_file.export_to_xml()
    # Instantiate a Settings object, set all runtime parameters, and export to XML
    settings_file = openmc.Settings()
    ##########################SETTINGS
    settings_file.batches = 120
    settings_file.inactive = 20
    settings_file.particles = 1000
    ###########################
    settings_file.temperature = {'range':(250,2500)}
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
    plot_file.export_to_xml()
    print("%s seconds" % (time.time() - start_time))
      

