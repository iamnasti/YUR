def determ_el(el):
    if el[0:2] == '11':
        return 'Na'
    elif el[0:2] == '12':
        return 'Mg'
    elif el[0:2] == '13':
        return 'Al'
    elif el[0:2] == '14':
        return 'Si'
    elif el[0:2] == '20':
        return 'Ca'
    elif el[0:2] == '22':
        return 'Ti'
    elif el[0:2] == '23':
        return 'V'
    elif el[0:2] == '24':
        return 'Cr'
    elif el[0:2] == '25':
        return 'Mn'
    elif el[0:2] == '26':
        return 'Fe'
    elif el[0:2] == '27':
        return 'Co'
    elif el[0:2] == '28':
        return 'Ni'
    elif el[0:2] == '29':
        return 'Cu'
    elif el[0:2] == '40':
        return 'Zr'
    elif el[0:2] == '41':
        return 'Nb'
    elif el[0:2] == '42':
        return 'Mo'
    elif el[0:2] == '43':
        return 'Tc'
    elif el[0:2] == '51':
        return 'Sb'
    elif el[0:2] == '74':
        return 'W'
    elif el[0:1] == '1':
        return 'H'
    elif el[0:1] == '4':
        return 'Be'
    elif el[0:1] == '5':
        return 'B'
    elif el[0:1] == '6':
        return 'C'
    elif el[0:1] == '7':
        return 'N'
    elif el[0:1] == '8':
        return 'O'

def translaterlist(el):
      elnuclides = ['O', 'Na', 'Mo', 'Mn', 'Zr', 'Mg', 'Al', 'Ca', 'Ti', 'V', 'Co', 'Cu', 'H', 'Si', 'C', 'B10', 'B11']
      outelnuclides = ['O16', 'Na23', 'Mo0', 'Mn55',  'Zr0', 'Mg0', 'Al27', 'Ca0', 'Ti0', 'V0', 'Co59', 'Cu0', 'H1',  'Si0', 'C0', 'B10', 'B11']
      excl = [ 'Fe', 'Cr', 'Ni', 'Nb', 'N', 'Sb']
      if (el in elnuclides):
          return [(outelnuclides[elnuclides.index(el)], 1.0)]
      elif (el in excl):
          e = []
          ee = openmc.Element(el).expand(1, 'ao')
          for v in ee:
                e.append((v[0], v[1]))
          return e
      else:
          return []
#
def get_unite_list(nuclist, ellist):
    element_from = []
    element_ind = []
    element_val = []
    for j, e in enumerate(ellist):
        for ee in translaterlist(e):
            if (ee[0] in nuclist):
                element_from.append(j)
                element_ind.append(nuclist.index(ee[0]))
                element_val.append(ee[1])
            else:
                element_from.append(j)
                nuclist.append(ee[0])
                element_ind.append(len(nuclist) - 1)
                element_val.append(ee[1])
    return nuclist,element_from, element_ind, element_val
