from create import *


def translate_function(inp_dict, names, mode_nuclide):
    out_dict = {}
    for k, v in inp_dict.items():
        mat = openmc.Material(name=str(k))
        mat.set_density('sum')
        for name, value in zip(names, v):
            if (mode_nuclide):
                mat.add_nuclide(name, value)
            else:
                try:
                    mat.add_element(name, value)
                except ValueError:
                    mat.add_nuclide(name, value)
        out_dict[k] = mat
    return out_dict


if __name__ == "__main__":
    start_time = time.time()
    n_arg = argv[1]
    m_arg = argv[2]

    nuclide = xml_f(m_arg)
    element = name_el(n_arg)

    dictions = run(n_arg)

    check(nuclide, "check_nuclide.txt")
    diction_1_out = translate_function(dictions[0], nuclide, True)
    diction_2_out = translate_function(dictions[1], element, False)
    diction_3_out = translate_function(dictions[2], element, False)
    print("%s seconds" % (time.time() - start_time))

# use this for check
# check(diction_2_out, "chec.txt")
# check(diction_1_out, "chec_1.txt")
# check(diction_3_out,"chec_3.txt")
