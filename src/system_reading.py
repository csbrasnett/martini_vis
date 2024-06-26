
from vermouth.gmx import read_itp
from vermouth.forcefield import ForceField
from .topology import input_topol_reader
def _misc_file_reader(lines):
    """
    Try to handle miscellaneous files which can't be read into the force field by read_itp

    Parameters
    ----------
    lines: list
        list of lines from the input .itp file


    Returns
    -------
    defines: dict
        dictionary of definition: parameters for bonded terms that have been externally defined
    others: list
        lines from an input itp that don't have #define in them

    if others == lines, None is returned, because the input file doesn't contain anything interesting
    """
    defines = {i.split()[1]: i.split(';')[0].split()[2:5] for i in lines if '#define' in i}
    others = [line for line in lines if '#define' not in line]
    if others == lines:
        return None
    else:
        return defines, others

def system_reading(topology):

    """
    read a .top file's contents into a ForceField
    """

    # get the topology file
    print(f"Reading input {topology}")
    topol_lines = input_topol_reader(topology)

    # for each molecule in the system, read in the itp
    d = {}
    for i, j in enumerate(topol_lines['core_itps']):
        with open(j) as f:
            d[i] = f.readlines()

    # read the molecules into the forcefield
    ff = ForceField('martini3001')

    system_defines = {}
    for i, j in enumerate(d.keys()):
        try:
            read_itp(d[j], ff)
        except OSError:
            '''
            if we can't read the file into the system directly, we have something that isn't strictly a molecule
            most likely its the force field definition file (eg. martini_v3.0.0.itp) but we can't be sure
            a common one is something with #defines in for generic molecule bonded terms
            '''
            misc_result = _misc_file_reader(d[j])
            if misc_result is None:
                # if [ defaults ] is found, then it's martini_v3.0.0.itp or similar. ignore it.
                if "[ defaults ]" not in [k.strip() for k in d[j]]:
                    print(f"Error reading {topol_lines['core_itps'][i]}. Will ignore and exclude from output system.")
            else:
                # this means we've separated things out successfully and can read the actual itp content now
                # tracking system_defines is useful for when we sort the #TODOs below.
                for key, value in misc_result[0].items():
                    system_defines[key] = value
                read_itp(misc_result[1], ff)

    return ff, topol_lines, system_defines
