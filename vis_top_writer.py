#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: chris
"""

from vermouth.gmx import read_itp, write_molecule_itp
from vermouth.forcefield import ForceField
import os
import argparse
from argparse import RawTextHelpFormatter
import copy
from pathlib import Path
import networkx as nx


def misc_file_reader(lines):
    '''
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
    '''
    defines = {i.split()[1]: i.split(';')[0].split()[2:5] for i in lines if '#define' in i}
    others = [line for line in lines if '#define' not in line]
    if others == lines:
        return None
    else:
        return defines, others


def go_writer(ff, molname, go_bonds):
    '''
    write a go network only topology for a particular molecule

    Parameters
    ----------
    ff: vermouth forcefield
        vermouth force field containing the input syste
    molname: str
        name of the molecule to separate out
    go_bonds: list
        list of go  bonds to write out

    Returns
    -------
    None
    '''
    # remove all interactions from the molecule
    for interaction_type in list(ff.blocks[molname].interactions):
        del ff.blocks[molname].interactions[interaction_type]

    # add the elastic network bonds back in
    for bond in go_bonds:
        ff.blocks[molname].add_interaction('bonds', [bond[0], bond[1]], list(bond[2:]))

    # write the file out
    mol_out = ff.blocks[molname].to_molecule()
    mol_out.meta['moltype'] = molname + '_go'

    header = [f'Elastic network topology for {molname}', 'NOT FOR SIMULATIONS']

    with open(molname + '_go.itp', 'w') as outfile:
        write_molecule_itp(mol_out, outfile=outfile, header=header)


def en_writer(ff, molname, en_bonds):
    '''
    write an elastic network only topology for a particular molecule

    Parameters
    ----------
    ff: vermouth forcefield
        vermouth force field containing the input syste
    molname: str
        name of the molecule to separate out
    en_bonds: list
        list of elastic network bonds to write out

    Returns
    -------
    None
    '''
    # remove all interactions from the molecule
    for interaction_type in list(ff.blocks[molname].interactions):
        del ff.blocks[molname].interactions[interaction_type]

    # add the elastic network bonds back in
    edges = []
    for bond in en_bonds:
        edges.append(bond.atoms)
        ff.blocks[molname].add_interaction('bonds', bond.atoms, bond.parameters)

    #make a graph, look at the degrees of the nodes
    G = nx.Graph()
    G.add_nodes_from(ff.blocks[molname].nodes)
    G.add_edges_from(edges)
    degrees = [G.degree[node] for node in G.nodes]
    # handle the points where more EN bonds have been written than VMD can handle (12)
    if any([i>12 for i in [G.degree[node] for node in G.nodes]]):

        print(f"There are atoms in {molname} which have > 12 elastic network bonds."
              " Some will be removed and recorded for posterity")

        over_limit_ind = []
        over_limit_num = []
        for i, j in enumerate(degrees):
            if j > 12:
                over_limit_ind.append(i)
                over_limit_num.append(j-12)

        l0 = []
        l1 = []
        l2 = []
        for interaction in ff.blocks[molname].interactions['bonds']:
            cond0 = (interaction.atoms[0] in over_limit_ind)
            cond1 = (interaction.atoms[1] in over_limit_ind)
            if cond0 and not cond1:
                l0.append([interaction, interaction.atoms[0]])
            elif cond1 and not cond0:
                l1.append([interaction, interaction.atoms[1]])
            elif cond0 and cond1:
                l2.append([interaction, interaction.atoms[0], interaction.atoms[1]])
                # print(interaction)
        X = [i[0] for i in l0]+[i[0] for i in l1]+[i[0] for i in l2]+[i[0] for i in l2]
        Y = [i[1] for i in l0]+[i[1] for i in l1]+[i[1] for i in l2]+[i[2] for i in l2]
        Z = [list(x) for x in sorted(zip(Y,X), key=lambda pair: pair[0])]

        # make sure we have a consistent set of indices to remove edges from
        s0 = set([x[0] for x in Z])
        s1 = set(over_limit_ind)

        assert s0.difference(s1) == s1.difference(s0)

        removed = []
        sorted_atoms = [x[0] for x in Z]
        sorted_interactions = [x[1] for x in Z]
        print_err = True
        for index, value in enumerate(over_limit_ind):
            inds = [i for i, x in enumerate(sorted_atoms) if x == value]
            r_inds = iter(range(len(inds)))
            target = over_limit_num[index]
            target_interactions = [sorted_interactions[i] for i in inds]

            while target > 0:
                try:
                    i = next(r_inds)
                    removed.append(target_interactions[i])
                    G.remove_edge(target_interactions[i].atoms[0],
                                  target_interactions[i].atoms[1]
                                  )
                    ff.blocks[molname].remove_interaction('bonds', target_interactions[i].atoms)
                    target -= 1
                except nx.exception.NetworkXError:
                    if print_err == True:
                        print('Something went wrong while removing excess elastic network bonds. This is a placeholder statement while things are fixed.')
                    print_err = False
                    pass

    try:
        assert not any([i > 12 for i in [G.degree[node] for node in G.nodes]])
    except AssertionError:
        print("Couldn't remove all bonds necessary, probably can't network in VMD")

    with open(molname + '_surplus_en.txt', 'w') as extra_en:
        extra_en.write(f'Elastic network bonds removed from {molname}_en.itp\n')
        extra_en.write('This is for noting in visualisation, not for simulation\n\n')
        extra_en.write(f'These bonds will be missing if you load {molname}_en.itp in vmd\n')
        extra_en.write("having been present in your simulation. If you're inspecting your\n")
        extra_en.write('elastic network because you suspect some error because of it, bear this in mind.\n')

        extra_en.write('   i    j func b0 kb\n')

        for i in removed:
            extra_en.writelines(f'{i.atoms[0]:4d} {i.atoms[1]:4d} '+
                                f'{i.parameters[0]:1s} {i.parameters[1]:5s} {i.parameters[2]:5s}'+
                                '\n')

    # write the file out
    mol_out = ff.blocks[molname].to_molecule()
    mol_out.meta['moltype'] = molname + '_en'

    header = [f'Elastic network topology for {molname}', 'NOT FOR SIMULATIONS']

    with open(molname + '_en.itp', 'w') as outfile:
        write_molecule_itp(mol_out, outfile=outfile, header=header)


def topol_writing(topol_lines, ext='vis', W_include=None):
    '''

    Write new .top file based on the input one

    Parameters
    ----------
    topol_lines: dict
        lines from the input topology file split up into different keys, as per input_topol_reader
    ext: str
        the extension to use in looking for new file names to add to the output .top file
    W_include
        if W_include is not None (ie. args.system has been given something)
        then the line for water is not written out in the .top file

    Returns
    -------
    None
    '''
    # this will write everything and is gross but it should work
    # because everything will be included
    new_topol_head = []
    # make a new header for the .top file to write the absolute paths
    for i in topol_lines['core_itps']:
        new_topol_head.append(f'#include "{os.path.abspath(i)}"\n')
    # get the new vis files to write for the topology header
    vis_files = [i for i in os.listdir(os.getcwd()) if f'_{ext}' in i]
    for i in vis_files:
        new_topol_head.append(f'#include "{os.path.abspath(i)}"\n')

    # correct the [ molecules ] directive of the top file to correct for the new molecule names.
    topol_rest_vis = []
    original_mols = [i.split(f'_{ext}')[0] for i in vis_files]
    if W_include is not None:
        original_mols.remove('W')
    for i in topol_lines['molecules']:
        if any(mol in i['name'] for mol in original_mols):
            topol_rest_vis.append(i['name'] + f'_{ext}\t' + i['n_mols'] + '\n')
    # combine the sections of the vis.top and write it out.
    vis_topol = new_topol_head + topol_lines['system'] + topol_rest_vis

    with open(f'{ext}.top', 'w') as f:
        f.writelines(vis_topol)


def index_writer(system):
    '''
    Write a .ndx file for a system without water
    Parameters
    ----------
    system: str
        the .gro file to read and write the non-water index for

    Returns
    -------

    '''
    if os.path.splitext(system)[1] != '.gro':
        raise TypeError('Must provide a file in .gro format')

    print("Writing a waterless index file. Here're some helpful commands for reference:\n"
          "\tgmx trjconv -f traj_comp.xtc -s topol.tpr -pbc mol -n index.ndx -e 0 -o vis.gro\n"
          "\tgmx trjconv -f traj_comp.xtc -s topol.tpr -pbc mol -n index.ndx -o vis.xtc"
          )

    # read the file, take the number of the line if it doesn't have water in
    with open(f"{system}", 'r') as file:
        lines = file.readlines()
    # this skips the header and footer of the file
    l = [j for (j, k) in enumerate(lines[2:-1], start=1) if k[10:15].strip() != "W"]
    # split the lines every 12th index as gromacs requires
    lines_out = [l[x:x + 12] for x in range(0, len(l), 12)]

    # write the index file
    with open('index.ndx', 'w') as fout:
        fout.write('[ not water ]\n')
        for lineout in lines_out:
            fout.write(' '.join(map(str, lineout)) + ' \n')


def input_topol_reader(file):
    inclusions = []
    molecules = []
    system = []
    go = []
    # read the topology file, find the header lines which are not core martini itps
    with open(file) as f:
        mols = False
        sys = False
        for line in f.readlines():
            if '#include' in line:
                if 'go' not in line:
                    inclusions.append(line.split('"')[1])
                else:
                    go.append(line.split('"')[1])
            if mols == True:
                molecules.append({'name': line.split()[0],
                                  'n_mols': line.split()[1]})
            if 'molecules' in line:
                mols = True
                sys = False
                system.append(line)
            if 'system' in line:
                sys = True
            if sys == True:
                system.append(line)
    topol_lines = {'core_itps': inclusions,
                   'system': system,
                   'molecules': molecules}

    if len(go) > 0:
        topol_lines['go'] = go

    return topol_lines


if __name__ == '__main__':

    msg = ('Write out simplified molecule topologies to use cg_bonds-v5.tcl\n\n'
           'output .itp files will have the same name as the input ones, but\n'
           'only consist of [ moleculetype ], [ atoms ], and [ bonds ] directives.\n'
           'We also write a vis.top with correct paths and names for the new files.\n'
           'In VMD, first source cg_bonds.tcl, then run:\n'
           '\tcg_bonds -top vis.top\n'
           'and now you should be able to see your system with the CG bonds.'
           )

    parser = argparse.ArgumentParser(description=msg, formatter_class=RawTextHelpFormatter)
    parser.add_argument("-p", dest="topology", type=Path, help="input .top file used", default="topol.top")
    parser.add_argument("-f", dest="system", type=Path,
                        help=("Gromacs .gro file for which to write a non-water index file. "
                              "Equivalent to an index group of !W in gmx make_ndx. "
                              "Giving this option will automatically exclude W from your output vis.top")
                        )
    parser.add_argument("-vs", default=True, action="store_false", dest='virtual_sites',
                        help=("Don't write bonds between virtual sites and their constructing atoms. "
                              " (Bonds are written by default. Specify this flag if you don't want them written.)")
                        )
    parser.add_argument("-el", default=False, action="store_true", dest='elastic',
                        help="Write elastic network of input proteins to separate files"
                        )
    parser.add_argument("-ef", default=700, dest='en_force', type=float,
                        help="Force constant used for elastic network. Default = 700, standard for Martini 3."
                        )
    parser.add_argument("-go", default=False, action="store_true",
                        help="Go network options")
    parser.add_argument("-gf", type=Path, dest='go_path',
                        help="Nonbonded parameter itp file for your go network")

    args = parser.parse_args()

    # get files in the present directory
    files = os.listdir(os.getcwd())

    # get the topology file
    top = os.getcwd() + f'/{args.topology}'

    print(f"Reading input {args.topology}")
    topol_lines = input_topol_reader(top)

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
            misc_result = misc_file_reader(d[j])
            if misc_result is None:
                print(f"Error reading {topol_lines['core_itps'][i]}. Maybe you have a go model but didn't specify -go?")
                pass
            else:
                # this means we've separated things out successfully and can read the actual itp content now
                for key, value in misc_result[0].items():
                    system_defines[key] = value
                read_itp(misc_result[1], ff)

    # iterate over the molecules to make visualisation topologies
    keep = ['bonds', 'constraints', 'pairs', 'virtual_sitesn',
            'virtual_sites2', 'virtual_sites3']
    print("Writing visualisable topology files")
    for molname, block in ff.blocks.items():
        # write vis topols for molecules we're actually interested in
        try:
            assert molname in [i['name'] for i in topol_lines['molecules']]
        except AssertionError:
            continue

        # delete the interactions which are not bonds
        for interaction_type in list(block.interactions):
            if interaction_type not in keep:
                del block.interactions[interaction_type]
        # remove meta (ie. the #IFDEF FLEXIBLE) from the bonds
        for bond in block.interactions['bonds']:
            bond.meta.clear()
        if args.elastic:
            en_bonds = []
            for bond in list(block.interactions['bonds']):
                # these conditions are taken from the martini3001 forcefield, may cause upset in future
                # TODO monitor these conditions for future force fields
                # TODO work out how to insert the #defines statements into the molecule proper. Can't assign.
                if len(bond.parameters) == 1:
                    continue
                    # bond.parameters = system_defines[bond.parameters[0]]
                    # print(molname, bond.parameters, system_defines[bond.parameters[0]])

                cond0 = abs(float(bond.parameters[2]) - args.en_force) < 0.1  # elastic network proper
                cond1 = abs(float(bond.parameters[1]) - 0.970) < 0.1  # long beta elastic
                cond2 = abs(float(bond.parameters[1]) - 0.640) < 0.1  # short beta elastic
                if any([cond0, cond1, cond2]):
                    en_bonds.append(bond)
                    block.remove_interaction('bonds', bond.atoms)
            ff_en_copy = copy.deepcopy(ff)
            en_writer(ff_en_copy, molname, en_bonds)

        # this should then keep any constraints which don't have IFDEF statements
        # eg. alpha helices are described by constraints without these.
        # however, the remove_interactions function doesn't work atm.
        for bond in list(block.interactions['constraints']):
            if bond.meta.get('ifndef'):
                block.remove_interaction('constraints', bond.atoms)

        # rewrite pairs as bonds for visualisation
        for bond in block.interactions['pairs']:
            block.add_interaction('bonds', bond.atoms, bond.parameters[:2] + ['10000'])

        del block.interactions['pairs']

        # make bonds between virtual sites and each of the constructing atoms
        go_dict = {}
        for vs_type in ['virtual_sitesn', 'virtual_sites2', 'virtual_sites3']:
            if args.virtual_sites:
                for vs in block.interactions[vs_type]:
                    site = vs.atoms[0]
                    constructors = vs.atoms[1:]
                    # this avoids pointless bonds between a virtual site directly on top of
                    # its singular constructing atom
                    block.nodes[site]['mass'] = 1
                    if args.go:
                        # make a dictionary of atype: node index
                        # this is for later so the 'bond' can be drawn properly.
                        # assert that this site has the name CA to check its a go site and not another VS.
                        aname = block.nodes[site]['atomname']
                        atype = block.nodes[site]['atype']
                        if aname == 'CA':
                            go_dict[atype] = site
                    else:
                        for constructor in constructors:
                            # completely arbitrary parameters, the bond just needs to exist
                            block.add_interaction('bonds', [site, constructor],
                                                  ['1', '1', '1000'])
            if block.interactions[vs_type]:
                del block.interactions[vs_type]

        if args.go:
            if args.go_path:
                go_nb_file = args.go_file
            else:
                go_nb_file = "go_nbparams.itp"

            if not os.path.isfile(go_nb_file):
                raise FileNotFoundError("Gō nonbonded itp does not exist. Specify using -gf")

            with open(go_nb_file, "r") as f:
                nb_lines = f.readlines()
            nb_lines = [line.split(';')[0].split(' ') for line in nb_lines if '[' not in line]
            bonds_list = [[go_dict[i[0]], go_dict[i[1]], '1', i[3], '1000'] for i in nb_lines]
            ff_go_copy = copy.deepcopy(ff)
            go_writer(ff_go_copy, molname, bonds_list)

        # write out the molecule with an amended name
        mol_out = block.to_molecule()
        mol_out.meta['moltype'] = molname + '_vis'

        header = [f'Visualisation topology for {molname}', 'NOT FOR SIMULATIONS']

        with open(molname + '_vis.itp', 'w') as outfile:
            write_molecule_itp(mol_out, outfile=outfile, header=header)

    if args.elastic:
        topol_writing(topol_lines, 'en', W_include=args.system)
    if args.go:
        topol_writing(topol_lines, 'go', W_include=args.system)
    topol_writing(topol_lines, W_include=args.system)

    if args.system is not None:
        index_writer(args.system)

    print('All done!')
