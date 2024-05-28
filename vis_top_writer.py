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
    #remove all interactions from the molecule
    for interaction_type in list(ff.blocks[molname].interactions):
        del ff.blocks[molname].interactions[interaction_type]

    #add the elastic network bonds back in
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
    #remove all interactions from the molecule
    for interaction_type in list(ff.blocks[molname].interactions):
        del ff.blocks[molname].interactions[interaction_type]

    #add the elastic network bonds back in
    for bond in en_bonds:
        ff.blocks[molname].add_interaction('bonds', bond.atoms, bond.parameters)

    # write the file out
    mol_out = ff.blocks[molname].to_molecule()
    mol_out.meta['moltype'] = molname + '_en'

    header = [f'Elastic network topology for {molname}', 'NOT FOR SIMULATIONS']

    with open(molname + '_en.itp', 'w') as outfile:
        write_molecule_itp(mol_out, outfile=outfile, header=header)

def topol_writing(topol_core_itps, topol_rest, ext='vis'):
    '''

    Write new .top file based on the input one

    Parameters
    ----------
    topol_core_itps: list
        itp files #included in the input .top file with 'martini' in their names
    topol_rest: list
        itp files #included in the input .top file without 'martini' in their names
    ext: str
        the extension to use in looking for new file names to add to the output .top file

    Returns
    -------
    None
    '''
    # make a new header for the .top file to write the absolute paths
    new_topol_head = []
    for i in topol_core_itps:
        new_topol_head.append(f'#include "{os.path.abspath(i)}"\n')

    # get the new vis files to write for the topology header
    vis_files = [i for i in os.listdir(os.getcwd()) if f'_{ext}' in i]
    for i in vis_files:
        new_topol_head.append(f'#include "{os.path.abspath(i)}"\n')

    # correct the [ molecules ] directive of the top file to correct for the new molecule names.
    topol_rest_vis = []
    original_mols = [i.split(f'_{ext}')[0] for i in vis_files]
    for i in topol_rest:
        if any(mol in i for mol in original_mols):
            topol_rest_vis.append(i.split()[0] + f'_{ext}\t' + i.split()[1] + '\n')
        else:
            if len(i.split()) > 1:
                if args.system and i.split()[0] == 'W':
                    pass
                else:
                    topol_rest_vis.append(i)
            else:
                topol_rest_vis.append(i)
    # combine the sections of the vis.top and write it out.
    vis_topol = new_topol_head + topol_rest_vis

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
    if os.path.splitext(system)[1]!='.gro':
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
    parser.add_argument("-p", dest="topology", type=str, help="input .top file used", default="topol.top")
    parser.add_argument("-f", dest="system", type=str,
                        help=("Gromacs .gro file for which to write a non-water index file. "
                              "Equivalent to an index group of !W in gmx make_ndx. "
                              "Giving this option will automatically exclude W from your output vis.top")
                        )
    parser.add_argument("-s", default=True, action="store_false", dest='virtual_sites',
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

    args = parser.parse_args()

    #assume user error and be forgiving
    if args.en_force and not args.elastic:
        if not args.go:
            args.elastic = True

    #get files in the present directory
    files = os.listdir(os.getcwd())
    
    #get the topology file
    top = os.getcwd() + f'/{args.topology}'

    print(f"Reading input {args.topology}")

    #read the topology file, find the header lines which are not core martini itps
    with open(top) as f:
        topol_lines = f.readlines()

    topol_core_itps = [a.split('"')[1] for a in topol_lines if "martini" in a]
    topol_rest = [a for a in topol_lines if "#include" not in a]
    if not args.go:
        topol_editable_itps = [a.split('"')[1] for a in topol_lines if "#include" in a and "martini" not in a]
    else:
        topol_editable_itps = [a.split('"')[1] for a in topol_lines if "#include" in a and "martini" not in a and "go" not in a]
        topol_go_itps = [a.split('"')[1] for a in topol_lines if "go" in a]

    #for each molecule in the system, read in the itp
    d = {}
    for i,j in enumerate(topol_editable_itps):
        with open(j) as f:
            d[i] = f.readlines()
    
    #read the molecules into the forcefield
    ff = ForceField('martini3001')
    
    for i in d.keys():
        read_itp(d[i], ff)
    
    #iterate over the molecules to make visualisation topologies 
    keep = ['bonds', 'constraints', 'pairs', 'virtual_sitesn',
            'virtual_sites2', 'virtual_sites3']
    print("Writing visualisable topology files")
    for molname, block in ff.blocks.items():
        #delete the interactions which are not bonds
        for interaction_type in list(block.interactions):
            if interaction_type not in keep:
                del block.interactions[interaction_type]
        # remove meta (ie. the #IFDEF FLEXIBLE) from the bonds
        for bond in block.interactions['bonds']:
            bond.meta.clear()
        if args.elastic:
            en_bonds = []
            for bond in list(block.interactions['bonds']):
                cond0 = abs(float(bond.parameters[2])-args.en_force)<0.1 #elastic network proper
                cond1 = abs(float(bond.parameters[1]) - 0.970)< 0.1 # long beta elastic
                cond2 = abs(float(bond.parameters[1]) - 0.640)< 0.1 # short beta elastic
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

        #rewrite pairs as bonds for visualisation
        for bond in block.interactions['pairs']:
            block.add_interaction('bonds', bond.atoms, bond.parameters[:2] + ['10000'])
        
        del block.interactions['pairs']

        #make bonds between virtual sites and each of the constructing atoms
        go_dict = {}
        for vs_type in ['virtual_sitesn', 'virtual_sites2', 'virtual_sites3']:
            if args.virtual_sites:
                for vs in block.interactions[vs_type]:
                    site = vs.atoms[0]
                    constructors = vs.atoms[1:]
                    for constructor in constructors:
                        # completely arbitrary parameters, the bond just needs to exist
                        block.add_interaction('bonds', [site, constructor],
                                              ['1', '1', '10000'])
                    if args.go:
                        #make a dictionary of atype: node index
                        #this is for later so the 'bond' can be drawn properly.
                        #assert that this site has the name CA to check its a go site and not another VS.
                        aname = block.nodes[site]['atomname']
                        atype = block.nodes[site]['atype']
                        if aname == 'CA':
                            go_dict[atype] = site

            del block.interactions[vs_type]

        if args.go:
            if "go_nbparams.itp" not in topol_go_itps:
                print('need go_nbparams.itp to use this')
            with open("go_nbparams.itp", "r") as f:
                nb_lines = f.readlines()
            nb_lines = [line.split(';')[0].split(' ') for line in nb_lines if '[' not in line]
            bonds_list = [(go_dict[i[0]], go_dict[i[1]], '1', i[3], '1000') for i in nb_lines]
            ff_go_copy = copy.deepcopy(ff)
            go_writer(ff_go_copy, molname, bonds_list)

        #write out the molecule with an amended name
        mol_out = block.to_molecule()
        mol_out.meta['moltype'] = molname+'_vis'

        header = [f'Visualisation topology for {molname}', 'NOT FOR SIMULATIONS']

        with open(molname+'_vis.itp', 'w') as outfile:
            write_molecule_itp(mol_out, outfile=outfile, header = header)

    topol_writing(topol_core_itps, topol_rest)
    if args.elastic:
        topol_writing(topol_core_itps, topol_rest, 'en')
    if args.go:
        topol_writing(topol_core_itps, topol_rest, 'go')

    if args.system:
        index_writer(args.system)

    print('All done!')
