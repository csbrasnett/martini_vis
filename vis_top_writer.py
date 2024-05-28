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

def en_remove(ff, molname, en_bonds):
    for interaction_type in list(ff.blocks[molname].interactions):
        del ff.blocks[molname].interactions[interaction_type]

    for bond in en_bonds:
        ff.blocks[molname].add_interaction('bonds', bond.atoms, bond.parameters)
    mol_out = ff.blocks[molname].to_molecule()
    mol_out.meta['moltype'] = molname + '_en'

    header = [f'Elastic network topology for {molname}', 'NOT FOR SIMULATIONS']

    with open(molname + '_en.itp', 'w') as outfile:
        write_molecule_itp(mol_out, outfile=outfile, header=header)


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

    args = parser.parse_args()

    #assume user error and be forgiving
    if args.en_force and not args.elastic:
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
    topol_editable_itps = [a.split('"')[1] for a in topol_lines if "#include" in a and "martini" not in a]
    
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
            ff_copy = copy.deepcopy(ff)
            en_remove(ff_copy, molname, en_bonds)

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
        for vs_type in ['virtual_sitesn', 'virtual_sites2', 'virtual_sites3']:
            if args.virtual_sites:
                for vs in block.interactions[vs_type]:
                    site = vs.atoms[0]
                    constructors = vs.atoms[1:]
                    for constructor in constructors:
                        # completely arbitrary parameters, the bond just needs to exist
                        block.add_interaction('bonds', [site, constructor],
                                              ['1', '1', '10000'])
            del block.interactions[vs_type]

        #write out the molecule with an amended name
        mol_out = block.to_molecule()
        mol_out.meta['moltype'] = molname+'_vis'
        
        header = [f'Visualisation topology for {molname}', 'NOT FOR SIMULATIONS']
        
        with open(molname+'_vis.itp', 'w') as outfile:
            write_molecule_itp(mol_out, outfile=outfile, header = header)
    
    #make a new header for the .top file to write the absolute paths
    new_topol_head = []
    for i in topol_core_itps:
        new_topol_head.append(f'#include "{os.path.abspath(i)}"\n')

    #get the new vis files to write for the topology header
    vis_files = [i for i in os.listdir(os.getcwd()) if '_vis' in i]
    for i in vis_files:
        new_topol_head.append(f'#include "{os.path.abspath(i)}"\n')
    
    #correct the [ molecules ] directive of the top file to correct for the new molecule names.
    topol_rest_vis = []
    original_mols = [i.split('_vis')[0] for i in vis_files]
    for i in topol_rest:
        if any(mol in i for mol in original_mols):
            topol_rest_vis.append(i.split()[0]+'_vis\t' + i.split()[1]+'\n')
        else:
            if len(i.split())>1:
                if args.system and i.split()[0] == 'W':
                    pass
                else:
                    topol_rest_vis.append(i)
            else:
                topol_rest_vis.append(i)
    #combine the sections of the vis.top and write it out.
    vis_topol = new_topol_head + topol_rest_vis
    
    with open('vis.top', 'w') as f:
        f.writelines(vis_topol)

    if args.system:
        print("Writing a waterless index file. Here're some helpful commands for reference:\n"
              "\tgmx trjconv -f traj_comp.xtc -s topol.tpr -pbc mol -n index.ndx -e 0 -o vis.gro\n"
              "\tgmx trjconv -f traj_comp.xtc -s topol.tpr -pbc mol -n index.ndx -o vis.xtc"
              )
        #read the file, take the number of the line if it doesn't have water in
        with open(f"{args.system}",'r') as file:
            lines = file.readlines()
        #this skips the header and footer of the file
        l = [j for (j,k) in enumerate(lines[2:-1], start = 1) if k[10:15].strip()!="W"]
        #split the lines every 12th index as gromacs requires
        lines_out = [l[x:x+12] for x in range(0, len(l),12)]

        #write the index file
        with open('index.ndx', 'w') as fout:
            fout.write('[ not water ]\n')
            for lineout in lines_out:
                fout.write(' '.join(map(str, lineout))+' \n')
    print('All done!')
