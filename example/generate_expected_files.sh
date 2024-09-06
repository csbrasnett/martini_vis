#!/bin/bash

source /your/martini_vis/venv/bin/activate
source /your/gromacs/installation/bin/GMXRC

wget https://raw.githubusercontent.com/marrink-lab/martini-forcefields/main/martini_forcefields/regular/v3.0.0/gmx_files/martini_v3.0.0.itp -O martini.itp

wget https://files.rcsb.org/download/1ubq.pdb

grep "^ATOM" 1ubq.pdb > 1UBQ_clean.pdb

martinize2 -f 1UBQ_clean.pdb -o topol.top -x 1UBQ_cg.pdb -dssp /usr/local/bin/mkdssp -p backbone -ff martini3001 -elastic -ef 700.0 -el 0.5 -eu 0.9 -ea 0 -ep 0 -scfix -cys auto -maxwarn 1

gmx editconf -f 1UBQ_cg.pdb -c -d 5 -bt dodecahedron -o newbox.gro

martini_vis -p topol.top -el -ef 700 -vf
