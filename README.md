# martini_vis
![Martinizing a protein (Ubiquitin, PBD: 1UBQ) and visualising its elastic network.
left to right: atomistic representation of the protein. 
Martini representation of the protein overlaid on the atomistic one, showing the direct backbone and side chain.
Martini representation of the protein overlaid on the atomistic one, showing all bonds, with the elastic network in black.
Martini representation of the protein, showing all bonds.](/grain/chris/martini_vis/image.png "Visualising elastic networks")

Scripts to aid the visualisation of coarse-grained martini trajectories.

Martini_vis centres on the `vis_top_writer.py` script, which uses [vermouth](https://github.com/marrink-lab/vermouth-martinize) to stably
rewrite your input topology files as ones that can be used for visualisation in vmd.

This builds on previous work that produced the `cg_bonds-v5.tcl` script, which reads martini topology information into vmd. At some point
`cg_bonds-v5.tcl` wasn't working with the latest models, so I wrote a script that should sort everything out. 

If the solution here isn't working for you, please open an issue!

## Dependencies

The only non-standard library used to run `vis_top_writer.py` is [Vermouth](https://github.com/marrink-lab/vermouth-martinize). 
Please ensure you have Vermouth installed in your Python environment before running the script.

## Usage

1) Run `./vis_top_writer.py` with your .top file you used to run a simulation. This will produce:
   1) Edited *_vis.itp files for all the non-standard (e.g. water, ions) molecules in your system described in the input .top file.  
      * NB. by default, virtual sites will be rewritten as "real", with bonds between the sites and their constructing atoms. 
      To stop this, use the `-s` flag when running `vis_top_writer.py`.
   2) `vis.top`, a new .top file for your system and the visualisable topologies. `cg_bonds-v5.tcl` requires absolute paths to your itps, which is solved by running the script.
   3) Optionally by providing the .gro file you plan to visualise, you can write an index file without containing your system without water to use in processing your trajectory. 
   Something like `gmx trjconv -f traj_comp.xtc -s topol.tpr -n index.ndx -pbc mol -o vis.xtc` will write new trajectory using the index file provided. As there is only one index group,
   no further interaction with `trjconv` is required. 
      * NB. if you use this option, then `vis.top` will not contain an entry for the waters in your system at all.
2) Load your simulation into vmd.
3) `source cg_bonds-v5.tcl` in vmd.
4) Load your visualisable topologies using `cg_bonds-v5 -top vis.top`.
5) Visualise your simulation with bonds in Martini!

## Notes on using `vis_top_writer.py`

### Input
We assume the input topology looks like something as follows:

```
#include "martini_v3.0.0.itp"
#include "martini_v3.0.0_solvents_v1.itp"
#include "martini_v3.0.0_ions_v1.itp"
#include "my_protein.itp"
#include "something_else.itp"

[ system ]
system name

[ molecules ]
Protein 10
W         1000
something_else   10
NA               10
CL               10
```

i.e. you have :
1) a set of input topology files from Martini with "martini" in their names
2) a set of input topology files that are specific to your system, and __don't__ have "martini" anywhere in their names

`vis_top_writer.py` will maintain things that are generic to martini systems (like molecules that are defined in set 1: water, ions, etc.) and will make visualisable topologies for anything it finds in set 2 (like your protein and other molecules)

### Output

Running `./vis_top_writer.py -p topol.top -f frame.gro` on the above system, together with a .gro file that you want an index file for
will produce the following output:

```
#include "/absolute/path/to/file/martini_v3.0.0.itp"
#include "/absolute/path/to/file/martini_v3.0.0_solvents_v1.itp"
#include "/absolute/path/to/file/martini_v3.0.0_ions_v1.itp"
#include "/absolute/path/to/file/my_protein_vis.itp"
#include "/absolute/path/to/file/something_else_vis.itp"

[ system ]
system name

[ molecules ]
Protein_vis 10
something_else_vis   10
NA               10
CL               10
```
So the differences are:
1) The input topologies are:
   * Now written with absolute paths to the directory they were found in
   * For non-default martini topologies (ie. the molecules of interest in your system)
   new visualisable topologies have been written and included
2) The [ molecules ] directive of the file:
   * Lists the new names of the visualising topologies
   * No longer has an entry for the water in the system, because a water-less index file was written.

