# martini_vis

Scripts to help visualise coarse-grained martini topologies

## Usage

1) Run `./vis_top_writer.py` with your .top file you used to run a simulation. This will produce:
    * edited *_vis.itp files for all the non-standard (e.g. water, ions) molecules in your system described in the input .top file.  
    * `vis.top`, a new .top file for your system and the visualisable topologies.
2) Load your simulation into vmd.
3) `source cg_bonds-v5.tcl` in vmd.
4) Load your visualisable topologies using `cg_bonds-v5 -top vis.top`.
5) Visualise your simulation with bonds in Martini!

# Notes on using `vis_top_writer.py`

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
