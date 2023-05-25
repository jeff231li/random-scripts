# Charging a Cyclodextrin
This python module allows one to assign partial charges on a whole cyclodextrin molecule efficiently using a fragment-based approach. The method was adapted from [wutobias](https://github.com/wutobias/collection/blob/master/python/am1bcc-glyco-fragmenter.py)'s code. Basically, the user specifies the monomer fragments and the code will break the cyclodextrin molecule into smaller pieces according to the matching fragments. The code will cap the broken bonds with a methyl and then assign AM1-BCC charges. The smaller pieces are stitched back to the cyclodextrin without the methyl caps and normalizes the charge. 

I have included two examples here,
1) [01-charging_cyclodextrin_derivatves.ipynb](01-charging_cyclodextrin_derivatves.ipynb) - charging a beta-cyclodextrin with one monomer as a derivative (2,6-di-O-methyl)
2) [02-charging_cyclodextrin_linked_dimers.ipynb](02-charging_cyclodextrin_linked_dimers.ipynb) - charging a covalently-linked beta-cyclodextrin dimer (thiol linker).

### Dependencies for the code:
* RDKit
* OpenFF-Toolkit
* AmberTools
