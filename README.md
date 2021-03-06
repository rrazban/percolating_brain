# percolating brain
Percolation theory is applied to diffusion MRI scans from human individuals.

`process/` contains the script to extract white matter tract lengths and densities from diffusion MRI. 
Sample outputs from `process/dmri2adjacency_matrix.py` are provided in `sample_outputs/`
* note that raw diffusion MRI scans across individuals are not provided in this repository and must be obtained directly from the respective databank

`analyze/` contains scripts to study the percolation probability curves generated by targeted attack 
on tract length or density extracted from diffusion MRI.

`simulation/` containts scripts to assess the mechanism underlying the percolation probability curves.
