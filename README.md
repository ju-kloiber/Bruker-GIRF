# Bruker-GIRF
Repository for loading gradient shapes from Bruker and computing GIRFs

# Compute GIRF
This module can be used to compute GIRF from Bruker data. An example can be found in the subfolder /examples. After importing necessary modules, a few parameters and variables need to be specified. These parameters are stored in a dictionary and include:

* experiment: name of the current experiment. Used to load Bruker data and store computed outputs. 
* folder: subfolder containing the Bruker data (jcamp files).
* files: array containing all files names used to compute the GIRF. Each file represents one gradient amplitude. 
* amps: array containing the gradient amplitudes used in the experiment. First amplitude corresponds to first file in "files" etc.
* p: slew rate in mT/m/ms used in the experiment.
* tDwell: dwell time of the acquistion in ms.
* nPoints: number of points in the jcamp files.
* gamma: gyromagnetic ratio in MHz/T.
* BW: bandwidth in kHz used for the raised cosine filter.
* α: roll-off factor of the raised cosine filter. 

The gradients can be computed from the k-space trajectories stored in the jcamp files. The function load_k() loads the k-space trajectories in all three spatial directions as well as the time grid. The gradient are obtained via the function compute_G().

The nominal gradients are specified by the amplitudes and slew rate. The function construct_nominal_G() can be used to construct triangular input shapes.

The GIRF is computed via the function compute_girf(). It returns the GIRF in time and frequency domain as well as the corresponding timepoints and frequencies. 

The function raised_cosine() constructs a raised cosine filter on the given frequency interval. It can be applied by pointwise multiplication with the GIRF. 

# Apply GIRF
When a GIRF and nominal gradients are available, GIRF-corrected gradients can be computed. The function apply_girf() applies the GIRF to the nominal gradients. 

