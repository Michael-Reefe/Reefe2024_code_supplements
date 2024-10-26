Description of contents:

The analysis code repository can be found at https://github.com/Michael-Reefe/Loki.jl
- This contains the full source code, an example jupyter notebook to test the code on, and a README file with 
  information on the installation, dependencies, and expected outputs of the code.
- ***IMPORTANT***: The analysis done in the paper was performed using v2.0 of the code. To make sure the code provided here
  runs, make sure you change the code to the 'v2.0' branch by running `git checkout v2.0` after cloning the repository.
- First install julia 1.10 or later from https://julialang.org/downloads/. Then follow the code
  installation instructions under Loki.jl/README.md to install the Loki.jl package within its own project.
  The scripts included in this folder will assume that it has been installed this way.

The example notebook 'example.ipynb' is included under the examples/ directory of the github repository.
- Once the code has been installed, this example notebook can be used to test that it's working.
  Just open up the jupyter notebook with the julia kernel and try to run it. Note that you will need to provide a 
  sample reduced JWST data file yourself for the code to run on.
- To run the notebook, you can start julia, activate the Loki.jl project, and use the following commands:
    ```
    julia> using IJulia
    julia> notebook()
    ```
- Then navigate to the notebook and run it.  Alternatively, the VSCode extensions for jupyter + julia work well together.
- This is a short, quick example on a small data set and is only meant to show the general capability of the
  code. On my 2020 M1 Macbook Pro, it takes about 1.5 minutes to run all of the code blocks.

To run the code on our actual data and reproduce the resultant 2D flux & velocity maps seen in the manuscript:

- Download the data from JWST program ID 2439 from MAST (https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html)
- Note: you will need to download the stage 1 UNCALIBRATED data and perform the data reduction using the attached
  notebooks to match our results.
- Place the uncalibrated fits files in a folder, splitting them up by the source data and background data (which 
  you can tell from the target name in the fits headers)
- First run MRS_FlightNB1.ipynb to perform the main data reduction, then run Correct_MIRI_Cubes.ipynb to apply the 
  scaling factor correction.  Once the data is fully reduced, you are ready to run our analysis scripts.

The main code driver file is 'fit_phoenix_single_channel.jl'. This file will run fits to each individual spaxel for
  a single channel of the MIRI/MRS data, separating the host and AGN components. Before running this file, there are
  a few things to know, outlined below.
- You will need to edit the file paths accordingly to point to where you have placed your reduced data files.
- 'ext_template_phoenix.txt' contains the template for the extinction profile shape as a function of wavelength that we 
  use. This template was generated from a previous fit to the integrated spectrum over all channels with the 'extinction_curve'
  option set to "decompose".
- 'PHOENIX.EXTINCTION.FROM_H2.BLUR1.CSV' contains a map of the optical depth at 9.7 um (tau_9.7) in each spaxel in the
  channel 3 cube. This template was generated from a previous iteration of fits by looking at the ratio of the H2 S(4) and 
  S(3) lines.
- Finally, to run the file and fit the channel 3 data, use the command `julia --project=Loki.jl fit_phoenix_single_channel.jl 3`. 
- To run fits for other channels, replace the '3' with '1', '2', or '4'. Refer to the README file in the Loki repository 
  for information about the output folders and files. You must run channel 3 before the other channels since the channel 3
  WCS parameters created in the input cubes are used to reproject PHOENIX.EXTINCTION.FROM_H2.BLUR1.CSV onto other
  channels.
- The channel 3 fits will likely take a few hours to complete. On my 2020 M1 Macbook Pro, it took ~4 hours with 4 parallel processes.

There is an additional code driver, 'fit_phoenix_multi_channel.jl' that will perform individual spaxel fits on the combined
  channel 2+3+4 cube. Most of the above notes also apply here, including the requirement that the individual channel 3 fit
  be performed first to generate the Phoenix_Cluster_BCG.channel3.rest_frame.fits file for reprojecting the optical depth map.
  To run this, use `julia --project=Loki.jl fit_phoenix_multi_channel.jl`. This took about 2 hours for me (shorter than the
  single-channel runs because global minimization on each spaxel has been disabled).

Finally, there is 'fit_phoenix_multi_channel_aperture.jl', which performs a fit to the combined channel 2+3+4 cube integrated
  over the northern aperture. This fit locks the AGN PSF template amplitudes to the values found in the individual spaxel fits,
  so it requires 'fit_phoenix_multi_channel.jl' to be run first. This will run 100 bootstrap iterations, which can be sped up by
  allowing multithreading by adding the threads option, i.e.:  `julia --project=Loki.jl --threads=4 fit_phoenix_multi_channel_aperture.jl`.


Feel free to email me with any inquiries:
mreefe (at) mit (dot) edu

