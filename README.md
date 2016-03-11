Call Networks Analysis Through GLM Point Process Modelling
==========================================================
By Dan Stowell.
Includes a modified copy of "code_GLM" by Jonathan Pillow.


Software Requirements
=====================
Here we give the software we used to run this code. Where version numbers are given in brackets this
means it's the version we used, not necessarily the exact requirement.

* Octave (3.8.1)
   * It should work with Matlab too (untested)

* Python 2.7.x
   * With the following modules:
      * numpy (1.8.2)
      * scipy (0.13.3)
      * matplotlib (1.3.1)
      * subprocess32

* pdflatex, and the "tikz" latex package (to draw social network plots)

Tested on Ubuntu Linux 14.04.

It should run on Linux and Mac fine. Untested on Windows.


How to Run the Code for the Tests Used in the Paper
===================================================

1. Download the data files: you can do this by running the following command:

         python download_data.py

      If you have any problems downloading, these URLs give direct access:

      zf4f: https://figshare.com/articles/Zebra_finch_group_calling_zf4f/1613791

      gill et al: https://figshare.com/s/73cbdf96ab156a0f0a69

2. Prepare the mex code (only needs doing once).

         cd code_GLM/tools_mexcode
         octave --eval "initialize_mexcode"
         cd ../..

3. Invoke the main Octave script, which will run the GLM analysis multiple times on subsets of the zf4f data.
   THIS WILL TAKE A WHILE TO RUN! On my laptop, to analyse the zf4f takes about half an hour, then to analyse
   the Gill data takes overnight.

         cd code_GLM
         source runme.octave
         cd ..

4. Invoke the main Python script, which will run GLM and cross-correlation analysis on real and simulated data, and make various plots.

         cd callnets
         python callnets.py

    (By default the synthetic ABC test only makes a tiny amount of data, since otherwise the test
    takes hours. You can tweak the duration directly in the code, to create plots more like the paper.)

    You can then also run the auxiliary scripts which do some of the tests and figures in the paper:

         python nonlins.py  # just plots what the nonlinearities look like
         python oddsratios.py  # calculates odds-ratios between softplus and exp models
         python analyseadj.py  # analyses the predictability of one segment from the previous

How to Adapt the Code for Your Own Use
======================================

To run the code on your own data, have a look at the Octave/Matlab code in `code_GLM/testscripts/zf4f_glm_each.m`, which simply iterates through a set of datafiles calling `dofit_fromcsv_GLM_zf4f()` for each one. You can simply call `dofit_fromcsv_GLM_zf4f()` yourself, have a look at its parameters. If you specify a `csvoutdir` it outputs data files which are useful for inspecting the results (as is done by `callnets.py`).

The data format we use for call data is simple CSV (Comma Separated Values), with three columns:

1. start time, in seconds
2. duration, in seconds (not used! you can just fix this as zero)
3. which individual called, an integer from 1 to K

We do not include a header row in the CSV. Have a look at one of the CSV files in the `data` folder to see examples of the format.

Once the code has run, it generates output CSVs and plots. If you look in the `outcsv` folder you should see various files, including summary statistics. The file with name ending in `_kernels_discret.csv` holds the discretised kernel shapes, which you should be able to load and plot/analyse very easily.

Copyright and Licence
=====================
The code in the "code_GLM" folder is mostly (c) 2010 Jonathan Pillow with additions (c) 2015 Dan Stowell,
and it is licensed under the GPL3+ as shown in /code_GLM/license/LICENSE.txt

The code in the "callnets" folder is (c) 2015 Dan Stowell,
and it is licensed under the GPL2+.


Changes made to code_GLM
========================
The "downstream" changes I have made to code_GLM include:

* Adding a regularisation parameter to the fitting (added to `MLfit_GLM()`)
* Adding a parametric softplus nonlinearity function
* Changing parameters such as the size and shape of the basis kernels for the call-analysis case
* Adding convenience functions for analysing CSV call data, and for sampling from a fitted model
* Small change to `makeStimRows()` to allow for size-1 vectors.

