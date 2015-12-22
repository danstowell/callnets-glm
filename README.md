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

Tested on Ubuntu Linux 14.04.


How to Run the Code
===================

1. Download the zf4f data files:

             TODO **********************************************************

2. Prepare the mex code (only needs doing once).

         cd code_GLM/tools_mexcode
         octave
         initialize_mexcode

3. Invoke the main Octave script, which will run the GLM analysis multiple times on subsets of the zf4f data.

         cd code_GLM
         source runme.octave

4. Invoke the main Python script, which will run GLM and cross-correlation analysis on real and simulated data, and make various plots.

         cd callnets
         python callnets.py

    You can then also run the auxiliary scripts which do some of the tests and figures in the paper:

         python nonlins.py  # just plots what the nonlinearities look like
         python oddsratios.py  # calculates odds-ratios between softplus and exp models
         python analyseadj.py  # analyses the predictability of one segment from the previous


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

