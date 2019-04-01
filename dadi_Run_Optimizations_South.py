import sys
import os
import numpy
import dadi
import pylab
from datetime import datetime
import Optimize_Functions
import Demographics1D

'''
Usage: python dadi_Run_Optimizations.py

This is meant to be a general use script to run dadi to fit any model on an
afs/jsfs with one to three populations. The user will have to edit information about
their allele frequency spectrum and provide a custom model. The instructions are
annotated below, with a #************** marking sections that will have to be edited.
Several examples of how to use various arguments to control optimizations are shown.

This script must be in the same working directory as Optimize_Functions.py, which
contains all the functions necessary.

If you'd like to use this script for larger sets of models already available, please
look on the other repositories to see how to import models from external model scripts.

General workflow:
 The optimization routine runs a user-defined number of rounds, each with a user-defined
 or predefined number of replicates. The starting parameters are initially random, but after
 each round is complete the parameters of the best scoring replicate from that round are
 used to generate perturbed starting parameters for the replicates of the subsequent round.
 The arguments controlling steps of the optimization algorithm (maxiter) and perturbation
 of starting parameters (fold) can be supplied by the user for more control across rounds.
 The user can also supply their own set of initial parameters, or set custom bounds on the
 parameters (upper_bound and lower_bound) to meet specific model needs. This flexibility
 should allow these scripts to be generally useful for model-fitting with any data set.

Outputs:
 For each model run, there will be a log file showing the optimization steps per replicate
 and a summary file that has all the important information. Here is an example of the output
 from a summary file, which will be in tab-delimited format:

 Model	Replicate	log-likelihood	AIC	chi-squared	theta	optimized_params(nu1, nu2, m, T)
 sym_mig	Round_1_Replicate_1	-1684.99	3377.98	14628.4	383.04	0.2356,0.5311,0.8302,0.182
 sym_mig	Round_1_Replicate_2	-2255.47	4518.94	68948.93	478.71	0.3972,0.2322,2.6093,0.611
 sym_mig	Round_1_Replicate_3	-2837.96	5683.92	231032.51	718.25	0.1078,0.3932,4.2544,2.9936
 sym_mig	Round_1_Replicate_4	-4262.29	8532.58	8907386.55	288.05	0.3689,0.8892,3.0951,2.8496
 sym_mig	Round_1_Replicate_5	-4474.86	8957.72	13029301.84	188.94	2.9248,1.9986,0.2484,0.3688

Notes/Caveats:
 The likelihood and AIC returned represent the true likelihood only if the SNPs are
 unlinked across loci. For ddRADseq data where a single SNP is selected per locus, this
 is true, but if SNPs are linked across loci then the likelihood is actually a composite
 likelihood and using something like AIC is no longer appropriate for model comparisons.
 See the discussion group for more information on this subject.

Citations:
 If you use these scripts for your work, please cite the following publication:
    Portik, D.M., Leach, A.D., Rivera, D., Blackburn, D.C., Rdel, M.-O.,
    Barej, M.F., Hirschfeld, M., Burger, M., and M.K.Fujita. 2017.
    Evaluating mechanisms of diversification in a Guineo-Congolian forest
    frog using demographic model selection. Molecular Ecology 26: 52455263.
    doi: 10.1111/mec.14266

-------------------------
Written for Python 2.7
Python modules required:
-Numpy
-Scipy
-dadi
-------------------------

Daniel Portik
daniel.portik@gmail.com
https://github.com/dportik
Updated May 2018
'''

#===========================================================================
# Import data to create joint-site frequency spectrum
#===========================================================================

#**************
snps = "../final_snps_for_dadi.recode.vcf.data"

#Create python dictionary from snps file
dd = dadi.Misc.make_data_dict(snps)

#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=["South"]

#**************
#projection sizes, in ALLELES not individuals
proj = [30]

#Convert this dictionary into folded AFS object
#[polarized = False] creates folded spectrum object
fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = True)

#print some useful information about the afs or jsfs
print "\n\n============================================================================\nData for site frequency spectrum\n============================================================================\n"
print "projection", proj
print "sample sizes", fs.sample_sizes
sfs_sum = numpy.around(fs.S(), 2)
print "Sum of SFS = ", sfs_sum, '\n', '\n'

#================================================================================
# Here is an example of using a custom model within this script
#================================================================================
'''
 We will use a function from the Optimize_Functions.py script:

 Optimize_Routine(fs, pts, outfile, model_name, func, rounds, param_number, fs_folded=True, reps=None, maxiters=None, folds=None, in_params=None, in_upper=None, in_lower=None, param_labels=" ")

   Mandatory Arguments =
    fs:  spectrum object name
    pts: grid size for extrapolation, list of three values
    outfile:  prefix for output naming
    model_name: a label to slap on the output files; ex. "no_mig"
    func: access the model function from within 'dadi_Run_Optimizations.py' or from a separate python model script, ex. after importing Models_2D, calling Models_2D.no_mig
    rounds: number of optimization rounds to perform
    param_number: number of parameters in the model selected (can count in params line for the model)
    fs_folded: A Boolean value (True or False) indicating whether the empirical fs is folded (True) or not (False).

   Optional Arguments =
     reps: a list of integers controlling the number of replicates in each of the optimization rounds
     maxiters: a list of integers controlling the maxiter argument in each of the optimization rounds
     folds: a list of integers controlling the fold argument when perturbing input parameter values
     in_params: a list of parameter values
     in_upper: a list of upper bound values
     in_lower: a list of lower bound values
     param_labels: list of labels for parameters that will be written to the output file to keep track of their order
'''



'''
Now let's use the function to run an optimization routine for our data and this model.
This would perform three rounds of optimizations.
'''
#create a prefix to label the output files
prefix = "South_misid"
#make sure to define your extrapolation grid size
pts = [60,70,80]

reps = [20,30,50]
maxiters = [10,15,25]
folds = [3,2,1]

#Remember the order for mandatory arguments as below
#Optimize_Routine(fs, pts, outfile, model_name, func, rounds, param_number, fs_folded)

func_anc = dadi.Numerics.make_anc_state_misid_func(Demographics1D.two_epoch)
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "two_epoch", func_anc, 3, 3, fs_folded=False, reps = reps, maxiters = maxiters, folds = folds, param_labels = "nu, T, p_misid")

func_anc = dadi.Numerics.make_anc_state_misid_func(Demographics1D.growth)
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "growth", func_anc, 3, 3, fs_folded=False, reps = reps, maxiters = maxiters, folds = folds, param_labels = "nu, T, p_misid")

func_anc = dadi.Numerics.make_anc_state_misid_func(Demographics1D.three_epoch)
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "three_epoch", func_anc, 3, 5, fs_folded=False, reps = reps, maxiters = maxiters, folds = folds, param_labels = "nuB, nuF, TB, TF, p_misid")

upper = [50,30,30,30]
func_anc = dadi.Numerics.make_anc_state_misid_func(Demographics1D.bottlegrowth)
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "bottlegrowth", func_anc, 3, 4, fs_folded=False, reps = reps, maxiters = maxiters, folds = folds, param_labels = "nuB, nuF, T, p_misid", in_upper=upper)

