import sys
import os
import numpy
import dadi
import pylab
from datetime import datetime
import Optimize_Functions
import basic_2DModels

'''
Usage: python dadi_Run_2D_Set.py

This is a modified version of the 'dadi_Run_Optimizations.py' script in which
we run optimizations for 2D comparisons for a large set of models that have been
made available as part of published works. These models are stored in the
Models_2D.py script, and will be called directly here. The user can delete or
comment out models to analyze a subset of the models available.
This script must be in the same working directory as Optimize_Functions.py, which
contains all the functions necessary, as well as the  Models_2D.py script, which
has all the model definitions.

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

Citations:
 If you use these scripts or the main diversification models for your work, please
 cite the following publication:
    Portik, D.M., Leache, A.D., Rivera, D., Blackburn, D.C., Rodel, M.-O.,
    Barej, M.F., Hirschfeld, M., Burger, M., and M.K. Fujita. 1017.
    Evaluating mechanisms of diversification in a Guineo-Congolian forest
    frog using demographic model selection. Molecular Ecology 26: 5245-5263.
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
Updated May 1018

Modified from the dadi_pipeline original version.
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
pop_ids=["North", "South"]

#**************
#projection sizes, in ALLELES not individuals
proj = [38,30]

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
# Calling external 2D models from the Models_2D.py script
#================================================================================
'''
 We will use a function from the Optimize_Functions.py script for our optimization routines:

 Optimize_Routine(fs, pts, outfile, model_name, func, rounds, param_number, fs_folded=True, reps=None, maxiters=None, folds=None, in_params=None, in_upper=None, in_lower=None, param_labels=" ")

   Mandatory Arguments =
    fs:  spectrum object name
    pts: grid size for extrapolation, list of three values
    outfile:  prefix for output naming
    model_name: a label to help label the output files; ex. "no_mig"
    func: access the model function from within 'moments_Run_Optimizations.py' or from a separate python model script, ex. after importing Models_2D, calling Models_2D.no_mig
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

#create a prefix based on the population names to label the output files
#ex. Pop1_Pop2
prefix = "2Dfs"

#**************
#make sure to define your extrapolation grid size (based on your projections)
pts = [80,90,100]

#**************
#Set the number of rounds here
rounds = 4

#define the lists for optional arguments
#you can change these to alter the settings of the optimization routine
reps = [100,80,60,50]
maxiters = [100,60,40,30]
folds = [3,2,1,1]

#**************
#Indicate whether your frequency spectrum object is folded (True) or unfolded (False)
fs_folded = False


"""
Models with size changes in split populations only
"""

func_anc = dadi.Numerics.make_anc_state_misid_func(basic_2DModels.three_epoch_after_split_nomig)
params = [1.34,0.23,1.47,7.57,0.67,0.85,0.001,0.009,0.01,0.009]
lower = [0.01,0.01,0.01,0.01,0.01,0.01,0.0001,0.0001,0.001,0.0001]
upper = [30,30,30,30,30,30,10,10,10,1]
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "three_epoch_after_split_nomig", func_anc, rounds, 10, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a,nu2a,nu1b,nu2b,nu1c,nu2c,T1,T2,T3,p_misid", in_params=params, in_upper=upper, in_lower=lower)


func_anc = dadi.Numerics.make_anc_state_misid_func(basic_2DModels.three_epoch_after_split_secmig)
params = [8.55,0.25,1.89,0.99,0.086,15.85,10.84,0.20,0.29,0.002,0.006]
lower = [0.001,0.01,0.01,0.01,0.001,0.01,0.01,0.01,0.001,0.0001,0.0001]
upper = [30,30,30,30,30,30,20,10,10,10,1]
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "three_epoch_after_split_secmig", func_anc, rounds, 11, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a,nu2a,nu1b,nu2b,nu1c,nu2c,m,T1,T2,T3,p_misid", in_params=params, in_upper=upper, in_lower=lower)


func_anc = dadi.Numerics.make_anc_state_misid_func(basic_2DModels.three_epoch_after_split_mig)
params = [0.45,0.58,1.54,4.38,0.93,1.37,8.17,0.99,1.03,0.07,0.03]
lower = [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.001,0.001]
upper = [30,30,30,30,30,30,20,10,10,10,1]
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "three_epoch_after_split_mig", func_anc, rounds, 11, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a,nu2a,nu1b,nu2b,nu1c,nu2c,m,T1,T2,T3,p_misid", in_params=params, in_upper=upper, in_lower=lower)
