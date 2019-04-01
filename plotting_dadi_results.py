import dadi
import pylab
from numpy import array

# first let's get my real data in here
dd = dadi.Misc.make_data_dict('../pop_structure/dadi/final_snps_for_dadi.recode.vcf.data')

North_fs = dadi.Spectrum.from_data_dict(dd, ['North'], [38], polarized = True)
South_fs = dadi.Spectrum.from_data_dict(dd, ['South'], [30], polarized = True)

# import the models I used
from dadi import Demographics1D

# then I can use the functions from the dadi pipeline to fit the models with the optimised parameters
import sys
sys.path.append('/Users/denisemartini/dadi_pipeline/Plotting')
import Plotting_Functions

# plotting for the NI population, bottlegrowth model:
fs = North_fs
pts = [40,50,60]
prefix = "North"
emp_params = [29.9791,1.0686,0.4496]

model_fit = Plotting_Functions.Fit_Empirical(fs, pts, prefix, "bottlegrowth", Demographics1D.bottlegrowth, emp_params, fs_folded=False)

Plotting_Functions.Plot_1D(fs, model_fit, prefix, "NI_bottlegrowth")

# it is not a great fit for the high frequency shared alleles.
# I am afraid there might be a problem with the unfolded spectra
# (ancestral misidentification error), I might have to fold the FS.
# Because to apply a correction to that error I would need the flanking
# bases of both ref and outgroup and that would be messy I think.
