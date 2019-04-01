import sys
sys.path.append('/Users/denisemartini/dadi_pipeline/Plotting')
sys.path.append('/Users/denisemartini/dadi')

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

# Now doing the same for the SI pop and I am afraid I will see the same problem.

fs = South_fs
pts = [40,50,60]
prefix = "South"
emp_params = [0.1314,1.6206,0.1034,0.4908]

model_fit = Plotting_Functions.Fit_Empirical(fs, pts, prefix, "three_epoch", Demographics1D.three_epoch, emp_params, fs_folded=False)

Plotting_Functions.Plot_1D(fs, model_fit, prefix, "SI_three_epoch")

# I will try to fix this with the ancestral misidentification correction
# that comes with the latest version of dadi. Basically, you add an extra
# parameter to the model, that accounts for ancestral misidentification.
# If this fixes my models great, otherwise I will have to fold the spectrum.

# Let's try some plotting again, after adding the correction for ancestral misidentification.

fs = South_fs
pts = [60,70,80]
prefix = "South_misid"
emp_params = [47.5337,1.6926,0.6433,0.0433]

func_anc = dadi.Numerics.make_anc_state_misid_func(Demographics1D.bottlegrowth)

model_fit = Plotting_Functions.Fit_Empirical(fs, pts, prefix, "bottlegrowth", func_anc, emp_params, fs_folded=False)

Plotting_Functions.Plot_1D(fs, model_fit, prefix, "SI_bottlegrowth")

# Much better but still not quite a proper fit. It might have to do with the fact that
# the nuB parameter gets so close to the upper boundaries once again.

# Expecting the same for the NI:
fs = North_fs
pts = [60,70,80]
prefix = "North_misid"
emp_params = [49.6807,1.3588,0.7463,0.0435]

func_anc = dadi.Numerics.make_anc_state_misid_func(Demographics1D.bottlegrowth)

model_fit = Plotting_Functions.Fit_Empirical(fs, pts, prefix, "bottlegrowth", func_anc, emp_params, fs_folded=False)

Plotting_Functions.Plot_1D(fs, model_fit, prefix, "NI_bottlegrowth")

# Same here. I do hope it is because of the boundaries. I will try fitting again
# with increased upper boundaries. If that keeps failing I might have to
# give up and fold the spectrum.

# Trying the other models best fit, out of curiosity.
emp_params = [2.9529,1.1872,0.0449]
func_anc = dadi.Numerics.make_anc_state_misid_func(Demographics1D.two_epoch)
model_fit = Plotting_Functions.Fit_Empirical(fs, pts, prefix, "two_epoch", func_anc, emp_params, fs_folded=False)
Plotting_Functions.Plot_1D(fs, model_fit, prefix, "NI_two_epoch")

emp_params = [5.6064,6.7774,0.0445]
func_anc = dadi.Numerics.make_anc_state_misid_func(Demographics1D.growth)
model_fit = Plotting_Functions.Fit_Empirical(fs, pts, prefix, "growth", func_anc, emp_params, fs_folded=False)
Plotting_Functions.Plot_1D(fs, model_fit, prefix, "NI_growth")

emp_params = [0.4734,2.2299,0.4451,0.9992,0.0439]
func_anc = dadi.Numerics.make_anc_state_misid_func(Demographics1D.three_epoch)
model_fit = Plotting_Functions.Fit_Empirical(fs, pts, prefix, "three_epoch", func_anc, emp_params, fs_folded=False)
Plotting_Functions.Plot_1D(fs, model_fit, prefix, "NI_three_epoch")

# the bottlegrowth is definitely the best fit...
# After another round...

emp_params = [12.2993,1.1869,0.7961,0.1057,0.0494]
func_anc = dadi.Numerics.make_anc_state_misid_func(Demographics1D.three_epoch)
model_fit = Plotting_Functions.Fit_Empirical(fs, pts, prefix, "three_epoch", func_anc, emp_params, fs_folded=False)
Plotting_Functions.Plot_1D(fs, model_fit, prefix, "NI_three_epoch")

emp_params = [6.9479,1.3472,0.7791,0.0455]
func_anc = dadi.Numerics.make_anc_state_misid_func(Demographics1D.bottlegrowth)
model_fit = Plotting_Functions.Fit_Empirical(fs, pts, prefix, "bottlegrowth", func_anc, emp_params, fs_folded=False)
Plotting_Functions.Plot_1D(fs, model_fit, prefix, "NI_bottlegrowth")

fs = South_fs
pts = [60,70,80]
prefix = "South_misid"
emp_params = [35.2712,1.9514,0.8587,0.0528]
func_anc = dadi.Numerics.make_anc_state_misid_func(Demographics1D.bottlegrowth)
model_fit = Plotting_Functions.Fit_Empirical(fs, pts, prefix, "bottlegrowth", func_anc, emp_params, fs_folded=False)
Plotting_Functions.Plot_1D(fs, model_fit, prefix, "SI_bottlegrowth")

# I really don't see it approaching a good fit. Time to fold the FS.

# What actually happens if I fold the FS?
fs = dadi.Spectrum.from_data_dict(dd, ['North','South'], [38,30], polarized = False)

dadi.Plotting.plot_single_2d_sfs(fs, vmin=0.1)

# Still a lot of rare alleles...and it gets a bit asymmetrical it seems.
# quickly saving this figure:
fig = pylab.figure()
dadi.Plotting.plot_single_2d_sfs(fs, vmin=0.1)
fig.savefig('../pop_structure/dadi/dadi2Dfs_folded.png', dpi=300, bbox_inches='tight')

fs = dadi.Spectrum.from_data_dict(dd, ['North'], [38], polarized = False)
dadi.Plotting.plot_1d_fs(fs)

fs = dadi.Spectrum.from_data_dict(dd, ['South'], [30], polarized = False)
dadi.Plotting.plot_1d_fs(fs)

# So, now that I ran the models with the folded spectrum, I definitely
# get better likelihoods, we are still between three epoch and bottlegrowth,
# in both cases with a population explosion first followed by a decrease.
# But will they actually fit?

fs = dadi.Spectrum.from_data_dict(dd, ['North'], [38], polarized = False)
pts = [60,70,80]
prefix = "North_folded"
emp_params = [29.9482,1.7453,0.9815]
model_fit = Plotting_Functions.Fit_Empirical(fs, pts, prefix, "bottlegrowth", Demographics1D.bottlegrowth, emp_params, fs_folded=True)
Plotting_Functions.Plot_1D(fs, model_fit, prefix, "NI_bottlegrowth")

emp_params = [21.8656,0.5254,0.7354,0.0467]
model_fit = Plotting_Functions.Fit_Empirical(fs, pts, prefix, "three_epoch", Demographics1D.three_epoch, emp_params, fs_folded=True)
Plotting_Functions.Plot_1D(fs, model_fit, prefix, "NI_three_epoch")

# Looks like the three epoch is a better fit, the residuals have a smaller range,
# but it might not be completely optimized yet.
# Let's take a look at the SI. Here the likelihoods are much closer.

fs = dadi.Spectrum.from_data_dict(dd, ['South'], [30], polarized = False)
pts = [60,70,80]
prefix = "South_folded"
emp_params = [29.9944,1.9569,0.8404]
model_fit = Plotting_Functions.Fit_Empirical(fs, pts, prefix, "bottlegrowth", Demographics1D.bottlegrowth, emp_params, fs_folded=True)
Plotting_Functions.Plot_1D(fs, model_fit, prefix, "SI_bottlegrowth")

emp_params = [9.3137,0.6909,0.6636,0.0386]
model_fit = Plotting_Functions.Fit_Empirical(fs, pts, prefix, "three_epoch", Demographics1D.three_epoch, emp_params, fs_folded=True)
Plotting_Functions.Plot_1D(fs, model_fit, prefix, "SI_three_epoch")
# Probably the best fit so far, ever...still isn't perfect but it would be
# enough to go on with I think. What is bugging me is using the folded FS really,
# because the signal I am getting is exactly the same and I really am
# wondering if what I need isn't simply more replicates to find a really
# good fit there as well. We definitely were close to these parameters
# there as well. I will try one last time, giving it very many iterations.
