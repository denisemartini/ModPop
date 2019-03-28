import dadi

dd = dadi.Misc.make_data_dict('../pop_structure/dadi/')

fs = dadi.Spectrum.from_data_dict(dd, ['YRI','CEU'], [20,20])

import pylab

dadi.Plotting.plot_single_2d_sfs(fs, vmin=0.1)

from numpy import array
import sys
sys.path.append('/Users/denisemartini/gutenkunstlab-dadi-f2f4b565089a/examples/YRI_CEU')

import demographic_models

data = dadi.Spectrum.from_file('/Users/denisemartini/gutenkunstlab-dadi-f2f4b565089a/examples/YRI_CEU/YRI_CEU.fs')
ns = data.sample_sizes
pts_l = [40,50,60]

func = dadi.Demographics2D.split_mig
func = demographic_models.prior_onegrow_mig

# Now let's optimize parameters for this model.

# The upper_bound and lower_bound lists are for use in optimization.
# Occasionally the optimizer will try wacky parameter values. We in particular
# want to exclude values with very long times, very small population sizes, or
# very high migration rates, as they will take a long time to evaluate.
# Parameters are: (nu1F, nu2B, nu2F, m, Tp, T)
upper_bound = [100, 100, 100, 10, 3, 3]
lower_bound = [1e-2, 1e-2, 1e-2, 0, 0, 0]

# This is our initial guess for the parameters, which is somewhat arbitrary.
p0 = [2,0.1,2,1,0.2,0.2]
# Make the extrapolating version of our demographic model function.
func_ex = dadi.Numerics.make_extrap_log_func(func)

# Perturb our parameters before optimization. This does so by taking each
# parameter a up to a factor of two up or down.
p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,
                              lower_bound=lower_bound)
# Do the optimization. By default we assume that theta is a free parameter,
# since it's trivial to find given the other parameters. If you want to fix
# theta, add a multinom=False to the call.
# The maxiter argument restricts how long the optimizer will run. For real
# runs, you will want to set this value higher (at least 10), to encourage
# better convergence. You will also want to run optimization several times
# using multiple sets of intial parameters, to be confident you've actually
# found the true maximum likelihood parameters.
print('Beginning optimization ************************************************')
popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l,
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(p0), maxiter=3)
# The verbose argument controls how often progress of the optimizer should be
# printed. It's useful to keep track of optimization process.
print('Finshed optimization **************************************************')

# These are the actual best-fit model parameters, which we found through
# longer optimizations and confirmed by running multiple optimizations.
# We'll work with them through the rest of this script.
popt = [1.881, 0.0710, 1.845, 0.911, 0.355, 0.111]
print('Best-fit parameters: {0}'.format(popt))

# Calculate the best-fit model AFS.
model = func_ex(popt, ns, pts_l)
# Likelihood of the data given the model AFS.
ll_model = dadi.Inference.ll_multinom(model, data)
print('Maximum log composite likelihood: {0}'.format(ll_model))
# The optimal value of theta given the model.
theta = dadi.Inference.optimal_sfs_scaling(model, data)
print('Optimal value of theta: {0}'.format(theta))

# Plot a comparison of the resulting fs with the data.
import pylab
pylab.figure(1)
dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range=3,
                                    pop_ids =('YRI','CEU'))

# Save the figure
pylab.savefig('YRI_CEU.png', dpi=50)

# Estimate parameter uncertainties using the Godambe Information Matrix, to
# account for linkage in the data. To use the GIM approach, we need to have
# spectra from bootstrapping our data.  Let's load the ones we've provided for
# the example.
# (We're using Python list comprehension syntax to do this in one line.)
all_boot = [dadi.Spectrum.from_file('/Users/denisemartini/gutenkunstlab-dadi-f2f4b565089a/examples/YRI_CEU/bootstraps/{0:02d}.fs'.format(ii))
            for ii in range(100)]
uncerts = dadi.Godambe.GIM_uncert(func_ex, pts_l, all_boot, popt, data,
                                  multinom=True)
# uncert contains the estimated standard deviations of each parameter, with
# theta as the final entry in the list.
print('Estimated parameter standard deviations from GIM: {0}'.format(uncerts))

# For comparison, we can estimate uncertainties with the Fisher Information
# Matrix, which doesn't account for linkage in the data and thus underestimates
# uncertainty. (Although it's a fine approach if you think your data is truly
# unlinked.)
uncerts_fim = dadi.Godambe.FIM_uncert(func_ex, pts_l, popt, data, multinom=True)
print('Estimated parameter standard deviations from FIM: {0}'.format(uncerts_fim))

print('Factors by which FIM underestimates parameter uncertainties: {0}'.format(uncerts/uncerts_fim))

# What if we fold the data?
# These are the optimal parameters when the spectrum is folded. They can be
# found simply by passing data.fold() to the above call to optimize_log.
popt_fold =  array([1.907,  0.073,  1.830,  0.899,  0.425,  0.113])
uncerts_folded = dadi.Godambe.GIM_uncert(func_ex, pts_l, all_boot, popt_fold,
                                         data.fold(), multinom=True)
print('Folding increases parameter uncertainties by factors of: {0}'.format(uncerts_folded/uncerts))

# Let's do a likelihood-ratio test comparing models with and without migration.
# The no migration model is implemented as
# demographic_models.prior_onegrow_nomig
func_nomig = demographic_models.prior_onegrow_nomig
func_ex_nomig = dadi.Numerics.make_extrap_log_func(func_nomig)
# These are the best-fit parameters, which we found by multiple optimizations
popt_nomig = array([ 1.897,  0.0388,  9.677,  0.395,  0.070])
model_nomig = func_ex_nomig(popt_nomig, ns, pts_l)
ll_nomig = dadi.Inference.ll_multinom(model_nomig, data)

# Since LRT evaluates the complex model using the best-fit parameters from the
# simple model, we need to create list of parameters for the complex model
# using the simple (no-mig) best-fit params.  Since evalution is done with more
# complex model, need to insert zero migration value at corresponding migration
# parameter index in complex model. And we need to tell the LRT adjust function
# that the 3rd parameter (counting from 0) is the nested one.
p_lrt = [1.897,  0.0388,  9.677, 0, 0.395,  0.070]

adj = dadi.Godambe.LRT_adjust(func_ex, pts_l, all_boot, p_lrt, data,
                              nested_indices=[3], multinom=True)
D_adj = adj*2*(ll_model - ll_nomig)
print('Adjusted D statistic: {0:.4f}'.format(D_adj))

# Because this is test of a parameter on the boundary of parameter space
# (m cannot be less than zero), our null distribution is an even proportion
# of chi^2 distributions with 0 and 1 d.o.f. To evaluate the p-value, we use the
# point percent function for a weighted sum of chi^2 dists.
pval = dadi.Godambe.sum_chi2_ppf(D_adj, weights=(0.5,0.5))
print('p-value for rejecting no-migration model: {0:.4f}'.format(pval))


## These were the examples from the dadi distribution, now let's see what happens with my data.
import dadi
import pylab
from numpy import array

dd = dadi.Misc.make_data_dict('../pop_structure/dadi/final_snps_for_dadi.recode.vcf.data')
# I expected these two commands to run very slowly, but they were actually quite fast
# this is just trying out with the standard projection commands
fs = dadi.Spectrum.from_data_dict(dd, ['North','South'], [20,20], polarized = True)
# let's see how many segregating sites I have left like that:
S = fs.S()
round(S)
# 70432, that sounds like a lot, but it is almost 20,000 snps thrown, so we can probably do better.
# I know that the missingness in my dataset should be at ~30%. So let's try projecting down to 70% of the alleles.
fs = dadi.Spectrum.from_data_dict(dd, ['North','South'], [68,60], polarized = True)
S = fs.S()
round(S)
# 56242, so that is a lot worse...maybe we can find a good spot in between?
fs = dadi.Spectrum.from_data_dict(dd, ['North','South'], [38,34], polarized = True)
S = fs.S()
round(S)
# 71542 (38,34)
# 72100 (36,32)
# 72552 (34,30)
# trying all combinations is a bit silly, I should check what maximises each population in turn:
seg_sites = {}
for x in range(24, 52):
            fs =  dadi.Spectrum.from_data_dict(dd, ['North'], [x], polarized=True)
            s = fs.S()
            seg_sites[x] = round(s)
            print(x, seg_sites[x])

# highest value for the North pop is definitely 38, even though 40 and 36 are close.
seg_sites = {}
for x in range(20, 46):
            fs =  dadi.Spectrum.from_data_dict(dd, ['South'], [x], polarized=True)
            s = fs.S()
            seg_sites[x] = round(s)
            print(x, seg_sites[x])
# and 30/32 is the best value for the South pop, so it looks like I was close. Anyway, let's put it together.
fs = dadi.Spectrum.from_data_dict(dd, ['North','South'], [38,30], polarized = True)
S = fs.S()
round(S)
# 72831 is the end number and it sounds good.

# Let's take a look at the spectrum.
dadi.Plotting.plot_single_2d_sfs(fs, vmin=0.1)

# This image is telling me that there are loads of alleles that are at low frequency in both populations
# and really not many that are fixed in one or the other (low frequency in one pop and high in the other)
# but for a small patch of near fixation in the South pop. But those might be really few snps.
# I want to save this fs to file.
fs.to_file('../pop_structure/dadi/dadi2D.fs')

# Save the figure
fig = pylab.figure()
dadi.Plotting.plot_single_2d_sfs(fs, vmin=0.1)
fig.savefig('../pop_structure/dadi/dadi2Dfs.png', dpi=300, bbox_inches='tight')
# note for self: you need to run the above three lines all together, or hydrogen
# evaluates them separately and does not actually print out the figure
