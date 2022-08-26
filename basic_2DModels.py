"""
Two population demographic models.
From dadi distribution, dadi_pipeline and my own additions.
"""
import numpy

from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum


"""
Basic models with no size changes
"""


def snm(notused, ns, pts):
    """
    ns = (n1,n2)

    Standard neutral model, populations never diverge.
    """
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def no_mig(params, ns, pts):
    """
    Split into two populations, no migration.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations)
    """
    nu1, nu2, T = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def sym_mig(params, ns, pts):
    """
    Split into two populations, with symmetric migration.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations)
    m: Migration rate between populations (2*Na*m)
    """
    nu1, nu2, m, T = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


"""
Models with size changes in the ancestral population
"""


def three_epoch(params, ns, pts):
    """
    params = (nuB,nuF,TB,TF)
    ns = (n1,n2)

    Two instantaneous size changes, populations never diverge.

    nuB: Ratio of population size after first instantanous change to ancient
         population size
    nuF: Ratio of population size after second instantanous change
    TB: Time in the past at which first instantaneous change happened
       (in units of 2*Na generations)
    TF: Time in the past at which second instantaneous change happened
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nuB,nuF,TB,TF = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)

    phi = Integration.one_pop(phi, xx, TB, nuB)
    phi = Integration.one_pop(phi, xx, TF, nuF)

    phi = PhiManip.phi_1D_to_2D(xx, phi)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def three_epoch_split(params, ns, pts):
    """
    params = (nuB,nuF,TB,TF,Ts,nu1,nu2)
    ns = (n1,n2)

    Two instantaneous size changes followed by split.

    nuB: Ratio of population size after first instantanous change to ancient
         population size
    nuF: Ratio of population size after second instantanous change
    TB: Time in the past at which first instantaneous change happened
       (in units of 2*Na generations)
    TF: Time in the past at which second instantaneous change happened
    Ts: Time in the past at which the two populations split.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nuB,nuF,TB,TF,Ts,nu1,nu2 = params
    return three_epoch_split_mig((nuB,nuF,0,TB,TF,Ts,nu1,nu2), ns, pts)


def three_epoch_split_mig(params, ns, pts):
    """
    params = (nuB,nuF,m,TB,TF,Ts,nu1,nu2)
    ns = (n1,n2)

    Two instantaneous size changes followed by split with
    migration.

    nuB: Ratio of population size after first instantanous change to ancient
         population size
    nuF: Ratio of population size after second instantanous change
    m: Migration rate between the two populations (2*Na*m).
    TB: Time in the past at which first instantaneous change happened
       (in units of 2*Na generations)
    TF: Time in the past at which second instantaneous change happened
    Ts: Time in the past at which the two populations split.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nuB,nuF,m,TB,TF,Ts,nu1,nu2 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)

    phi = Integration.one_pop(phi, xx, TB, nuB)
    phi = Integration.one_pop(phi, xx, TF, nuF)

    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def two_epoch(params, ns, pts):
    """
    params = (nu,T)
    ns = (n1,n2)

    Instantaneous size change, populations never diverge.

    nu: Ratio of population size after instantanous change to ancient
         population size
    T: Time in the past at which instantaneous change happened
       (in units of 2*Na generations)
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu,T = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)

    phi = Integration.one_pop(phi, xx, T, nu)

    phi = PhiManip.phi_1D_to_2D(xx, phi)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def two_epoch_split(params, ns, pts):
    """
    params = (nu,T,Ts,nu1,nu2)
    ns = (n1,n2)

    Two instantaneous size changes followed by split.

    nu: Ratio of population size after instantanous change to ancient
         population size
    T: Time in the past at which instantaneous change happened
       (in units of 2*Na generations)
    Ts: Time in the past at which the two populations split.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu,T,Ts,nu1,nu2 = params
    return two_epoch_split_mig((nu,0,T,Ts,nu1,nu2), ns, pts)


def two_epoch_split_mig(params, ns, pts):
    """
    params = (nu,m,T,Ts,nu1,nu2)
    ns = (n1,n2)

    Two instantaneous size changes followed by split with
    migration.

    nu: Ratio of population size after instantanous change to ancient
         population size
    T: Time in the past at which instantaneous change happened
       (in units of 2*Na generations)
    m: Migration rate between the two populations (2*Na*m).
    Ts: Time in the past at which the two populations split.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu,m,T,Ts,nu1,nu2 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)

    phi = Integration.one_pop(phi, xx, T, nu)

    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


"""
Models with size changes in both ancestral and split population
"""


def beforeafter_split_nomig(params, ns, pts):
    """
    params = (nu,T,Ts,nu1a,nu2a,Tb,nu1b,nu2b)
    ns = (n1,n2)

    Size change in the ancestral population, followed by split and another size
    change, absence of gene flow.

    nu: Ratio of population size after instantanous change to ancient
         population size
    T: Time in the past at which instantaneous change happened
       (in units of 2*Na generations)
    Ts: Time in the past at which the two populations split.
    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    Tb: The scaled time between the second size change and present.
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu,T,Ts,nu1a,nu2a,Tb,nu1b,nu2b = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)

    phi = Integration.one_pop(phi, xx, T, nu)

    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, Ts, nu1a, nu2a)
    phi = Integration.two_pops(phi, xx, Tb, nu1b, nu2b)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def beforeafter_split_secmig(params, ns, pts):
    """
    params = (nu,m,T,Ts,nu1a,nu2a,Tb,nu1b,nu2b)
    ns = (n1,n2)

    Size change in the ancestral population, followed by split with gene flow,
    then another size change and isolation.

    nu: Ratio of population size after instantanous change to ancient
         population size
    T: Time in the past at which instantaneous change happened
       (in units of 2*Na generations)
    m: Migration rate between the two populations (2*Na*m).
    Ts: Time in the past at which the two populations split.
    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    Tb: The scaled time between the second size change and present.
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu,m,T,Ts,nu1a,nu2a,Tb,nu1b,nu2b = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)

    phi = Integration.one_pop(phi, xx, T, nu)

    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, Ts, nu1a, nu2a, m12=m, m21=m)
    phi = Integration.two_pops(phi, xx, Tb, nu1b, nu2b, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def beforeafter_split_mig(params, ns, pts):
    """
    params = (nu,m,T,Ts,nu1a,nu2a,Tb,nu1b,nu2b)
    ns = (n1,n2)

    Size change in the ancestral population, followed by split with gene flow,
    then another size change and gene flow.

    nu: Ratio of population size after instantanous change to ancient
         population size
    T: Time in the past at which instantaneous change happened
       (in units of 2*Na generations)
    m: Migration rate between the two populations (2*Na*m).
    Ts: Time in the past at which the two populations split.
    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    Tb: The scaled time between the second size change and present.
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu,m,T,Ts,nu1a,nu2a,Tb,nu1b,nu2b = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)

    phi = Integration.one_pop(phi, xx, T, nu)

    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, Ts, nu1a, nu2a, m12=m, m21=m)
    phi = Integration.two_pops(phi, xx, Tb, nu1b, nu2b, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


"""
Models with size changes in split populations only
"""


def three_epoch_after_split_nomig(params, ns, pts):
    """
    Split followed by two size changes with no gene flow.
    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1: The scaled time between the split and the first size change (in units of 2*Na generations).
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    T2: The scale time between first and second size change.
    nu1c: Size of population 1 after time interval.
    nu2c: Size of population 2 after time interval.
    T3: The scaled time between the second size change and present.
    """
    nu1a, nu2a, nu1b, nu2b, nu1c, nu2c, T1, T2, T3 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1a, nu2a)

    phi = Integration.two_pops(phi, xx, T2, nu1b, nu2b)

    phi = Integration.two_pops(phi, xx, T3, nu1c, nu2c)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def three_epoch_after_split_secmig(params, ns, pts):
    """
    Split with gene flow, followed by size change with continuous symmetrical gene flow,
    then isolation after second size change.
    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1: The scaled time between the split and the first size change (in units of 2*Na generations).
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    T2: The scale time between first and second size change.
    nu1c: Size of population 1 after time interval.
    nu2c: Size of population 2 after time interval.
    T3: The scaled time between the second size change and present.
    m: Migration between pop 2 and pop 1.
    """
    nu1a, nu2a, nu1b, nu2b, nu1c, nu2c, m, T1, T2, T3 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1a, nu2a, m12=m, m21=m)

    phi = Integration.two_pops(phi, xx, T2, nu1b, nu2b, m12=m, m21=m)

    phi = Integration.two_pops(phi, xx, T3, nu1c, nu2c, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def three_epoch_after_split_mig(params, ns, pts):
    """
    Split with gene flow, followed by two size changes with continuous symmetrical gene flow.
    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    T1: The scaled time between the split and the first size change (in units of 2*Na generations).
    nu1b: Size of population 1 after time interval.
    nu2b: Size of population 2 after time interval.
    T2: The scale time between first and second size change.
    nu1c: Size of population 1 after time interval.
    nu2c: Size of population 2 after time interval.
    T3: The scaled time between the second size change and present.
    m: Migration between pop 2 and pop 1.
    """
    nu1a, nu2a, nu1b, nu2b, nu1c, nu2c, m, T1, T2, T3 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1a, nu2a, m12=m, m21=m)

    phi = Integration.two_pops(phi, xx, T2, nu1b, nu2b, m12=m, m21=m)

    phi = Integration.two_pops(phi, xx, T3, nu1c, nu2c, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
