#!/usr/bin/env python3
#
# Pre-computes the PT tables used to model P(k).
# Right now this uses Zeldovich.
#
import numpy as np
import astropy.cosmology
import sys

import zeldovich as Z




def growth_factor(cc,zz):
    """Approximate the growth factor, normalized to 1 today."""
    afid = 1.0/(1.0+zz)
    zval = 1/np.logspace(np.log10(afid),0.0,100)-1.0
    Dz   = np.exp(-np.trapz(cc.Om(zval)**0.55,x=np.log(1/(1+zval))))
    return(Dz)




def save_tables(cosmo,zlist,pkfn,zfid):
    """Computes tables for each z in zlist.
       'cosmo' is an astropy Cosmology instance,
       'zlist' is a list of redshifts at which to make the tabes,
       'pkfn' is a file to read P(k) from (k and P, in Mpc/h units),
       'zfid' is the fiducial redshift at which P(k) is evaluated.
       Saves the result as a series of (named) text files."""
    # Read P(k) and set up some basic cosmology values.
    pk   = np.loadtxt(pkfn)
    hub  = cosmo.H0.value / 100.
    Dfid = growth_factor(cosmo,zfid)
    # Now loop through each redshift and make a table for that redshift
    # saving it to a named file.  This can be used later to quickly set
    # up a Zeldovich (or PT) instance.
    for z in zlist:
        iz  = int(np.rint(100*z))
        Dz  = growth_factor(cosmo,z)/Dfid
        zel = Z.Zeldovich(pk[:,0],pk[:,1]*Dz**2)
        zel.make_table(kmin=0.01,kmax=1.0,nk=100)
        zel.write_table("PT_table_z{:03d}.txt".format(iz))
    #






if __name__=="__main__":
    if len(sys.argv)!=3:
        raise RuntimeError("Usage: {:s} ".format(sys.argv[0])+\
                           "<Pk-fname> <zfid>")
    else:
        zlist = [0.25,0.5,1.0,1.5,2.0]
        pkfn  = sys.argv[1]
        zfid  = float(sys.argv[2])
        # Assume these is ascii text file with two columns: k,P(k),
        # being the linear power spectrum at z.
        save_tables(astropy.cosmology.Planck15,zlist,pkfn,zfid)
    #
