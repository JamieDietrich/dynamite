# Grab all stars and known parameters, including planets, from ExEP Precursor Science Target List
## Get stellar values uniformly from Caleb Harada's SPORES work

# For stars with no known planets, create simulated systems via Monte Carlo iterations
## First draw from giant planet occurrence rate
## For each draw, create Kepler-like small planet system out to stability limit with giant planet
## Test both 0% and 80% truncation at 100-300 days
## Find likelihood of <1.7 R_Earth planet (or <5 M_Earth planet?) within stellar HZ
## Find likelihood of giant planet within stellar + tides HZ

# For stars with known planets, skip the above

# For all systems
## Test stability and likelihood of an additional rocky temperate planet in each MC iteration - likely very low
## Test systems with giant planets for different values of likelihood upper limit for Earth-like moon at stellar + tides HZ
## Combine with CHZ metric from Austin

import math
import traceback
import numpy as np
import scipy.stats as spst
import matplotlib.pyplot as plt
import scipy.interpolate as spint

from dynamite_test import *
from datetime import datetime, timedelta

class hwo_sims:

    def __init__(self):
        """Runs the simulations."""

        try:
            # Get stellar parameters from new database
            self.HDname = np.loadtxt("spores_catalog_v2.1.0.csv", delimiter=',', dtype=str, skiprows=162, usecols=(5))
            self.steff, self.srad, self.smass, self.sage = np.loadtxt("spores_catalog_v2.1.0.csv", delimiter=',', dtype=float, unpack=True, skiprows=162, usecols=(26, 32, 34, 221))                
            self.dyn = dynamite(False)

            # Create/use planetary systems and test for additional planets and CHZ value

            ## Use giant planet occurrence rate power law from Bergsten et al. (in prep)
            self.gpor_params = [0.963, 0.25, 2.7, 1, 0.1, 0.05, -0.665, 0.477, 0.81, 1.25]

            ## Set up period ratio grid for small planets
            self.PRgrid = np.logspace(0,1)

            with np.errstate(divide='ignore'):
                Dgrid = np.log(2.*(self.PRgrid**(2./3.)-1.)/(self.PRgrid**(2./3.)+1.))

            Dgrid[0] = -4
            self.cdfPR = spst.norm(-0.9, 0.41).cdf(Dgrid)

            ## Multiplicity-dependent inclination and eccentricity distributions from He et al. (2020)
            self.inc_sigmas = [0.84, 0.85, 0.86, 0.86, 0.85, 0.84, 0.81, 0.79, 0.77]
            self.ecc_sigmas = [0.641, 0.670, 0.689, 0.701, 0.704, 0.685, 0.666, 0.632, 0.612, 0.587]

            ## Stellar Habitable Zone Limit Interpolation
            """
            turbet_grid = [[2300 + i * 100 for i in range(41)], [0.5, 1, 4.5, 10],
                           [0.772, 0.774, 0.777, 0.779, 0.782, 0.784, 0.787, 0.790, 0.793, 0.796, 0.799, 0.803, 0.806, 0.810, 0.814, 0.817, 0.820, 0.824, 0.828, 0.832, 0.838, 0.843, 0.850, 0.858, 0.866, 0.875, 0.884, 0.894, 0.903, 0.912, 0.922, 0.932, 0.943, 0.962, 0.988, 1.014, 1.040, 1.071, 1.101, 1.130, 1.156],
                           [0.772, 0.774, 0.777, 0.779, 0.782, 0.784, 0.787, 0.790, 0.793, 0.796, 0.800, 0.806, 0.831, 0.830, 0.824, 0.827, 0.831, 0.836, 0.843, 0.850, 0.859, 0.867, 0.875, 0.885, 0.895, 0.905, 0.916, 0.929, 0.941, 0.953, 0.971, 0.992, 1.014, 1.036, 1.066, 1.096, 1.126, 1.164, 1.203, 1.242, 1.281],
                           [0.772, 0.774, 0.778, 0.781, 0.783, 0.784, 0.789, 0.792, 0.795, 0.800, 0.816, 0.845, 0.858, 0.854, 0.849, 0.856, 0.863, 0.872, 0.884, 0.896, 0.910, 0.923, 0.939, 0.957, 0.976, 0.995, 1.016, 1.046, 1.076, 1.105, 1.134, 1.186, 1.241, 1.297, 1.354, 1.455, 1.565, 1.709, 2, 2, 2],
                           [0.775, 0.778, 0.781, 0.783, 0.786, 0.788, 0.791, 0.795, 0.800, 0.808, 0.838, 0.881, 0.880, 0.876, 0.875, 0.890, 0.905, 0.921, 0.942, 0.964, 0.988, 1.011, 1.037, 1.075, 1.114, 1.155, 1.196, 1.240, 1.325, 1.412, 1.499, 1.588, 2, 2, 2, 2, 2, 2, 2, 2, 2]]

            self.interp = spint.RegularGridInterpolator((turbet_grid[1], turbet_grid[0]), [turbet_grid[t] for t in range(2, 6)],)
            """
            ## Run simulations with Monte Carlo iterations for each stellar system
            self.mc_iter = 100000
            res = self.run_new_mp(self.run_sims, np.arange(len(self.srad)), ())
            results = []

            for j in range(len(res[0])):
                res1 = []
                
                for i in range(len(res)):
                    res1.append(res[i][j])

                results.append(np.array(res1))

            oes, ces = results

            with open("values.txt", 'w') as f:
                f.write("System\tOptimistic\tConservative\n")
                
                for i in range(len(oes)):
                    f.write(self.HDname[i] + "\t" + str(oes[i]/self.mc_iter) + "\t" + str(ces[i]/self.mc_iter) + "\n")

            plt.hist(oes/self.mc_iter)
            plt.savefig("optimistic_eta_hist.png")
            plt.clf()
            plt.hist(ces/self.mc_iter)
            plt.savefig("conservative_eta_hist.png")
            plt.close()
            """
            print("The stars in this list have an optimistic estimated eta_Earth of", round((np.mean(oes))/self.mc_iter, 4), "+/-", round(np.std(oes)/self.mc_iter, 4), "and a conservative estimated eta_Earth of", round(np.mean(ces)/self.mc_iter, 4), "+/-", round(np.std(ces)/self.mc_iter, 4))
            print("The star with the highest estimated optimistic eta_Earth of", round(max(oes)/self.mc_iter, 4),"is", self.HDname[np.where(oes == max(oes))[0][0]], "and the star with the highest estimated conservative eta_Earth of", round(max(ces)/self.mc_iter, 4), "is", self.HDname[np.where(ces == max(ces))[0][0]])
            exit()
            exit()
            """
        except Exception as e:
            print(e)
            exit()

        exit()

    def cdf_draw(self, grid, val, cdf):
        """Draws from the cdf of a distribution"""

        if len(np.where(val - cdf < 0)[0]) == 0:
            return grid[-1]

        else:
            return grid[np.where(val - cdf < 0)[0][0]]

    def K3(self, val, meas, M):
        """Calculates semi-major axis in au using period in days or vice versa (with mass in solar masses)"""

        if meas == "a":
            return (4*math.pi**2*(val*1.5e11)**3/(6.67e-11*M*2e30))**(1/2)/86400

        else:
            return (6.67e-11*M*2e30*(val*86400)**2/(4*math.pi**2))**(1/3)/1.5e11

    def otegi_mr(self, measurement, predict):
        """Uses a power-law MR to predict mass in Earth values from radius in Earth values or vice versa"""

        if predict == "radius":
            if measurement > 5 and measurement < 120:
                return 0.7*measurement**0.63

            elif measurement > 120:
                return 14.3*(measurement/120)**-0.044

            return 1.03*measurement**0.29

        elif predict == "mass":
            M_r = 0.9*measurement**3.45
            M_v = 1.74*measurement**1.58
            
            if M_r > 5:
                return M_v

            return M_r

    def dynamite_prediction(self, args, P, KP, i):
        """Runs the necessary functions in DYNAMITE to get a prediction for one MC instance"""

        name, per, rad, mas, inc, ecc, R_star, M_star = args
        stable = "NO"
        GMfp213 = (self.dyn.G*M_star*self.dyn.M_sun/(4*math.pi**2))**(1/3)
        em = np.median(ecc)
        il, Pinc, cdfi, inew, ib = self.dyn.syssim_incs(inc, per, GMfp213, R_star, False)
        el, Pecc, cdfe = self.dyn.syssim_eccs(len(ecc))
            
        if not KP:
            R, PR, cdfR = self.dyn.syssim_rads(rad)
            PP, deltas, cdfP = self.dyn.create_pers(self.dyn.create_epos_fP(P, per), per, mas, [], ecc, P, R, cdfR, il, Pinc, Pecc, em, cdfe, M_star, GMfp213)

        else:
            pass

        if sum(cdfP) == 0:
            return (0, 0, 0)
        
        try:
            scnt = 0
            
            while stable != "YES" and scnt < 1000:
                res = self.dyn.mc_test(0, (P, PP, cdfP, per, deltas, R, PR, cdfR, inew, ib, il, cdfi, el, Pecc, em, cdfe, [(p, 0, 0) for p in per], [(r, 0, 0) for r in rad], [(m, 0, 0, "Mass") for m in mas], [(i, 0, 0) for i in inc], [(e, 0, 0) for e in ecc], GMfp213, R_star, M_star, [], False))
                stable = res[4]
                scnt += 1

            if scnt == 1000:
                    print("COULDN'T FIND STABILITY FOR ITERATION", i, "OF SYSTEM", name, "AFTER 1000 TRIES, RETURNING NO NEW PLANET", datetime.now())
                    return (0, 0, 0)

        except:
            print(name, "HAD DYNAMITE ERROR FOR ITERATION " + str(i) + ", RETURNING NO NEW PLANET,", datetime.now())
            traceback.print_exc()
            return (0, 0, 0)

        return (res[7], res[8], self.otegi_mr(res[8], "mass"))

    def run_new_mp(self, func, arr, mp_args, domp=True):
        """Runs new multiprocessing code needed for iOS users."""

        if domp:
            with mp.Pool(processes=mp.cpu_count()) as pool:
                args = [(i, mp_args) for i in arr]
                results = pool.starmap(func, args)

        else:
            results = []
            
            for c, i in enumerate(arr):
                results.append(func(i, mp_args))

        return results

    def run_sims(self, i, args):
        """Runs the simulated systems with multiprocessing."""

        start = datetime.now()
        print("Starting", self.HDname[i], str(i+1) + "/" + str(len(self.srad)), start)
        
        # calculate HZ radii for star
        krungh_iHZ = (1.0385 + 1.2456e-4*(self.steff[i]-5780) + 1.4612e-8*(self.steff[i]-5780)**2 - 7.6345e-12*(self.steff[i]-5780)**3 - 1.7511e-15*(self.steff[i]-5780)**4)
        """
        if self.sage[i] < 10 and self.sage[i] > 0.5 and self.steff[i] < 6300 and self.steff[i] > 2300:
            turbet_iHZ = self.interp([self.sage[i], self.steff[i]])[0]**(-0.5)*self.srad[i]*(self.steff[i]/5772)**2

        else:
            turbet_iHZ = (krungh_iHZ - 0.065)**(-0.5)*self.srad[i]*(self.steff[i]/5772)**2
        """
        turbet_iHZ = (krungh_iHZ - 0.065)**(-0.5)*self.srad[i]*(self.steff[i]/5772)**2
        krungh_iHZ = krungh_iHZ**(-0.5)*self.srad[i]*(self.steff[i]/5772)**2
        kmaxgh_oHZ = (0.3507 + 5.9578e-5*(self.steff[i]-5780) + 1.6707e-9*(self.steff[i]-5780)**2 - 3.0058e-12*(self.steff[i]-5780)**3 - 5.1925e-16*(self.steff[i]-5780)**4)**(-0.5)*self.srad[i]*(self.steff[i]/5772)**2
        kemars_oHZ = (0.3207 + 5.4471e-5*(self.steff[i]-5780) + 1.5275e-9*(self.steff[i]-5780)**2 - 2.1709e-12*(self.steff[i]-5780)**3 - 3.8282e-16*(self.steff[i]-5780)**4)**(-0.5)*self.srad[i]*(self.steff[i]/5772)**2
        ciHZ = max(turbet_iHZ, krungh_iHZ)
        oiHZ = min(turbet_iHZ, krungh_iHZ)

        # create values for Earth-like planets and predictions
        elpo = 0
        elpc = 0
        predo = 0
        predc = 0
        current = datetime.now()
        
        try:    
            while (current - start).seconds < 43200:
                data = self.dyn.run_single(self.HDname[i], False)

                if data == None:    # no known planets in the system
                    for j in range(self.mc_iter):
                        sys_stable = False
                        
                        while not sys_stable:   # checking for full system stability after creation of system
                            ### create giant planets
                            pls1, pls1s, bp, bps, bpv, bpvs, pls2, pls2s, f0, sf = self.gpor_params
                            
                            pl1 = np.random.normal(pls1, pls1s)
                            pl2 = np.random.normal(pls2, pls2s)
                            bpx = np.random.normal(bp, bps)
                            bpy = 10**np.random.normal(bpv, bpvs)

                            while bpx < 0 or bpy < 0:
                                bpx = np.random.normal(bp, bps)
                                bpy = 10**np.random.normal(bpv, bpvs)
                            
                            logbins = np.logspace(np.log10(bpx*10**(-1.5)), np.log10(bpx*10**(1.5)), 51)
                            gporx = np.logspace(-1.5, 1.5, 5)*bpx
                            gpory = np.array([bpy*10**(pl1*-1.5), bpy*10**(pl1*-0.75), bpy, bpy*10**(pl2*0.75), bpy*10**(pl2*1.5)])
                            pdf_gp = 10**(np.interp(np.log10(logbins), np.log10(gporx), np.log10(gpory)))
                            intl_gp = np.cumsum(pdf_gp)
                            cdf_gp = intl_gp/intl_gp[-1]
                            inum_gp = math.floor(f0*sf)
                            frac_gp = f0*sf - inum_gp

                            if np.random.rand() < frac_gp:
                                num_gp = inum_gp + 1

                            else:
                                num_gp = inum_gp

                            if num_gp > 0:
                                per_gp = np.zeros(num_gp)
                                mas_gp = np.zeros(num_gp)
                                ecc_gp = np.zeros(num_gp)
                            
                                for k in range(num_gp):
                                    per_gp[k] = self.K3(self.cdf_draw(logbins, np.random.rand(), cdf_gp), "a", self.smass[i])
                                    mas_gp[k] = (np.random.rand()*4.5+0.5)*316
                                    ecc_gp[k] = spst.rayleigh.rvs(scale=0.03)
                    
                                l1, l2, l3 = zip(*sorted(zip(per_gp, mas_gp, ecc_gp)))
                                per_gp = np.array(l1)
                                mas_gp = np.array(l2)
                                ecc_gp = np.array(l3)

                                ### stability analysis for more than one planet
                                if num_gp > 1:
                                    deltas_gp = np.zeros(num_gp - 1)
                                    gp_stable = False

                                    while not gp_stable:                               
                                        for k in range(1, num_gp):
                                            a1 = self.K3(per_gp[k-1], "P", self.smass[i])
                                            a2 = self.K3(per_gp[k], "P", self.smass[i])
                                            deltas_gp[k-1] = 2*(a2*(1 - ecc_gp[k]) - a1*(1 + ecc_gp[k-1]))/(a2 + a1) * ((mas_gp[k] + mas_gp[k-1])/(1e6*self.smass[i]))**(-1/3)

                                        dl = spst.lognorm.rvs(spst.norm.rvs(0.40, 0.02, 1), loc=0, scale=np.exp(spst.norm.rvs(1.97, 0.03, 1)))
                                        ds = 0

                                        for k in range(len(deltas_gp)):
                                            if deltas_gp[k] > dl:
                                                ds += 1
                                        
                                        if ds == num_gp - 1:
                                            gp_stable = True
                                            break

                                        else:                        
                                            for k in range(num_gp):
                                                per_gp[k] = self.K3(self.cdf_draw(logbins, np.random.rand(), cdf_gp), "a", self.smass[i])
                                                mas_gp[k] = (np.random.rand()*4.5+0.5)*316
                                                ecc_gp[k] = spst.rayleigh.rvs(scale=0.03)

                                            l1, l2, l3 = zip(*sorted(zip(per_gp, mas_gp, ecc_gp)))
                                            per_gp = np.array(l1)
                                            mas_gp = np.array(l2)
                                            ecc_gp = np.array(l3)

                            rad_gp = spst.norm.rvs(14, 1.5, num_gp)

                            ### create small planets
                            trunc_period = 730
                            
                            ### allow for truncation (Millholland, He, & Zink 2022), determine if actually occurring and if so at what period
                            if np.random.rand() < 0.8:
                                trunc_period = np.random.rand()*200+100

                            sp_edge = min(trunc_period, (min(per_gp)/2 if num_gp > 0 else 730))
                            P_sp = np.arange(0.5, sp_edge + 0.01, 0.01)

                            if len(P_sp) == 0:
                                num_sp = 0

                            else:
                                ### find first planet period from Mulders et al. (2018) broken power law
                                pdf_psp1 = [((k/12)**1.6 if k < 12 else (k/12)**-0.9) for k in P_sp]
                                intl_psp1 = np.cumsum(pdf_psp1)
                                cdf_psp1 = intl_psp1/intl_psp1[-1]
                                perind = self.cdf_draw(P_sp, np.random.rand(), cdf_psp1)
                                per_sp = [perind]

                                ### add additional planets from Mulders et al. (2018) period ratios
                                while perind < sp_edge:
                                    spi = self.cdf_draw(self.PRgrid, np.random.rand(), self.cdfPR)*perind
                                    perind = spi

                                    if perind < sp_edge:
                                        per_sp.append(spi)

                                num_sp = len(per_sp)
                                    
                                ### create radius and mass from Mulders et al. (2018) and Bergsten et al. (2022) radii occurrence rates, He et al. (2019) intra-system radii, and Otegi et al. (2020) M-R relationship
                                rad_sp = np.zeros(num_sp)
                                mas_sp = np.zeros(num_sp)
                                ecc_sp = np.zeros(num_sp)
                            
                                R_sp = np.arange(0.5, 6.1, 0.1)
                                pdf_ssp = [((k/3.3)**-0.5 if k < 3.3 else (k/3.3)**-6) for k in R_sp]
                                intl_ssp = np.cumsum(pdf_ssp)
                                cdf_ssp = intl_ssp/intl_ssp[-1]
                                sys_rad = self.cdf_draw(R_sp, np.random.rand(), cdf_ssp)
                                cdf_rsp = spst.lognorm.cdf(R_sp, 0.3, scale=sys_rad)

                                for k in range(len(per_sp)):
                                    rad_sp[k] = self.cdf_draw(R_sp, np.random.rand(), cdf_rsp)
                                    mas_sp[k] = self.otegi_mr(rad_sp[k], "mass")
                                    ecc_sp[k] = spst.lognorm.rvs(self.ecc_sigmas[num_sp if num_sp < len(self.ecc_sigmas) else -1], scale=0.031*((num_sp+1)/5)**-1.74)

                                ### stability analysis for more than one planet
                                deltas_sp = np.zeros(num_sp - 1)
                                sp_stable = False

                                while num_sp > 1 and not sp_stable:
                                    for k in range(1, num_sp):
                                        a1 = self.K3(per_sp[k-1], "P", self.smass[i])
                                        a2 = self.K3(per_sp[k], "P", self.smass[i])
                                        deltas_sp[k-1] = 2*(a2*(1 - ecc_sp[k]) - a1*(1 + ecc_sp[k-1]))/(a2 + a1) * ((mas_sp[k] + mas_sp[k-1])/(1e6*self.smass[i]))**(-1/3)

                                    dl = spst.lognorm.rvs(spst.norm.rvs(0.40, 0.02, 1), loc=0, scale=np.exp(spst.norm.rvs(1.97, 0.03, 1)), size=1)
                                    ds = 0

                                    for k in range(len(deltas_sp)):
                                        if deltas_sp[k] > dl:
                                            ds += 1
                                
                                    if ds == num_sp - 1:
                                        sp_stable = True
                                        break

                                    else:
                                        perind = self.cdf_draw(P_sp, np.random.rand(), cdf_psp1)
                                        per_sp = [perind]

                                        ### add additional planets from Mulders et al. (2018) period ratios
                                        while perind < sp_edge:
                                            spi = self.cdf_draw(self.PRgrid, np.random.rand(), self.cdfPR)*perind
                                            perind = spi

                                            if perind < sp_edge:
                                                per_sp.append(spi)

                                        num_sp = len(per_sp)
                                        rad_sp = np.zeros(num_sp)
                                        mas_sp = np.zeros(num_sp)
                                        ecc_sp = np.zeros(num_sp)
                                    
                                        for k in range(len(per_sp)):
                                            rad_sp[k] = self.cdf_draw(R_sp, np.random.rand(), cdf_rsp)
                                            mas_sp[k] = self.otegi_mr(rad_sp[k], "mass")
                                            ecc_sp[k] = spst.lognorm.rvs(self.ecc_sigmas[num_sp if num_sp < len(self.ecc_sigmas) else -1], scale=0.031*((num_sp+1)/5)**-1.74)

                                        if num_sp <= 1:
                                            break

                                        else:
                                            deltas_sp = np.zeros(num_sp - 1)

                            ### check stability between inner and outer planets
                            if num_sp > 0 and num_gp > 0:
                                a1 = self.K3(per_sp[-1], "P", self.smass[i])
                                a2 = self.K3(per_gp[0], "P", self.smass[i])
                                delta_sg = 2*(a2*(1 - ecc_gp[0]) - a1*(1 + ecc_sp[-1]))/(a2 + a1) * ((mas_gp[0] + mas_sp[-1])/(1e6*self.smass[i]))**(-1/3)
                                dl = spst.lognorm.rvs(spst.norm.rvs(0.40, 0.02, 1), loc=0, scale=np.exp(spst.norm.rvs(1.97, 0.03, 1)), size=1)

                                if delta_sg > dl:
                                    sys_stable = True

                            else:
                                sys_stable = True

                        if num_sp > 0:
                            for k in range(num_sp):
                                if self.K3(per_sp[k], "P", self.smass[i]) < kemars_oHZ and self.K3(per_sp[k], "P", self.smass[i]) > oiHZ and rad_sp[k] < 1.8 and mas_sp[k] < 10:
                                    elpo += 1

                                    if self.K3(per_sp[k], "P", self.smass[i]) < kmaxgh_oHZ and self.K3(per_sp[k], "P", self.smass[i]) > ciHZ and rad_sp[k] < 1.5 and mas_sp[k] < 5:
                                        elpc += 1

                        current = datetime.now()

                        
                else:   # known planets in the system
                    """
                    per = [self.dyn.random_val_from_tup(k[0]) for k in data]
                    rad = [self.dyn.random_val_from_tup(k[1]) for k in data]
                    mas = [self.dyn.random_val_from_tup(k[2]) for k in data]
                    ecc = [self.dyn.random_val_from_tup(k[4]) for k in data]

                    incs = [k[3][0] for k in data if k[3][0] != "?"]
                    sys_inc = np.mean(incs)
                    num_p = len(data)
                    inc = [self.dyn.random_val_from_tup(k[3]) if k[3][0] != "?" else sys_inc + spst.lognorm.rvs(self.inc_sigmas[num_p - 1 if num_p <= len(self.inc_sigmas) else -1], loc=0, scale=1.1*((num_p+1)/5)**-1.73) for k in data]
                    res = (0, 0, 0)

                    #res = self.dynamite_prediction((self.HDname[i], per, rad, mas, inc, ecc, self.srad[i], self.smass[i]), np.arange(0.5, 730.01, 0.01), False, i)
                    """
                    data = self.dyn.run_single(self.HDname[i], True, False)

                    for j in range(len(data[0])):
                        if self.K3(data[0][j], "P", self.smass[i]) < kemars_oHZ and self.K3(data[0][j], "P", self.smass[i]) > oiHZ and data[3][j] < 1.8 and self.otegi_mr(data[3][j], "mass") < 10:
                            predo += 1

                            if self.K3(data[0][j], "P", self.smass[i]) < kmaxgh_oHZ and self.K3(data[0][j], "P", self.smass[i]) > ciHZ and data[3][j] < 1.5 and self.otegi_mr(data[3][j], "mass") < 5:
                                predc += 1

                    current = datetime.now()

                end = datetime.now()
                print("Finishing", self.HDname[i], str(i+1) + "/" + str(len(self.srad)), end, "DURATION", end - start)
                return elpo + predo, elpc + predc

            if (current - start).seconds >= 43200:
                print("NOT FINISHING", self.HDname[i], "stuck after 12 hours runtime.")
                return elpo + predo, elpc + predc
                
        except:
            print(self.HDname[i], "HAD MC ERROR, EXITING AT", datetime.now())
            traceback.print_exc()
            return elpo + predo, elpc + predc
        
if __name__ == '__main__':
    hwo_sims()
