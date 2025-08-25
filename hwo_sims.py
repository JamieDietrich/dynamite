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
from astropy import constants as const
from datetime import datetime, timedelta

class hwo_sims:

    def __init__(self):
        """Runs the simulations"""

        self.R_sun = const.R_sun.value
        self.M_sun = const.M_sun.value
        self.R_E = const.R_earth.value
        self.M_E = const.M_earth.value
        self.R_J = 10.973*self.R_E
        self.M_J = 317.8*self.M_E
        self.au = const.au.value
        self.s_B = const.sigma_sb.value
        self.G = const.G.value
        
        self.ap_arr = [1, 2, 3, 5]
        self.mp_arr = [0.414, 2, 5, 10]
        self.num_moons_mean = [2.7, 4, 7, 22.6]
        self.num_moons_std = [0.5, 1, 3, 5]
        self.moon_mass_mean = [[0.06, 0.01, 0.005, 0.001], [0.09, 0.03, 0.01, 0.005], [0.12, 0.08, 0.05, 0.007], [0.46, 0.38, 0.19, 0.08]]
        self.moon_mass_std = [[0.01, 0.005, 0.001, 0.0005], [0.015, 0.01, 0.005, 0.001], [0.02, 0.015, 0.01, 0.005], [0.1, 0.08, 0.07, 0.01]]
        self.e_mean = [0.13, 0.14, 0.15, 0.19]

        self.interpm = spint.RegularGridInterpolator((self.mp_arr, self.ap_arr), self.moon_mass_mean, bounds_error=False, fill_value=None)
        self.interps = spint.RegularGridInterpolator((self.mp_arr, self.ap_arr), self.moon_mass_std, bounds_error=False, fill_value=None)

        # Get stellar parameters from new database
        self.HDname = np.loadtxt("spores_catalog_v2.1.0.csv", delimiter=',', dtype=str, skiprows=1, usecols=(5))
        names = np.loadtxt("names.csv", delimiter=',', dtype=str)
        self.altnames = {}
        
        for i in names:
            if i[1] != '':
                self.altnames[i[0]] = i[1]

        self.steff, self.srad, self.smass, self.sage = np.loadtxt("spores_catalog_v2.1.0.csv", delimiter=',', dtype=float, unpack=True, skiprows=1, usecols=(26, 32, 34, 221))

        # Create/use planetary systems and test for additional planets and CHZ value
        self.per = []
        self.rad = []

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
        saved = False

        if saved:
            print(datetime.now(), "Loading files")
            kp = np.loadtxt("systems.txt", usecols=(1), dtype=str)
            elpo, predo, tnm, tmgpo, elpc, predc, tmgpc = np.loadtxt("data.txt", unpack=True)
            self.per, self.rad = np.loadtxt("planets.txt", unpack=True)
            
        else:
            try:
                self.dyn = dynamite(False)
                res = self.run_new_mp(self.run_sims, np.arange(len(self.srad)), ())
                results = []

                for j in range(len(res[0])-2):
                    res1 = []
                    
                    for i in range(len(res)):
                        res1.append(res[i][j])

                    results.append(np.array(res1))

                for i in range(len(res)):   
                    for j in range(len(res[i][-2])):
                        self.per.append(res[i][-2][j])
                        self.rad.append(res[i][-1][j])

                kp, elpo, predo, tnm, tmgpo, elpc, predc, tmgpc = results

                with open("values.txt", 'w') as f:
                    f.write("System\tKnown Planets?\tOptimistic Simulated\tOptimistic Predicted\tOptimistic Moon\tConservative Simulated\tConservative Predicted\tConservative Moon\tTotal Moon\n")
                    
                    for i in range(len(elpo)):
                        f.write(self.HDname[i] + "\t" + kp[i] + "\t" + str(elpo[i]/self.mc_iter) + "\t" + str(predo[i]/self.mc_iter) + "\t" + str(tmgpo[i]) + "\t" + str(elpc[i]/self.mc_iter) + "\t" + str(predc[i]/self.mc_iter) + "\t" + str(tmgpc[i]) + "\t" + str(tnm[i]) + "\n")

                np.savetxt("systems.txt", np.c_[self.HDname, [self.altnames[i] if i in self.altnames.keys() else "" for i in self.HDname], kp], fmt="%s", delimiter='\\')
                np.savetxt("data.txt", np.c_[elpo/self.mc_iter, predo/self.mc_iter, tnm, tmgpo, elpc/self.mc_iter, predc/self.mc_iter, tmgpc])
                np.savetxt("planets.txt", np.c_[self.per, self.rad])
                """
                print("The stars in this list have an optimistic estimated eta_Earth of", round((np.mean(oes))/self.mc_iter, 4), "+/-", round(np.std(oes)/self.mc_iter, 4), "and a conservative estimated eta_Earth of", round(np.mean(ces)/self.mc_iter, 4), "+/-", round(np.std(ces)/self.mc_iter, 4))
                print("The star with the highest estimated optimistic eta_Earth of", round(max(oes)/self.mc_iter, 4),"is", self.HDname[np.where(oes == max(oes))[0][0]], "and the star with the highest estimated conservative eta_Earth of", round(max(ces)/self.mc_iter, 4), "is", self.HDname[np.where(ces == max(ces))[0][0]])
                exit()
                exit()
                """
            except Exception as e:
                traceback.print_exc()
                exit()

        self.per = np.array(self.per).flatten()
        self.rad = np.array(self.rad).flatten()
        
        print(datetime.now(), "Making plots")
        plt.figure(figsize=(12, 8))
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', labelsize=18)
        plt.hist(elpo/self.mc_iter, bins=np.logspace(-4, 0, 21), label="Simulated")
        plt.hist(predo/self.mc_iter, bins=np.logspace(-4, 0, 21), label="Predicted")
        plt.hist(tmgpo/self.mc_iter, bins=np.logspace(-4, 0, 21), label="Moon")
        plt.xscale("log")
        plt.yscale("log")
        plt.legend(fontsize=22)
        plt.xlabel("Likelihood", fontsize=24)
        plt.ylabel("Counts", fontsize=24)
        plt.title("Optimistic estimates for temperate terrestrial planets", fontsize=28)
        plt.savefig("optimistic_eta_hist.png")
        plt.clf()
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', labelsize=18)
        plt.hist(elpc/self.mc_iter, bins=np.logspace(-4, 0, 21), label="Simulated")
        plt.hist(predc/self.mc_iter, bins=np.logspace(-4, 0, 21), label="Predicted")
        plt.hist(tmgpc/self.mc_iter, bins=np.logspace(-4, 0, 21), label="Moon")
        plt.legend(fontsize=22)
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("Likelihood", fontsize=24)
        plt.ylabel("Counts", fontsize=24)
        plt.title("Conservative estimates for temperate terrestrial planets", fontsize=28)
        plt.savefig("conservative_eta_hist.png")
        plt.clf()
        plt.hist(np.array(self.per), bins=np.logspace(-0.5, 5.5, 61))
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', labelsize=18)
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("Period (d)", fontsize=24)
        plt.ylabel("Counts", fontsize=24)
        plt.title("All planet periods (simulated and known)", fontsize=30)
        plt.savefig("sim_per_hist.png")
        plt.clf()
        plt.hist(np.array(self.per[np.where(self.rad < 6)]), bins=np.logspace(-0.5, 5.5, 61))
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', labelsize=18)
        plt.xlim(0.3, 1000)
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("Period (d)", fontsize=24)
        plt.ylabel("Counts", fontsize=24)
        plt.title("All small planet periods (simulated and known)", fontsize=30)
        plt.savefig("sim_per_sp_hist.png")
        plt.clf()
        plt.hist(np.array(self.per[np.where(self.rad > 9)]), bins=np.logspace(-0.5, 5.5, 61))
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', labelsize=18)
        plt.xlim(1, 5e5)
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("Period (d)", fontsize=24)
        plt.ylabel("Counts", fontsize=24)
        plt.title("All giant planet periods (simulated and known)", fontsize=30)
        plt.savefig("sim_per_gp_hist.png")
        plt.clf()
        plt.hist(np.array(self.rad), bins=np.logspace(-1, 1.5, 51))
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', labelsize=18)
        plt.xlim(0.4, 40)
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel(r"Planet Radius ($R_\oplus$)", fontsize=24)
        plt.ylabel("Counts", fontsize=24)
        plt.title("All planet radii (simulated and known)", fontsize=30)
        plt.savefig("sim_rad_hist.png")
        plt.close()
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

            elif measurement > 131.61:
                return 14.357*(measurement/131.61)**-0.04

            return 1.03*measurement**0.29

        elif predict == "mass":
            if measurement == 0:
                return 0
            
            M_r = 0.9*measurement**3.45
            M_v = 1.74*measurement**1.58
            M_j = 131.61*(measurement/14.357)**(-25)
            
            if M_r > 5:
                if M_j < 4131.4:
                    return M_j
                
                return M_v

            return M_r

    def make_moons(self, i, per, rad, mas, ecc):
        """Makes moons for giant planets"""

        tnm = 0
        tmgpo = 0
        tmgpc = 0

        for k in range(len(per)):
            if mas[k] > 0.414*self.M_J/self.M_E:
                ap = self.K3(per[k], "P", self.smass[i])
                nmm = np.interp(ap, self.ap_arr, self.num_moons_mean)
                nms = np.interp(ap, self.ap_arr, self.num_moons_std)
                num_moons = int(round(spst.norm.rvs(nmm, nms, 1)[0], 0))
                
                if num_moons > 0:
                    tnm += num_moons
                    mmm = self.interpm((mas[k]*self.M_E/self.M_J, ap))
                    
                    if mmm <= 0:
                        mmm = 0.001

                    rmm = self.otegi_mr(mmm, "radius")
                    mms = self.interps((mas[k]*self.M_E/self.M_J, ap))

                    if mms <= 0:
                        mms = 0.0005
                        
                    mem = np.interp(mas[k]*self.M_E/self.M_J, self.mp_arr, self.e_mean)
                    roche = rmm*self.R_E*(2*mas[k]/mmm)**(1/3)
                    hill = ap*self.au*(1-ecc[k])*(mmm/(3*(mmm + mas[k])))**(1/3)
                    
                    am = np.logspace(np.log10(roche), np.log10(0.4*hill), num_moons)
                    Mm = np.zeros(num_moons)
                    Rm = np.zeros(num_moons)
                    em = np.zeros(num_moons)
                    l = 0
                    m = 0
                    
                    while l < num_moons and m < 1000:
                        Mm = spst.norm.rvs(mmm, mms, 1)[0]
                        em = spst.norm.rvs(mem, 0.1, 1)[0]
                        
                        if Mm <= 0 or em >= 1:
                            pass

                        else:
                            if em < 0:
                                em = 0
                            
                            if Mm > 0.1 and em < 0.6:
                                Rm = self.otegi_mr(Mm, "radius")
                                cas = self.tidal_heat_flux(i, per[k], ecc[k], rad[k], mas[k], am[l], em, Rm, Mm)

                                if cas == 2:
                                    tmgpc += 1

                                elif cas == 1:
                                    tmgpo += 1
                            
                            l += 1
                            m += 1

                            if m == 1000:
                                break
                                print("WHILE LOOP ISSUE")

        return tnm, tmgpo, tmgpc
        
    def tidal_heat_flux(self, i, pp, ep, Rp, Mp, am, em, Rm, Mm):
        """Calculates the tidal heat flux on a moon due to its host planet"""
        
        AoE = 0
        AiE = 0.28
        ABE = 0.294
        ABJ = 0.503
        xm = 0.98
        fm = 4

        Mp = Mp*self.M_J
        Rp = Rp*self.R_J
        ap = self.K3(pp, "P", self.smass[i])*self.au
        Mm = Mm*self.M_E
        Rm = Rm*self.R_E
        g = self.G*Mm/Rm**2
        """
        pm = self.K3(am/au, "a", Mp/M_sun)*86400
        nm = math.sqrt(G*Mp/am**3)
        omega = 1/pm
        rho = Mm/(4*math.pi*Rm**3/3)
        temp = np.linspace(1200, 1600, 81)
        eta = np.median([1.6e5*np.exp(40000/t) for t in temp])
        mu = 5e10
        imk2 = 57*eta*omega/(4*rho*g*Rm*(1+(1+19*mu/(2*rho*g*Rm))**2*(eta*omega/mu)**2))
        hm = 10.5*imk2*Rm**5*nm**5*em**2/(G*4*math.pi*Rm**2)
        """
        k2 = 0.3
        ts = 638
        Zs = 3*self.G**2*k2*Mp**2*(Mp+Mm)*Rm**5*ts/am**9
        B = math.sqrt(1-em**2)
        f1 = 1+15.5*em**2+255*em**4/8+185*em**6/16+25*em**8/64
        f2 = 1+7.5*em**2+45*em**4/8+5*em**6/16
        f5 = 1+3*em**2+3*em**4/8
        hm = Zs*(f1-f2/f5)/(B*4*math.pi*Rm**2)

        rh = 1e14*Mm/self.M_E*math.exp(-self.sage[i]/1.368)/(4*math.pi*Rm**2)
        
        Lstar = 4*math.pi*(self.srad[i]*self.R_sun)**2*self.s_B*self.steff[i]**4
        Lp = 9.467e-12*self.G*(Mp/Rp)**2 + Lstar*(1-ABJ)/(4*math.pi*ap**2*math.sqrt(1-ep**2))
        Fsb = Lstar*(1-AoE)/(4*math.pi*ap**2*math.sqrt(1-ep**2))
        Fsrd = xm/4+math.pi*Rp**2*ABE/(fm*2*am**2)
        Fs = Fsb*Fsrd
        Fp = Lp*(1-AiE)/(4*math.pi*am**2*fm*math.sqrt(1-em**2))
        F = Fs + Fp + hm + rh

        o = 0.7344
        l = -4.065e4
        R = 8.31446262
        Pr = 610.616
        Tr = 273.13
        Pp = Pr*math.exp(l/(R*Tr))
        P0 = 104
        k0 = 0.055

        Frg = o*self.s_B*(l/(R*math.log(Pp/math.sqrt(2*P0*g/k0))))**4        

        marr = [0.1, 1, 5]
        cihz_arr = [[0.99, 1.107, 1.188], [1.209e-4, 1.332e-4, 1.433e-4], [1.404e-8, 1.58e-8, 1.707e-8], [-7.418e-12, -8.308e-12, -8.968e-12], [-1.713e-15, -1.931e-15, -2.084e-15]]
        cihz_coeff = [np.interp(Mm/self.M_E, marr, c) for c in cihz_arr]
        cihz_lim = cihz_coeff[0] + cihz_coeff[1]*(self.steff[i]-5780) + cihz_coeff[2]*(self.steff[i]-5780)**2 + cihz_coeff[3]*(self.steff[i]-5780)**3 + cihz_coeff[4]*(self.steff[i]-5780)**4
        cohz_lim = 0.356 + 6.171e-5*(self.steff[i]-5780) + 1.698e-9*(self.steff[i]-5780)**2 - 3.198e-12*(self.steff[i]-5780)**3 - 5.575e-16*(self.steff[i]-5780)**4
        oihz_lim = 1.776 + 2.136e-4*(self.steff[i]-5780) + 2.533e-8*(self.steff[i]-5780)**2 - 1.332e-11*(self.steff[i]-5780)**3 - 3.097e-15*(self.steff[i]-5780)**4
        oohz_lim = 0.32 + 5.547e-5*(self.steff[i]-5780) + 1.526e-9*(self.steff[i]-5780)**2 - 2.874e-12*(self.steff[i]-5780)**3 - 5.011e-16*(self.steff[i]-5780)**4
        Seff = F*cihz_coeff[0]/Frg

        if Seff > cohz_lim and Seff < cihz_lim:
            return 2

        elif Seff > oohz_lim and Seff < oihz_lim:
            return 1

        return 0

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
        """Runs new multiprocessing code needed for iOS users"""

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
        """Runs the simulated systems with multiprocessing"""

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
        tnm = 0
        tmgpo = 0
        tmgpc = 0
        pers = []
        rads = []
        
        current = datetime.now()
        
        try:    
            while (current - start).seconds < 43200:
                val1 = self.dyn.run_single(self.HDname[i], False, True)
                val2 = self.dyn.run_single(self.HDname[i][:-2], False, True)
                
                if val1 != None:
                    data1 = val1

                elif val2 != None:
                    data1 = val2

                else:
                    data1 = None

                if self.HDname[i] in self.altnames.keys():
                    val1 = self.dyn.run_single(self.altnames[self.HDname[i]], False, True)
                    val2 = self.dyn.run_single(self.altnames[self.HDname[i]][:-2], False, True)
                    val3 = self.dyn.run_single(self.altnames[self.HDname[i]] + " A", False, True)

                    if val1 != None:
                        data2 = val1

                    elif val2 != None:
                        data2 = val2

                    elif val3 != None:
                        data2 = val3
                        
                    else:
                        data2 = None

                else:
                    data2 = None

                if data1 == None and data2 == None:    # no known planets in the system
                    print(self.HDname[i], "not found in database, creating simulated systems including moons of giant planets")
                    
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
                                    mas_gp[k] = (np.random.rand()*9.5+0.5)*self.M_J/self.M_E
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
                                                mas_gp[k] = (np.random.rand()*4.5+0.5)*self.M_J/self.M_E
                                                ecc_gp[k] = spst.rayleigh.rvs(scale=0.03)

                                            l1, l2, l3 = zip(*sorted(zip(per_gp, mas_gp, ecc_gp)))
                                            per_gp = np.array(l1)
                                            mas_gp = np.array(l2)
                                            ecc_gp = np.array(l3)

                            rad_gp = [self.otegi_mr(k, "radius") for k in mas_gp]

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

                                    while rad_sp[k] > 6:
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

                                            while rad_sp[k] > 6:
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

                                pers.append(per_sp[k])
                                rads.append(rad_sp[k])
                        
                        if num_gp > 0:
                            tnmi, tmgpoi, tmgpci = self.make_moons(i, per_gp, rad_gp, mas_gp, ecc_gp)
                            tnm += tnmi
                            tmgpo += tmgpoi
                            tmgpc += tmgpci

                            for k in range(num_gp):
                                pers.append(per_gp[k])
                                rads.append(rad_gp[k])
                                
                        current = datetime.now()

                        
                else:   # known planets in the system

                    if data1 != None:
                        targ = self.HDname[i]

                    elif val1 != None:
                        targ = self.altnames[self.HDname[i]]

                    elif val2 != None:
                        targ = self.altnames[self.HDname[i]][:-2]

                    else:
                        targ = None
                        
                    data = self.dyn.run_single(targ, True, False)
                    
                    for k in range(len(data[-4])):
                        pers.append(data[-5][k][0])
                        rads.append(data[-4][k][0])

                    print("Checking", self.HDname[i], "for habitable moons of giant planets")

                    for j in range(len(data[0])):
                        pxp = [data[0][j], data[3][j], self.otegi_mr(data[3][j], "mass"), data[9][j]]

                        if self.K3(pxp[0], "P", self.smass[i]) < kemars_oHZ and self.K3(pxp[0], "P", self.smass[i]) > oiHZ and pxp[1] < 1.8 and pxp[2] < 10:
                            predo += 1

                            if self.K3(pxp[0], "P", self.smass[i]) < kmaxgh_oHZ and self.K3(pxp[0], "P", self.smass[i]) > ciHZ and pxp[1] < 1.5 and pxp[2] < 5:
                                predc += 1

                        msini = [data[-3][k][0] for k in range(len(data[-3]))]
                        inclim = np.arcsin(self.M_E*max(msini)/(13*self.M_J))
                        sysinc = np.random.rand()*(math.pi/2 - inclim) + inclim
                            
                        per = [self.dyn.random_val_from_tup(k) for k in data[-5]]
                        rad = [self.dyn.random_val_from_tup(k) for k in data[-4]]
                        mas = [self.dyn.random_val_from_tup(data[-3][k])/(np.sin(sysinc) if data[-2][k] == "?" else 1) for k in range(len(data[-3]))]
                        ecc = [self.dyn.random_val_from_tup(k) for k in data[-1]]

                        per.append(pxp[0])
                        rad.append(pxp[1])
                        mas.append(pxp[2])
                        ecc.append(pxp[3])

                        tnmi, tmgpoi, tmgpci = self.make_moons(i, per, rad, mas, ecc)
                        tnm += tnmi
                        tmgpo += tmgpoi
                        tmgpc += tmgpci

                    current = datetime.now()
                
                end = datetime.now()
                print("Finishing", self.HDname[i], str(i+1) + "/" + str(len(self.srad)), end, "DURATION", end - start)
                return "No" if data1 == None and data2 == None else "Yes", elpo, predo, tnm if tnm > 0 else 1, tmgpo, elpc, predc, tmgpc, pers, rads

            if (current - start).seconds >= 43200:
                print("NOT FINISHING", self.HDname[i], "stuck after 12 hours runtime")
                return "No" if data1 == None and data2 == None else "Yes", elpo, predo, tnm if tnm > 0 else 1, tmgpo, elpc, predc, tmgpc, pers, rads
                
        except:
            print(self.HDname[i], "HAD MC ERROR, EXITING AT", datetime.now())
            traceback.print_exc()
            return "No" if data1 == None and data2 == None else "Yes", elpo, predo, tnm if tnm > 0 else 1, tmgpo, elpc, predc, tmgpc, pers, rads
        
if __name__ == '__main__':
    hwo_sims()
