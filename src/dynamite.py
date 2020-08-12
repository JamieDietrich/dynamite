### DYNAmical Multi-planet Injection TEster (DYNAMITE) ###
### Main File ###
### Jeremy Dietrich ###
### jdietrich1@email.arizona.edu ###
### 2020 August 12 ###
### Version 1.3 ###
### Dietrich & Apai (2020), Astronomical Journal ###
### https://doi.org/10.3847/1538-3881/aba61d ###

import os
import ast
import sys
import math
import itertools
import numpy as np
from PPR import PPR
import scipy.stats as spst
from datetime import datetime
import astropy.constants as const
from dynamite_plots import dynamite_plots
from dynamite_targets import dynamite_targets
from mrexo import predict_from_measurement as pfm

class dynamite:

    def __init__(self, cfname="dynamite_config.txt"):
        """Runs the script"""

        print(datetime.now(), "Initiating DYNAMITE")
        np.random.seed(1)

        self.config_parameters = {}

        try:
            config_data = np.loadtxt(cfname, dtype=str, delimiter='::')

        except IOError:
            print("Error, configuration file not found!")
            exit()

        for i in range(len(config_data)):
            self.config_parameters[config_data[i, 0]] = config_data[i, 1] if config_data[i, 1].find("[") == -1 else ast.literal_eval(config_data[i, 1])

        targets_dict = dynamite_targets(self.config_parameters).get_targets(self.config_parameters["mode"], self.config_parameters["system"], self.config_parameters["radmax"], self.config_parameters["removed"])

        if len(targets_dict) == 0:
            print("Error: No targets selected!")
            exit()
       
        self.ppr = PPR((self, None))

        if self.config_parameters["saved"] == "False":
            if os.path.exists("table_" + self.config_parameters["mode"] + "_" + self.config_parameters["period"] + ".txt"):
                os.remove("table_" + self.config_parameters["mode"] + "_" + self.config_parameters["period"] + ".txt")

            if self.config_parameters["mode"] == "single":
                R_star, Rse, M_star, Mse, target, target_name = self.set_up(targets_dict, self.config_parameters["system"])
                targ_rem = []

                for i in range(len(target) - 1 if len(self.config_parameters["removed"]) > 0 else len(target)):
                    if target[i][2] not in self.config_parameters["removed"]:
                        targ_rem.append(target[i])

                if len(self.config_parameters["additional"][0]) > 0:
                    for i in self.config_parameters["additional"]:
                        x = []

                        for j in range(len(i) - 1):
                            if isinstance(i[j], tuple):
                                if self.config_parameters["use_mass"] == "True":
                                    x.append(self.mr_convert(i[j][0]))

                                else:
                                    x.append(i[j][0])

                            else:
                                x.append(i[j])

                        targ_rem.append(x)

                    for i in self.config_parameters["unconfirmed"]:
                        x = []

                        for j in range(len(i) - 1):
                            if isinstance(i[j], tuple):
                                if self.config_parameters["use_mass"] == "True":
                                    x.append(self.mr_convert(i[j][0]))

                                else:
                                    x.append(i[j][0])

                            else:
                                x.append(i[j])

                        targ_rem.append(x)

                targ_rem = np.array(targ_rem)
                target = targ_rem[targ_rem[:, 2].argsort()]
                data = self.run_monte_carlo(R_star, Rse, M_star, Mse, target, target_name) + ([target[i][2] for i in range(len(target))], [target[i][1] for i in range(len(target))])
                np.savez("saved_data.npz", data=data)

            else:
                targlist = []

                for tn in targets_dict.keys():
                    if self.config_parameters["mode"] == "all":
                        targlist.append(tn)

                    elif self.config_parameters["mode"] == "tess" and tn.find("TOI") != -1:
                        targlist.append(tn)

                    elif self.config_parameters["mode"] == "kepler" and tn.find("Kepler") != -1:
                        targlist.append(tn)

                    elif self.config_parameters["mode"] == "k2" and tn.find("K2") != -1:
                        targlist.append(tn)

                    elif self.config_parameters["mode"] == "test" and tn.find("test 3") != -1:
                        targlist.append(tn)

                Pk, P, PP, per, Rk, R, PR, ik, il, Pin, deltas, ratios, tdm, tdle, tdue, tpm, tple, tpue, targets, pers, rads = datavars = self.ppr.create_processes("mt_mc", (targets_dict, targlist), -len(targlist), self.process_data)

                np.savez("saved_data.npz", data=datavars)

        elif self.config_parameters["saved"] == "True":
            with np.load("saved_data.npz", allow_pickle=True) as data:
                Pk, P, PP, per, Rk, R, PR, ik, il, Pinc, deltas, ratios, tdm, tdle, tdue, tpm, tple, tpue, targets, pers, rads = data["data"]

        if self.config_parameters["plot"] == "True":
            if self.config_parameters["saved"] == "False" and self.config_parameters["mode"] == "single":
                Pk, P, PP, per, Rk, R, PR, ik, il, Pinc, deltas, ratios, tdm, tdle, tdue, tpm, tple, tpue, targets, pers, rads = data

            elif self.config_parameters["saved"] == "False" and self.config_parameters["mode"] != "single":
                Pk, P, PP, per, Rk, R, PR, ik, il, Pinc, deltas, ratios, tdm, tdle, tdue, tpm, tple, tpue, targets, pers, rads = datavars

            print(datetime.now(), "Creating Plots")
            plots = dynamite_plots(Pk, P, PP, per, Rk, R, PR, ik, il, Pinc, deltas, ratios, tdm, tdle, tdue, tpm, tple, tpue, targets, pers, rads, cfname)
            print(datetime.now(), "Finishing DYNAMITE")



    def mt_mc(self, targets_dict, targlist, i):
        """Runs the Monte Carlo code on multiple threads"""

        R_star, Rse, M_star, Mse, target, target_name = self.set_up(targets_dict, targlist[i])  

        try:
            target = target[target[:, 2].argsort()]

        except IndexError:
            print("ERROR IN DICT ON TARGET", target_name)

        data = self.run_monte_carlo(R_star, Rse, M_star, Mse, target, target_name) + ([target[i][2] for i in range(len(target))], [target[i][1] for i in range(len(target))])

        return {i:data}
    


    def process_data(self, data):
        """Processes data for the multithreading component"""

        return tuple([list(z) for z in zip(*itertools.chain([data[k] for k in data]))])



    def set_up(self, targets_dict, target):
        """Sets up target"""

        def get_arccos(star_pars, planet_pars):
            return round(np.arccos(planet_pars[0]/(self.K3(planet_pars[1], star_pars[2])/(star_pars[0]*const.R_sun.cgs.value)))*180/math.pi, 3)
           
        t = list(targets_dict[target])

        for x in range(len(t)):
            for y in range(len(t[x])):
                if isinstance(t[x][y], tuple):
                    if isinstance(t[x][y][0], str):
                        t[x][y] = locals()[t[x][y][0]](t[0],t[x][y][1])

                    else:
                        if t[x][y][1] == "Mass":
                            t[x][y] = self.mr_convert(t[x][y][0])

                        elif t[x][y][1] == "Radius":
                            t[x][y] = t[x][y][0]

        return t[0][0], t[0][1], t[0][2], t[0][3], np.array([t[i][:-1] for i in range(1, len(t))]), target


    def mr_convert(self, meas):
        """Runs conversion from mass to radius"""

        if self.config_parameters["mass_radius"] == "mrexo":
            return pfm(meas, predict="radius", dataset="kepler")[0]

        elif self.config_parameters["mass_radius"] == "otegi":
            return self.otegi_mr(meas, "radius")




    def process_inc_data(self, data):
        """Processes data for the multithreading component"""
        
        return tuple([list(itertools.chain(*i)) for i in zip(*itertools.chain([data[k] for k in data]))])



    def inc_test(self, il, incn, rylgh, j):
        """Tests the best system inclination."""

        ibs = []
        fib = []

        for k in range(len(incn)):
            test = 0

            for m in range(len(incn[k])):
                test += spst.rayleigh.pdf(abs(incn[k][m]-il[j]), rylgh)

            ibs.append(il[j])
            fib.append(test)
            
        return {j: (ibs, fib)}



    def run_monte_carlo(self, R_star, Rse, M_star, Mse, target, target_name):
        """Runs the Monte Carlo analysis."""

        inc = [target[i][0] for i in range(len(target))]
        rad = [target[i][1] for i in range(len(target))]
        per = [target[i][2] for i in range(len(target))]
        ratios = []

        for i in range(1, len(per)):
            ratios.append(per[i]/per[i-1])

        p0 = min(per)
        r1 = min(rad)
        r2 = max(rad)
        P = np.arange(0.5, 730.1, 0.1)

        print(datetime.now(), "Creating Period Distributions for", target_name)

        if self.config_parameters["period"] == "epos":
            PP, deltas = self.epos_pers(p0, per, rad, P, M_star)

        elif self.config_parameters["period"] == "syssim":
            PP, deltas = self.syssim_pers(per, rad, P, M_star)

        print(datetime.now(), "Creating Planet Radius Distributions for", target_name)

        if self.config_parameters["radius"] == "epos":
            R, PR, cdfR = self.epos_rads(r1, r2)

        elif self.config_parameters["radius"] == "syssim":
            R, PR, cdfR = self.syssim_rads(self.config_parameters["radtype"], rad)

        print(datetime.now(), "Creating Inclination Distributions for", target_name)
        il = np.linspace(0, 180.1, 1802)
        fi = np.zeros(len(il))
        rylgh = 2
        ibs = []
        fib = []
        incn = []

        for case in [[False] + list(t) for t in list(itertools.product([False,True], repeat=len(inc)-1))]:
            incn.append([180-inc[i] if case[i] else inc[i] for i in range(0, len(inc))])

        ibs, fib = self.ppr.create_processes("inc_test", (il, incn, rylgh), -len(il), self.process_inc_data)
        ib1 = ibs[np.where(fib == max(fib))[0][0]]

        if ib1 > 90:
            ib1 = 180 - ib1

        ib = ibs[np.where(fib == max(fib))[0][0]]

        if ib > 90:
            ib = 180 - ib

        inew = np.linspace(0, 10, 101)
        finew = spst.rayleigh.pdf(inew + ib, ib, rylgh)

        if len(inc) == 2:
            finew = finew*0.62

            for j in range(len(il)):
                fi[j] = np.sin(il[j]*math.pi/180)*19/75

            fi = fi*19/(75*np.trapz(fi, il))

        elif len(inc) == 3:
            finew = finew*0.81

            for j in range(len(il)):
                fi[j] = np.sin(il[j]*math.pi/180)*19/150

            fi = fi*19/(150*np.trapz(fi, il))

        i_ib = np.where(np.isclose(il,ib))[0][0]

        for j in range(len(inew)):
            fi[i_ib + j] += finew[j]

        cdfi = np.array([1 - math.exp(-(inew[j])**2/(2*(rylgh)**2)) for j in range(len(inew))])
        Pinc = fi/2

        print(datetime.now(), "Running Monte Carlo for", target_name)
        Pk = []
        Rk = []
        ik = []

        for k in range(int(self.config_parameters["MC_chain"])):
            for j in range(len(PP)):
                if np.random.rand() < PP[j]:
                    """
                    ib = ib1 + np.random.normal(0, 10, 1)

                    if ib > 180:
                        ib = 360 - ib

                    if ib < 0:
                        ib = -ib

                    inew = np.linspace(0, 10, 101)
                    finew = spst.rayleigh.pdf(inew + ib, ib, rylgh)

                    if len(inc) == 2:
                        finew = finew*0.62

                        for j in range(len(il)):
                            fi[j] = np.sin(il[j]*math.pi/180)*19/75

                        fi = fi*19/(75*np.trapz(fi, il))

                    elif len(inc) == 3:
                        finew = finew*0.81

                        for j in range(len(il)):
                            fi[j] = np.sin(il[j]*math.pi/180)*19/150

                        fi = fi*19/(150*np.trapz(fi, il))

                    i_ib = np.where(abs(il - ib) < 0.1)[0][0]

                    for j in range(len(inew)):
                        fi[i_ib + j] += finew[j]

                    cdfi = np.array([1 - math.exp(-(inew[j])**2/(2*(rylgh)**2)) for j in range(len(inew))])
                    Pinc = fi/2
                    """
                    Pk.append(P[j])
                    RR = np.random.rand()
                    Rk.append(R[np.where(RR - cdfR < 0)[0][0]])
                    ii = np.random.rand()
                    iso = np.random.rand()

                    if (len(inc) == 2 and iso > 0.38) or (len(inc) == 3 and iso > 0.19) or (len(inc) >= 4):
                        if len(np.where(ii - cdfi < 0)[0]) == 0:
                            ik.append(inew[-1] + ib)

                        else:
                            ik.append(inew[np.where(ii - cdfi < 0)[0][0]] + ib)

                    else:
                        ik.append(np.arccos(np.random.rand()*2 - 1)*180/math.pi)

        print(datetime.now(), "Calculating Best Fit Predictions for", target_name)

        Pis = np.hstack(np.array([np.arange(0.5, 1.001, 0.001), np.arange(1.01, 10.01, 0.01), np.arange(10.1, 100.1, 0.1), np.arange(101,731,1)]))

        if self.config_parameters["period"] == "epos":
            PPi, _ = self.epos_pers(p0, per, rad, Pis, M_star)

            if self.config_parameters["mode"] == "tess":
                low_gap_P = ["TOI 561", "TOI 431", "TOI 1238", "TOI 732", "TOI 696", "TOI 175", "TOI 663", "TOI 1469", "TOI 2096", "TOI 1260", "TOI 270", "TOI 396", "TOI 836", "TOI 411", "TOI 1130", "TOI 1269", "TOI 1453", "TOI 714", "TOI 1749", "TOI 125", "TOI 1445", "TOI 1438", "TOI 119", "TOI 763", "TOI 2084", "TOI 1136", "TOI 1803", "TOI 1064", "TOI 266", "TOI 178", "TOI 776", "TOI 1339", "TOI 214", "TOI 700", "TOI 1266", "TOI 1812", "TOI 553", "TOI 699", "TOI 1277"]
                interior = ["TOI 2095", "TOI 282"]
        
                if target_name in low_gap_P or target_name in interior:
                    PPz = PPi

                else:
                    Pz = np.where(PPi == 0)[0]
                    PPz = PPi[Pz[1]:Pz[-1]]

            else:
                PPz = PPi

            Pm = Pis[np.where(PPi == np.amax(PPz))[0][0]]
            PPm = np.amax(PPz)

        elif self.config_parameters["period"] == "syssim":
            PPi, _ = self.syssim_pers(per, rad, Pis, M_star)
            Pm = Pis[np.where(PPi == np.amax(PPi))[0][0]]
            PPm = np.amax(PPi)

        if len(Pis[np.where((Pis < Pm) & (PPi < 0.606*PPm))]) > 0:
            Ple = Pm - Pis[np.where((Pis < Pm) & (PPi < 0.606*PPm))][-1]

        else:
            Ple = 0.606*Pm

        if len(Pis[np.where((Pis > Pm) & (PPi < 0.606*PPm))]) > 0:
            Pue = Pis[np.where((Pis > Pm) & (PPi < 0.606*PPm))][0] - Pm

        else:
            Pue = Pm/0.606

        Rm = np.percentile(Rk, 50)
        Rle = Rm - np.percentile(Rk, 16)
        Rue = np.percentile(Rk, 84) - Rm
        tdm = (Rm*const.R_earth.cgs.value/(R_star*const.R_sun.cgs.value))**2*1e6
        tdle = 2*(Rm*const.R_earth.cgs.value/(R_star*const.R_sun.cgs.value))**2*math.sqrt((Rle/Rm)**2 + (Rse/R_star)**2)*1e6
        tdue = 2*(Rm*const.R_earth.cgs.value/(R_star*const.R_sun.cgs.value))**2*math.sqrt((Rue/Rm)**2 + (Rse/R_star)**2)*1e6
        ntrans = 0
        ntl = 0
        ntu = 0

        print(datetime.now(), "Running Transit Calculations for", target_name)

        for j in range(len(ik)):
            if math.cos(ik[j]*math.pi/180) < (R_star*const.R_sun.cgs.value + Rm*const.R_earth.cgs.value)/self.K3(Pm, M_star):
                ntrans += 1

            if math.cos(ik[j]*math.pi/180) < ((R_star - Rse)*const.R_sun.cgs.value + (Rm - Rle)*const.R_earth.cgs.value)/self.K3(Pm, (M_star + Mse)):
                ntl += 1

            if math.cos(ik[j]*math.pi/180) < ((R_star + Rse)*const.R_sun.cgs.value + (Rm + Rue)*const.R_earth.cgs.value)/self.K3(Pm, (M_star - Mse)):
                ntu += 1

        print(datetime.now(), "Writing out Best Values for", target_name)
        tpm = ntrans/len(ik)
        tple = max(1e-3, ntrans/len(ik) - ntl/len(ik)) if (tpm != 0 and tpm != 1) else ntrans/len(ik) - ntl/len(ik)
        tpue = max(1e-3, ntu/len(ik) - ntrans/len(ik)) if (tpm != 0 and tpm != 1) else ntu/len(ik) - ntrans/len(ik)
        Pm = round(Pm, 1) if Pm > 10 else round(Pm, 3) if Pm < 1 else round(Pm, 2)
        Pue = round(Pue, 1) if Pue > 10 else round(Pue, 3) if Pue < 1 else round(Pue, 2)
        Ple = round(Ple, 1) if Ple > 10 else round(Ple, 3) if Ple < 1 else round(Ple, 2)
        Rm = round(Rm, 3) if Rm < 1 else round(Rm, 2)
        Rue = round(Rue, 3) if Rue < 1 else round(Rue, 2)
        Rle = round(Rle, 3) if Rle < 1 else round(Rle, 2)
        tdm = int(round(tdm/10, 0))*10 if tdm > 1000 else int(round(tdm, 0))
        tdue = int(round(tdue/100, 0))*100 if tdue > 10000 else int(round(tdue/10, 0))*10 if tdue > 1000 else int(round(tdue, 0))
        tdle = int(round(tdle/10, 0))*10 if tdle > 1000 else int(round(tdle, 0))
        tpm = round(tpm, 3)
        tpue = round(tpue, 3)
        tple = round(tple, 3)
        target_values = [target_name, Pm, Ple, Pue, Rm, Rle, Rue]

        f = open("table_" + self.config_parameters["mode"] + "_" + self.config_parameters["period"] + ".txt", "a+")
        f.write(target_name + " & $" + str(Pm) + "^{" + str(Pue) + "}_{" + str(Ple) + "}$ & $" + str(Rm) + "^{" + str(Rue) + "}_{" + str(Rle) + "}$ & $" + str(R_star) + "\pm" + str(Rse) + "$ & $" + str(tdm) + "^{" + str(tdue) + "}_{" + str(tdle) + "}$ & $" + str(tpm) + "^{" + str(tpue) + "}_{" + str(tple) + "}$ \\\\\n")
        f.close()

        return Pk, P, PP, per, Rk, R, PR, ik, il, Pinc, deltas, ratios, tdm, tdle, tdue, tpm, tpue, tple, target_values



    def epos_pers(self, p0, per, rad, P, M_star):
        """Generates probability of periods using dimensionless spacing in period ratios from EPOS (Mulders et al. 2018)"""

        fP = np.zeros(len(P))
        fD = np.zeros(len(P))
        ind = 0

        for i in range(len(P)):
            if P[i] < p0:
                fP[i] = ((P[i]/12)**1.6 if P[i] < 12 else (P[i]/12)**-0.9)

            else:
                ind = i
                break

        logD = -0.9
        sigma = 0.41

        PRgrid = np.logspace(0,1)

        with np.errstate(divide='ignore'):
            Dgrid = np.log(2.*(PRgrid**(2./3.)-1.)/(PRgrid**(2./3.)+1.))

        Dgrid[0] = -4
        pdfP = spst.norm(logD, sigma).pdf(Dgrid)
        cdfP = spst.norm(logD, sigma).cdf(Dgrid)

        for i in range(ind, len(fP)):
            for j in range(len(per)):
                if P[i] > per[j]:
                    fP[i] = np.interp(P[i]/per[j], PRgrid, pdfP)

                    if j < len(per) - 1 and P[i] < per[j+1]:
                        fP[i] *= np.interp(per[j+1]/P[i], PRgrid, pdfP)

                    elif j < len(per) - 1 and P[i] >= per[j+1]:
                        fP[i] = 0
                    
        m = np.zeros(len(per))

        for k in range(len(per)):
            if self.config_parameters["mass_radius"] == "mrexo":
                m[k] = pfm(measurement=rad[k], predict='mass', dataset='kepler')[0]

            elif self.config_parameters["mass_radius"] == "otegi":
                m[k] = self.otegi_mr(rad[k], "mass")

        for i in range(len(P)):
            for k in range(len(per)):
                m1 = m[k]
                m2 = np.mean(m)
                dc = 8
                a1 = self.K3(P[i], M_star) if P[i] < per[k] else self.K3(per[k], M_star)
                a2 = self.K3(per[k], M_star) if P[i] < per[k] else self.K3(P[i], M_star)
                fD[i] = 2*(a2 - a1)/(a2 + a1) * ((m1 + m2)*const.M_earth.cgs.value/(3*M_star*const.M_sun.cgs.value))**(-1/3)

                if fD[i] < dc:
                    fP[i] = 0

        Du = np.arange(0, max(fD) + 1)
        fDu = np.zeros(len(Du))

        for i in range(len(fD)):
            j = int(fD[i])
            fDu[j] += fP[i]/np.trapz(fP, P)

        return fP/np.trapz(fP, P), fDu



    def syssim_pers(self, per, rad, P, M_star):
        """Generates probability of periods using clustered periods from He, Ford, and Ragozzine (2019)"""

        sigmap = 0.2
        Np = len(per)
        PR = [per[i]/per[i-1] for i in range(1,len(per))]
        bools = []

        for i in range(len(PR)):
            if PR[i] > spst.lognorm.ppf(0.95, Np*sigmap):
                bools.append(1)

            else:
                bools.append(0)

        def calc_bools_dict_value(bl):
            kt_planets = [p for p in range(1,len(bl) + 2)]
            kt_cluster = {}
            kt_cluster_count = 0

            for i in range(0, len(bl)):
                if bl[i] == 0:
                    if kt_planets[i] not in [v for sl in kt_cluster.values() for v in sl]:
                        kt_cluster_count +=  1
                        kt_cluster[kt_cluster_count] = [kt_planets[i], kt_planets[i+1]]

                    else:
                        found_in_key = [ k for (k,v) in kt_cluster.items() if kt_planets[i] in v][0]
                        kt_cluster[found_in_key].append(kt_planets[i+1])

                else:
                    if kt_planets[i] not in [v for sl in kt_cluster.values() for v in sl]:
                        kt_cluster_count += 1
                        kt_cluster[kt_cluster_count] = [kt_planets[i]]
                        kt_cluster_count += 1
                        kt_cluster[kt_cluster_count] = [kt_planets[i+1]]

                    else:
                        kt_cluster_count += 1
                        kt_cluster[kt_cluster_count] = [kt_planets[i+1]]

            return len(kt_cluster), list(kt_cluster.values())

        Nc, Np = calc_bools_dict_value(bools)

        if Nc == 1:
            Pcb = 0
            Pcbs = []
            fPcb = []

            for i in range(len(P)):
                test = 0

                for k in range(len(per)):
                    test += spst.lognorm.pdf(per[k]/P[i], len(per)*sigmap)

                Pcbs.append(P[i])
                fPcb.append(test)

            Pcb = round(Pcbs[np.where(fPcb == max(fPcb))[0][0]], 1)

        else:
            Pcb = []

            for i in range(Nc):
                Pcbs = []
                fPcb = []

                for k in range(len(P)):
                    test = 0

                    for j in Np[i]:
                        test += spst.lognorm.pdf(per[j - 1]/P[k], len(Np[i])*sigmap)

                    Pcbs.append(P[k])
                    fPcb.append(test)

                Pcb.append(Pcbs[np.where(fPcb == max(fPcb))[0][0]])

            Pcb = np.array([round(Pcb[i], 1) for i in range(len(Pcb))])

        Pip = np.linspace(math.sqrt(min(P)/max(P)), math.sqrt(max(P)/min(P)), 1001)
        fPip = []

        for i in range(Nc):
            fPip.append(spst.lognorm.pdf(Pip, (len(per) if Nc == 1 else len(Np[i]))*sigmap))

        fP = np.zeros(len(P))
        fD = np.zeros(len(P))

        for i in range(len(fP)):
            f = []

            for j in range(Nc):
                f.append(np.interp(P[i], (Pcb if Nc == 1 else Pcb[j])*Pip, fPip[j]))

            fP[i] = max(f)

        m = np.zeros(len(per))

        for k in range(len(per)):
            if self.config_parameters["mass_radius"] == "mrexo":
                m[k] = pfm(measurement=rad[k], predict='mass', dataset='kepler')[0]

            elif self.config_parameters["mass_radius"] == "otegi":
                m[k] = self.otegi_mr(rad[k], "mass")

        for i in range(len(P)):
            for k in range(len(per)):
                m1 = m[k]
                m2 = np.mean(m)
                dc = 8
                a1 = self.K3(P[i], M_star) if P[i] < per[k] else self.K3(per[k], M_star)
                a2 = self.K3(per[k], M_star) if P[i] < per[k] else self.K3(P[i], M_star)
                fD[i] = 2*(a2 - a1)/(a2 + a1) * ((m1 + m2)*const.M_earth.cgs.value/(3*M_star*const.M_sun.cgs.value))**(-1/3)

                if fD[i] < dc:
                    fP[i] = 0

        Du = np.arange(0, max(fD) + 1)
        fDu = np.zeros(len(Du))

        for i in range(len(fD)):
            j = int(fD[i])
            fDu[j] += fP[i]/np.trapz(fP, P)

        return fP/np.trapz(fP, P), fDu



    def epos_rads(self, r1, r2):
        """Generates probability of radius from uniform distribution between minimum and maximum planet radii in system (similar to same-R from Mulders et al. 2018)"""

        R = np.arange(r1, r2+0.01, 0.01)
        fR = np.ones(len(R))/(r2 - r1)
        cdfR = spst.uniform(r1, r2-r1).cdf(R)
    
        return R, fR/np.trapz(fR, R), cdfR



    def syssim_rads(self, Rcp, rad):
        """Generates probability of radius using powerlaw or clustered radii from He, Ford, and Ragozzine (2019)"""

        Rmin = float(self.config_parameters["radmin"])
        Rmax = float(self.config_parameters["radmax"])
        Rbreak = 3
        R = np.arange(Rmin, Rmax + 0.01, 0.01)
        pdfR = np.zeros(len(R))
        cdfR = np.zeros(len(R))

        if self.config_parameters["radtype"] == "powerlaw":
            R1 = -1
            R2 = -5

            for i in range(len(R)):
                pdfR[i] = (R[i]/Rbreak)**R1 if R[i] < Rbreak else (R[i]/Rbreak)**R2

            pdfR = pdfR/np.trapz(pdfR, R)

        elif Rcp == "clustered":
            sigmaR = 0.3
            pairwise = [p for p in itertools.combinations(np.arange(len(rad)), 2)]
            allR = []

            for i in range(len(rad) - 1):
                for j in range(i + 1, len(rad)):
                    if rad[j] > spst.lognorm.ppf(0.95, sigmaR, scale=np.exp(np.log(rad[i] + sigmaR**2))) or rad[i] > spst.lognorm.ppf(0.95, sigmaR, scale=np.exp(np.log(rad[j] + sigmaR**2))):
                        allR.append((i, j))

            def reachability_values(reachability_dict):
                return [v for pl in reachability_dict.values() for v in pl]

            def calc_rad_clusters(planet_numbers, pairwise, allR):
                sl_pairwise = [pw for pw in pairwise]

                for sl_pair in allR:
                    sl_pairwise.remove(sl_pair)

                # compute reachability

                reachability_dict = {}
                reachability_count = 0

                for node in planet_numbers:
                    changed = False

                    if node not in reachability_values(reachability_dict):
                        reachability_count +=  1
                        reachability_dict[reachability_count] = [node]
                        changed = True

                    while changed:
                        changed = False     

                        for sl_pair in sl_pairwise:
                            if sl_pair[0] in reachability_dict[reachability_count] and sl_pair[1] not in reachability_values(reachability_dict):
                                reachability_dict[reachability_count].append(sl_pair[1]) 
                                changed = True

                            if sl_pair[1] in reachability_dict[reachability_count] and sl_pair[0] not in reachability_values(reachability_dict):
                                reachability_dict[reachability_count].append(sl_pair[0])
                                changed = True

                return (len(reachability_dict), [ sorted([p+1 for p in cl]) for cl in reachability_dict.values()])

            Nc, Np = calc_rad_clusters(np.arange(len(rad)), pairwise, allR)

            if Nc == 1:
                Rcbs = []
                fRcb = []

                for i in range(len(R)):
                    test = 0

                    for j in range(len(rad)):
                        test += spst.lognorm.pdf(rad[j], sigmaR, scale=np.exp(R[i]))

                    Rcbs.append(R[i])
                    fRcb.append(test)

                Rcb = round(Rcbs[np.where(fRcb == max(fRcb))[0][0]], 2)

            else:
                Rcb = []

                for i in range(Nc):
                    Rcbs = []
                    fRcb = []

                    for k in range(len(R)):
                        test = 0

                        for j in Np[i]:
                            test += spst.lognorm.pdf(rad[j - 1], sigmaR, scale=np.exp(R[k]))

                        Rcbs.append(R[k])
                        fRcb.append(test)

                    Rcb.append(Rcbs[np.where(fRcb == max(fRcb))[0][0]])

                Rcb = np.array([round(Rcb[i], 2) for i in range(len(Rcb))])

            fRi = []

            for i in range(Nc):
                fRi.append(spst.lognorm.pdf(R, sigmaR, scale=np.exp(Rcb) if Nc == 1 else np.exp(Rcb[i])))

            for i in range(len(R)):
                pdfR[i] = max([fRi[j][i] for j in range(len(fRi))])

        for i in range(len(R)):
            cdfR[i] = np.trapz(pdfR[:i + 1], R[:i + 1])/np.trapz(pdfR, R)

        return R, pdfR/np.trapz(pdfR, R), cdfR



    def otegi_mr(self, measurement, predict):
        """Uses a power-law MR to predict mass in Earth values from radius in Earth values or vice versa"""

        if predict == "radius":
            R_r = 1.03*measurement**0.29
            R_v = 0.7*measurement**0.63
            
            if measurement > 25 or (measurement > 5 and self.config_parameters["otegi_rho"] == "volatile"):
                return R_v

            return R_r

        elif predict == "mass":
            M_r = 0.9*measurement**3.45
            M_v = 1.74*measurement**1.58

            if M_r > 25 or (M_r > 5 and self.config_parameters["otegi_rho"] == "volatile"):
                return M_v

            return M_r



    def K3(self, P, M):
        """Calculates semi-major axis in cm using period in days and mass in solar masses"""

        seconds_per_day = 86400
        
        return (const.G.cgs.value*(M*const.M_sun.cgs.value)*(P*seconds_per_day)**2/(4*math.pi**2))**(1/3)



if __name__ == '__main__':
    if len(sys.argv) > 1:
        dynamite(cfname=sys.argv[1])

    else:
        dynamite()
