### DYNAmical Multi-planet Injection TEster (DYNAMITE) ###
### Main File ###
### Jeremy Dietrich ###
### jdietrich1@email.arizona.edu ###
### 2020 October 28 ###
### Version 1.4 ###
### Dietrich & Apai (2020), AJ, 160, 107D ###
### https://iopscience.iop.org/article/10.3847/1538-3881/aba61d ###

import os
import ast
import sys
import math
import rebound
import itertools
import numpy as np
from PPR import PPR
import scipy.stats as spst
from datetime import datetime
import matplotlib.pyplot as plt
import astropy.constants as const
from dynamite_plots import dynamite_plots
from dynamite_targets import dynamite_targets
from mrexo import predict_from_measurement as pfm

class dynamite:

    def __init__(self, merged_data=None):
        """Runs the script"""

        cfname = "dynamite_config.txt"
        self.node_number = 1

        if len(sys.argv) >= 2:
            cfname = sys.argv[1]

        if len(sys.argv) == 3:
            self.node_number = int(sys.argv[2])

        self.seconds_per_day = 86400
        self.G = const.G.cgs.value
        self.au = const.au.cgs.value
        self.M_sun = const.M_sun.cgs.value
        self.R_sun = const.R_sun.cgs.value
        self.M_earth = const.M_earth.cgs.value
        self.R_earth = const.R_earth.cgs.value

        self.config_parameters = {}

        try:
            config_data = np.loadtxt(cfname, dtype=str, delimiter='::')

        except IOError:
            print("Error, configuration file not found!")
            exit()

        for i in range(len(config_data)):
            self.config_parameters[config_data[i, 0]] = config_data[i, 1] if config_data[i, 1].find("[") == -1 else ast.literal_eval(config_data[i, 1])

        if merged_data != None:
            nn = len(merged_data[-3])/4
            self.write_bf_pred_file(merged_data[0], merged_data[4], merged_data[7], merged_data[-2][0:int(len(merged_data[-2])/nn)], merged_data[-1][0:int(len(merged_data[-1])/nn)], merged_data[-4][0], merged_data[-3][0:4], True)
            return

        print(datetime.now(), "Initiating DYNAMITE")

        try:
            self.num_of_nodes = int(os.environ.get('SLURM_JOB_NUM_NODES'))

        except:
            self.num_of_nodes = 1
 
        interations = int(self.config_parameters["MC_chain"])
        self.interation_list = [interations // self.num_of_nodes + (1 if x < interations % self.num_of_nodes else 0) for x in range(self.num_of_nodes)]
        self.seed_start = 0

        for node in range(self.node_number - 1):
            self.seed_start += self.interation_list[node]

        targets_dict = dynamite_targets(self.config_parameters).get_targets(self.config_parameters["mode"], self.config_parameters["system"], self.config_parameters["radmax"], self.config_parameters["removed"])

        if len(targets_dict) == 0:
            print("Error: No targets selected!")
            exit()

        try:
            processes = int(os.environ.get('SLURM_CPUS_PER_TASK'))

        except:
            processes = None

        self.ppr = PPR((self, None), processes)

        if self.config_parameters["saved"] == "False":
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

                    if self.config_parameters["add_unconfirmed"] == "True":
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

                if self.num_of_nodes == 1:
                    np.savez("saved_data.npz", data=data)

                else:
                    np.savez("saved_data_" + str(self.node_number) + ".npz", data=data)

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

                Pk, P, PP, per, Rk, R, PR, ik, il, Pin, deltas, ratios, tdm, tdle, tdue, tpm, tple, tpue, targets, starvs, pers, rads = datavars = self.ppr.create_processes("mt_mc", (targets_dict, targlist), -len(targlist), self.process_data)

                if self.num_of_nodes == 1:
                    np.savez("saved_data.npz", data=datavars)

                else:
                    np.savez("saved_data_" + str(self.node_number) + ".npz", data=datavars)

        elif self.config_parameters["saved"] == "True":
            with np.load("saved_data.npz", allow_pickle=True) as data:
                Pk, P, PP, per, Rk, R, PR, ik, il, Pinc, deltas, ratios, tdm, tdle, tdue, tpm, tple, tpue, targets, starvs, pers, rads = data["data"]

        if self.config_parameters["plot"] == "True":
            if self.config_parameters["saved"] == "False" and self.config_parameters["mode"] == "single":
                Pk, P, PP, per, Rk, R, PR, ik, il, Pinc, deltas, ratios, tdm, tdle, tdue, tpm, tple, tpue, targets, starvs, pers, rads = data

            elif self.config_parameters["saved"] == "False" and self.config_parameters["mode"] != "single":
                Pk, P, PP, per, Rk, R, PR, ik, il, Pinc, deltas, ratios, tdm, tdle, tdue, tpm, tple, tpue, targets, starvs, pers, rads = datavars

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

            return round(np.arccos(planet_pars[0]/(((self.G*star_pars[2]*self.M_sun/(4*math*pi**2))**(1/3)*(planet_pars[1]*self.seconds_per_day)**(2/3))/(star_pars[0]*self.R_sun)))*180/math.pi, 3)
           
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
        ib = ibs[np.where(fib == max(fib))[0][0]]

        if ib > 90:
            ilb = ib
            ib = 180 - ib

        else:
            ilb = 180 - ib

        inew = np.linspace(0, ilb, int(ilb*10) + 1)
        rylghi = spst.rayleigh.pdf(inew, 0, rylgh)
        finew = [rylghi[np.where(np.isclose(inew, abs(ib - il[j])))[0][0]] for j in range(len(il) - 1)]
        finew.append(0.0)

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

        for j in range(len(finew)):
            fi[j] += finew[j]

        fi = fi/np.trapz(fi, il)
        cdfi = np.cumsum(fi)*(il[1]-il[0])
        Pinc = fi

        print(datetime.now(), "Creating Planet Radius Distributions for", target_name)

        if self.config_parameters["radius"] == "epos":
            R, PR, cdfR = self.epos_rads(r1, r2)

        elif self.config_parameters["radius"] == "syssim":
            R, PR, cdfR = self.syssim_rads(self.config_parameters["radtype"], rad)

        print(datetime.now(), "Creating Period Distributions for", target_name)

        if self.config_parameters["period"] == "epos":
            PP, deltas, cdfP = self.epos_pers(p0, per, rad, P, M_star)

        elif self.config_parameters["period"] == "syssim":
            PP, deltas, cdfP = self.syssim_pers(per, rad, P, M_star)

        print(datetime.now(), "Running Monte Carlo for", target_name)
        Pk = np.zeros(int(self.interation_list[self.node_number - 1]))
        Rk = np.zeros(len(Pk))
        ek = np.zeros(len(Pk))
        ik = np.zeros(len(Pk))
        ecc = np.zeros(len(per))
        GMfp213 = (self.G*M_star*self.M_sun/(4*math.pi**2))**(1/3)
        star_values = [R_star, Rse, M_star, Mse]

        Pk, Rk, ek, ik = self.ppr.create_processes("mc_test", (P, R, inew, ib, cdfP, cdfR, cdfi, per, rad, inc, ecc, GMfp213, M_star), -len(Pk), self.process_mc_data)

        tdm, tdle, tdue, tpm, tpue, tple, target_values = self.write_bf_pred_file(Pk, Rk, ik, per, rad, target_name, star_values, (True if self.num_of_nodes == 1 else False))

        return Pk, P, PP, per, Rk, R, PR, ik, il, Pinc, deltas, ratios, tdm, tdle, tdue, tpm, tpue, tple, target_values, star_values


    def write_bf_pred_file(self, Pk, Rk, ik, per, rad, target_name, star_values, write):
        """Writes out the best-fit prediction file."""

        R_star, Rse, M_star, Mse = star_values
        print(datetime.now(), "Calculating Best Fit Predictions for", target_name)
        Pis = np.hstack(np.array([np.arange(0.5, 1.001, 0.001), np.arange(1.01, 10.01, 0.01), np.arange(10.1, 100.1, 0.1), np.arange(101,731,1)]))

        if self.config_parameters["period"] == "epos":
            PPi, _, _ = self.epos_pers(min(per), per, rad, Pis, M_star)

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
        tdm = (Rm*self.R_earth/(R_star*self.R_sun))**2*1e6
        tdle = 2*(Rm*self.R_earth/(R_star*self.R_sun))**2*math.sqrt((Rle/Rm)**2 + (Rse/R_star)**2)*1e6
        tdue = 2*(Rm*self.R_earth/(R_star*self.R_sun))**2*math.sqrt((Rue/Rm)**2 + (Rse/R_star)**2)*1e6
        ntrans = 0
        ntl = 0
        ntu = 0

        print(datetime.now(), "Running Transit Calculations for", target_name)

        for j in range(len(ik)):
            if math.cos(ik[j]*math.pi/180) < (R_star*self.R_sun + Rm*self.R_earth)/self.K3(Pm, M_star):
                ntrans += 1

            if math.cos(ik[j]*math.pi/180) < ((R_star - Rse)*self.R_sun + (Rm - Rle)*self.R_earth)/self.K3(Pue, (M_star + Mse)):
                ntl += 1

            if math.cos(ik[j]*math.pi/180) < ((R_star + Rse)*self.R_sun + (Rm + Rue)*self.R_earth)/self.K3(Ple, (M_star - Mse)):
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

        if write:
            f = open("table_" + self.config_parameters["mode"] + "_" + self.config_parameters["period"] + ".txt", "w")
            f.write(target_name + " & $" + str(Pm) + "^{" + str(Pue) + "}_{" + str(Ple) + "}$ & $" + str(Rm) + "^{" + str(Rue) + "}_{" + str(Rle) + "}$ & $" + str(R_star) + "\pm" + str(Rse) + "$ & $" + str(tdm) + "^{" + str(tdue) + "}_{" + str(tdle) + "}$ & $" + str(tpm) + "^{" + str(tpue) + "}_{" + str(tple) + "}$ \\\\\n")
            f.close()

        return tdm, tdle, tdue, tpm, tpue, tple, target_values



    def mc_test(self, P, R, inew, ib, cdfP, cdfR, cdfi, per, rad, inc, ecc, GMfp213, M_star, k):
        """Runs MC in multi-threading."""

        np.random.seed(self.seed_start + k)
        add_iter = False
        Pmc = P[np.where(np.random.rand() - cdfP < 0)[0][0]]
        Rmc = R[np.where(np.random.rand() - cdfR < 0)[0][0]]
        emc = 0
        ii = np.random.rand()
        iso = np.random.rand()

        if (len(inc) == 2 and iso > 0.38) or (len(inc) == 3 and iso > 0.19) or (len(inc) >= 4):
            if len(np.where(ii - cdfi < 0)[0]) == 0:
                imc = inew[-1] + ib

            else:
                imc = inew[np.where(ii - cdfi < 0)[0][0]]

        else:
            imc = np.arccos(np.random.rand()*2 - 1)*180/math.pi

        permc = per[0]
        radmc = rad[0]

        for p in range(len(per) - 1):
            if Pmc > per[p]*np.sqrt(per[p+1]/per[p]):
                permc = per[p+1]
                radmc = rad[p+1]

        if self.config_parameters["mass_radius"] == "mrexo":
            Mmc = pfm(measurement=Rmc, predict='mass', dataset='kepler')[0]
            masmc = pfm(measurement=radmc, predict='mass', dataset='kepler')[0]

        elif self.config_parameters["mass_radius"] == "otegi": 
            Mmc = self.otegi_mr(Rmc, "mass")
            masmc = self.otegi_mr(radmc, "mass")
            
        a1 = GMfp213*(Pmc*self.seconds_per_day if Pmc < permc else permc*self.seconds_per_day)**(2/3)
        a2 = GMfp213*(permc*self.seconds_per_day if Pmc < permc else Pmc*self.seconds_per_day)**(2/3)
        D = 2*(a2 - a1)/(a2 + a1) * ((Mmc + masmc)*self.M_earth/(3*M_star*self.M_sun))**(-1/3)

        if D >= 8:
            add_iter = True

        else:
            add_iter = False
            """
            stable = True
            sim = rebound.Simulation()
            sim.units = ('day', 'au', 'Msun')
            sim.integrator = "mercurius"
            sim.add(m=M_star)
            m = np.zeros(len(per) + 1)
            
            for p in range(len(per)):
                if self.config_parameters["mass_radius"] == "mrexo":
                    m[p] = pfm(measurement=rad[p], predict='mass', dataset='kepler')[0]

                elif self.config_parameters["mass_radius"] == "otegi": 
                    m[p] = self.otegi_mr(rad[p], "mass")

                sim.add(m=m[p]*self.M_earth/self.M_sun, a=GMfp213*(per[p]*self.seconds_per_day)**(2/3)/self.au, e=ecc[p], inc=inc[p]*math.pi/180)

                if permc == per[p]:
                    mind = p
            
            m[-1] = Mmc
            sim.add(m=Mmc*self.M_earth/self.M_sun, a=GMfp213*(Pmc*self.seconds_per_day)**(2/3)/self.au, e=emc, inc=imc)
            sim.dt = per[0]/40
            amd = np.zeros((len(per) + 1, 3000))

            for it in range(3000):
                sim.integrate(per[0]*5000*(it+1)/3)
                l = sim.calculate_orbits()

                for j in range(1, len(l)):
                    try:
                        amd[j, it] = M_star*(m[j]*self.M_earth/self.M_sun)/(M_star + m[j]*self.M_earth/self.M_sun)*self.M_sun*math.sqrt(self.G*self.M_sun*(m[j]*self.M_earth/self.M_sun)*l[j].a*self.au)*(1-math.sqrt(1-l[j].e**2)*math.cos(l[j].inc*math.pi/180)) 

                    except ValueError:
                        stable = False

                        if j < mind:
                            ps = "PLANET " + str(j)

                        elif j > mind:
                            ps = "PLANET " + str(j-1)

                        else:
                            ps = "INJECTED PLANET"

                        print("SYSTEM ITERATION " + str(self.seed_start + k) + ":", ps, "EJECTED - EXITING INTEGRATION")
                        break

                for j in range(1, len(l) - 1):
                    if l[j].a*(1+l[j].e) >= l[j+1].a*(1-l[j+1].e):
                        stable = False

                        if j < mind:
                            ps1 = "PLANET " + str(j)
                            ps2 = ("PLANET " + str(j+1)) if j + 1 < mind else ("INJECTED PLANET")

                        elif j > mind:
                            ps1 = "PLANET " + str(j-1)
                            ps2 = "PLANET " + str(j)

                        else:
                            ps1 = "INJECTED PLANET"
                            ps2 = "PLANET " + str(j-1)

                        print("SYSTEM ITERATION " + str(self.seed_start + k) + ":", ps1, "AND", ps2, "CROSSED ORBITS - EXITING INTEGRATION")
                        break

                if stable == False:
                    break

            if stable:
                for j in range(len(amd)):
                    amdps = abs(np.fft.fft([amd[j,it] for it in range(len(amd[j]))]))**2
                    mps = max(amdps)

                    if sum([1 if amdps[it] > 0.05*mps else 0 for it in range(3000)]) > 0.01*len(amdps):
                        stable = False

                        if j < mind:
                            ps = "PLANET " + str(j+1)

                        elif j > mind:
                            ps = "PLANET " + str(j)

                        else:
                            ps = "INJECTED PLANET"

                        print("SYSTEM ITERATION " + str(k) + ":", ps, "UNSTABLE VIA SPECTRAL FRACTION")

            if stable:
                add_iter = True
            """
        if add_iter:
            return {k:(Pmc, Rmc, emc, imc)}

        else:
            return {k:None}


    def process_mc_data(self, data):
        """Processes data for the multithreading component"""

        return tuple([list(z) for z in zip(*itertools.chain([data[k] for k in data if data[k] != None]))])


    def epos_pers(self, p0, per, rad, P, M_star):
        """Generates probability of periods using dimensionless spacing in period ratios from EPOS (Mulders et al. 2018)"""

        fP = np.zeros(len(P))
        cdfP = np.zeros(len(P))
        fD = np.zeros(len(P))
        ind = 0
        logD = -0.9
        sigma = 0.41
        PRgrid = np.logspace(0,1)

        with np.errstate(divide='ignore'):
            Dgrid = np.log(2.*(PRgrid**(2./3.)-1.)/(PRgrid**(2./3.)+1.))

        Dgrid[0] = -4
        pdfPR = spst.norm(logD, sigma).pdf(Dgrid)
        cdfPR = spst.norm(logD, sigma).cdf(Dgrid)
        j = 0

        for i in range(len(P)):
            if P[i] < p0:
                fP[i] = ((P[i]/12)**1.6 if P[i] < 12 else (P[i]/12)**-0.9)*np.interp(p0/P[i], PRgrid, pdfPR)

            else:
                ind = i
                break

        for i in range(ind, len(fP)):
            if j < len(per) - 1 and P[i] > per[j+1]:
                j += 1

            fP[i] = np.interp(P[i]/per[j], PRgrid, pdfPR)

            if j != len(per) - 1:
                fP[i] *= np.interp(per[j+1]/P[i], PRgrid, pdfPR)
                    
        m = np.zeros(len(per))

        for k in range(len(per)):
            if self.config_parameters["mass_radius"] == "mrexo":
                m[k] = pfm(measurement=rad[k], predict='mass', dataset='kepler')[0]

            elif self.config_parameters["mass_radius"] == "otegi": 
                m[k] = self.otegi_mr(rad[k], "mass")

        m2 = np.mean(m)
        GMfp213 = (self.G*M_star*self.M_sun/(4*math.pi**2))**(1/3)
        dc = 8
        rats = [math.sqrt(per[k+1]/per[k]) for k in range(len(per) - 1)]
        k = 0

        for i in range(len(P)):
            if k < len(per) - 1 and P[i] > per[k]*rats[k] and P[i] < per[k+1]:
                k += 1

            a1 = GMfp213*(P[i]*self.seconds_per_day if P[i] < per[k] else per[k]*self.seconds_per_day)**(2/3)
            a2 = GMfp213*(per[k]*self.seconds_per_day if P[i] < per[k] else P[i]*self.seconds_per_day)**(2/3)
            fD[i] = 2*(a2 - a1)/(a2 + a1) * ((m[k] + m2)*self.M_earth/(3*M_star*self.M_sun))**(-1/3)
            fP[i] *= spst.norm.cdf(fD[i], loc=8)

        Du = np.arange(0, max(fD) + 1)
        fDu = np.zeros(len(Du))
        trap = np.trapz(fP, P)

        for i in range(len(fD)):
            j = int(fD[i])
            fDu[j] += fP[i]/trap

        for i in range(len(cdfP)):
            cdfP[i] = np.trapz(fP[:i+1], P[:i+1])/trap

        return fP/trap, fDu, cdfP



    def syssim_pers(self, per, rad, P, RR, M_star):
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

        m2 = np.mean(m)
        GMfp213 = (self.G*M_star*self.M_sun/(4*math.pi**2))**(1/3)
        Me3M13 = (self.M_earth/(3*M_star*self.M_sun))**(-1/3)
        dc = 8
        rats = [math.sqrt(per[k+1]/per[k]) for k in range(len(per) - 1)]
        k = 0

        for i in range(len(P)):
            if k < len(per) - 1 and P[i] > per[k]*rats[k] and P[i] < per[k+1]:
                k += 1

            a1 = GMfp213*(P[i]*self.seconds_per_day if P[i] < per[k] else per[k]*self.seconds_per_day)**(2/3)
            a2 = GMfp213*(per[k]*self.seconds_per_day if P[i] < per[k] else P[i]*self.seconds_per_day)**(2/3)
            fD[i] = 2*(a2 - a1)/(a2 + a1) * (m[k] + m2)*Me3M13

            #if fD[i] < dc:
                #fP[i] = 0

        Du = np.arange(0, max(fD) + 1)
        fDu = np.zeros(len(Du))
        trap = np.trapz(fP, P)

        for i in range(len(fD)):
            j = int(fD[i])
            fDu[j] += fP[i]/trap

        for i in range(len(fP)):
            cdfP[i] = np.trapz(fP[:i + 1], P[:i + 1])/trap

        return fP/trap, fDu



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
        
        return (self.G*M*self.M_sun*(P*self.seconds_per_day)**2/(4*math.pi**2))**(1/3)



if __name__ == '__main__':
    dynamite()
