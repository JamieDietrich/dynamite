### DYNAmical Multi-planet Injection TEster (DYNAMITE) ###
### Main File ###
### Jamie Dietrich ###
### jdietrich@asu.edu ###
### 2024 October 21 ###
### Version 3.0 ###
### Dietrich & Apai (2020), AJ, 160, 107D ###
### https://iopscience.iop.org/article/10.3847/1538-3881/aba61d ###
### Dietrich & Apai (2021), AJ, 161, 17D ###
### Dietrich, Apai, & Malhotra (2022), AJ, 163, 88D ###
### Basant, Dietrich, & Apai (2022), AJ, 164, 12B ###
### Basant, Dietrich, & Apai (2022), RNAAS, 6, 213 ###
### Dietrich (2024), AJ, 168, 119 ###

import os
import ast
import sys
import math
import socket
import rebound
import itertools
import numpy as np
import scipy.stats as spst
import multiprocessing as mp
import matplotlib.pyplot as plt

from PPR import PPR
from datetime import datetime
from astropy import constants as const
from scipy.signal import argrelextrema
from dynamite_plots import dynamite_plots
from dynamite_targets_db import dynamite_targets_db

class dynamite:

    def __init__(self, merged_data=None):
        """Runs the script"""

        self.cfname = "dynamite_config.txt"
        self.node_number = 1

        if len(sys.argv) >= 2:
            self.cfname = sys.argv[1]

        if len(sys.argv) == 3:
            self.node_number = int(sys.argv[2])

        self.seconds_per_day = 86400
        self.G = const.G.cgs.value
        self.au = const.au.cgs.value
        self.M_sun = const.M_sun.cgs.value
        self.R_sun = const.R_sun.cgs.value
        self.M_earth = const.M_earth.cgs.value
        self.R_earth = const.R_earth.cgs.value
        self.MES = self.M_earth/self.M_sun
        self.config_parameters = {}
        self.sdir = ""

        try:
            config_data = np.loadtxt(self.cfname, dtype=str, delimiter=':')

        except IOError:
            print("Error: configuration file not found!")
            exit()

        for i in range(len(config_data)):
            self.config_parameters[config_data[i, 0]] = config_data[i, 1] if config_data[i, 1].find("[") == -1 else ast.literal_eval(config_data[i, 1])

        self.startdatetime = str(datetime.now()).replace(" ", "_")
        print(self.startdatetime, "Initiating DYNAMITE on node", socket.gethostname())

        try:
            self.node_number = int(os.environ.get('SLURM_ARRAY_TASK_ID'))
            self.num_of_nodes = int(os.environ.get('SLURM_ARRAY_TASK_COUNT'))

        except:
            self.node_number = 1
            self.num_of_nodes = 1
 
        self.interations = int(self.config_parameters["MC_chain"])
        self.interation_list = [self.interations // self.num_of_nodes + (1 if x < self.interations % self.num_of_nodes else 0) for x in range(self.num_of_nodes)]
        self.seed_start = 0

        for node in range(self.node_number - 1):
            self.seed_start += self.interation_list[node]

        try:
            processes = int(os.environ.get('SLURM_CPUS_PER_TASK'))

        except:
            processes = None

        if sys.platform != "darwin" and os.cpu_count() > 16:
            self.ppr = PPR((self, None), processes)

        if merged_data != None:
            try:
                if len(merged_data[-5][0]) > 1:
                    for i in range(len(merged_data[-5])):
                        self.write_bf_pred_file((self.G*merged_data[-4][i][2]*self.M_sun/(4*math.pi**2))**(1/3), merged_data[0][i], merged_data[3][i], merged_data[6][i], merged_data[9][i], merged_data[-3][i], merged_data[-2][i], merged_data[-1][i], merged_data[11][i], merged_data[-5][i][0], merged_data[-4][i], True)

            except:
                self.write_bf_pred_file((self.G*merged_data[-4][2]*self.M_sun/(4*math.pi**2))**(1/3), merged_data[0], merged_data[3], merged_data[6], merged_data[9], merged_data[-3], merged_data[-2], merged_data[-1], merged_data[11], merged_data[-5][0], merged_data[-4], True)

            with np.load(self.sdir + "saved_data_" + self.config_parameters["system"].replace(" ", "_") if self.config_parameters["mode"] == "single" else self.config_parameters["mode"] + ".npz", allow_pickle=True) as data:
                Pk, P, PP, Rk, R, PR, ik, il, Pinc, ek, el, Pecc, deltas, tdm, tdle, tdue, tpm, tple, tpue, targets, starvs, pers, rads, mass, eccs = [data["arr_" + str(i)] for i in range(len(data))]

            if self.config_parameters["plot"] == "True":
                self.run_plots(data)

            return

        if self.config_parameters["saved"] == "False":
            if self.config_parameters["mode"] == "single":
                if isinstance(self.config_parameters["system"], str):
                    self.sdir = self.config_parameters["system"].replace("\"", "") + "/"
                    
                    if self.sdir != "" and not os.path.exists(self.sdir):
                        os.makedirs(self.sdir)
                        
                    data = self.run_single(self.config_parameters["system"].replace("\"", ""))

                    if self.config_parameters["plot"] == "True":
                        self.run_plots(data, self.config_parameters["system"].replace("\"", ""))

                elif isinstance(self.config_parameters["system"], list):
                    with open(self.sdir + "table_" + self.config_parameters["mode"] + "_" + self.config_parameters["period"] + "_run_" + self.startdatetime + ".txt", "w") as f:
                        f.write("Name --- Period --- Radius")
                        f.write(" --- Transit probability --- RV Limit")
                        f.write("\n")

                    for i in self.config_parameters["system"]:
                        self.sdir = self.config_parameters["system"][i].replace("\"", "") + "/"
                        
                        if self.sdir != "" and not os.path.exists(self.sdir):
                            os.makedirs(self.sdir)
                        
                        data = self.run_single(i)

                        if self.config_parameters["plot"] == "True":
                            self.run_plots(data, i)

            else:
                targets_dict = dynamite_targets_db(self.config_parameters["targets_db"]).get_targets(self.config_parameters["mode"], self.config_parameters["system"], self.config_parameters["radmax"], self.config_parameters["removed"])
                limits_dict = dynamite_targets_db(self.config_parameters["targets_db"]).get_limits(self.config_parameters["mode"], self.config_parameters["system"], self.config_parameters["radmax"])

                if targets_dict == None or len(targets_dict) == 0:
                    print("Error: No targets selected!")
                    exit()
                targlist = [tn for tn in targets_dict.keys()]

                if sys.platform != "darwin" and os.cpu_count() > 16:
                    Pk, P, PP, Rk, R, PR, ik, il, Pinc, ek, el, Pecc, deltas, tdm, tdle, tdue, tpm, tple, tpue, targets, starvs, pers, rads, mass, eccs = datavars = self.ppr.create_processes("mt_mc", (targets_dict, limits_dict, targlist), -len(targlist), self.process_data)

                else:
                    Pk, P, PP, Rk, R, PR, ik, il, Pinc, ek, el, Pecc, deltas, tdm, tdle, tdue, tpm, tple, tpue, targets, starvs, pers, rads, mass, eccs = datavars = self.run_new_mp(self.mt_mc, np.arange(len(targlist)), (targets_dict, limits_dict, targlist))

                if self.num_of_nodes == 1:
                    np.savez_compressed(self.sdir + "saved_data_" + self.config_parameters["mode"] + ".npz", *datavars)

                else:
                    np.savez_compressed(self.sdir + "saved_data_" + self.config_parameters["mode"] + "_" + str(self.node_number) + ".npz", *datavars)

                if self.config_parameters["plot"] == "True":
                    self.run_plots(datavars)

        elif self.config_parameters["saved"] == "True":
            if self.config_parameters["mode"] == "single":
                if isinstance(self.config_parameters["system"], str):
                    self.sdir = self.config_parameters["system"].replace("\"", "") + "/"

                    if self.sdir != "" and not os.path.exists(self.sdir):
                        os.makedirs(self.sdir)
                        
                    datavars = []

                    with np.load(self.sdir + "saved_data_" + self.config_parameters["system"].replace(" ", "_").replace("\"", "") + ".npz", allow_pickle=True) as data:
                        for i in range(len(data)):
                            datavars.append(data["arr_" + str(i)])

                    self.saved_writing(datavars, self.config_parameters["system"])

                    if sys.platform != "darwin" and os.cpu_count() > 16:
                        self.ppr.terminate_PPR()

                elif isinstance(self.config_parameters["system"], list):
                    with open(self.sdir + "table_" + self.config_parameters["mode"] + "_" + self.config_parameters["period"] + "_run_" + self.startdatetime + ".txt", "w") as f:
                        f.write("Name --- Period --- Radius")
                        f.write(" --- Transit probability --- RV Limit")
                        f.write("\n")

                    for i in self.config_parameters["system"]:
                        self.sdir = self.config_parameters["system"][i].replace("\"", "") + "/"

                        if self.sdir != "" and not os.path.exists(self.sdir):
                            os.makedirs(self.sdir)
                            
                        datavars = []

                        with np.load(self.sdir + "saved_data_" + i.replace(" ", "_") + ".npz", allow_pickle=True) as data:
                            for i in range(len(data)):
                                datavars.append(data["arr_" + str(i)])

                        self.saved_writing(datavars, i)

        print(datetime.now(), "Finishing DYNAMITE")

        if sys.platform != "darwin" and os.cpu_count() > 16:
            self.ppr.terminate_PPR()



    def saved_writing(self, datavars, targsys):
        """Writes the files and runs the plots from saved data."""

        Pk, P, PP, Rk, R, PR, ik, il, Pinc, ek, el, Pecc, deltas, tdm, tdle, tdue, tpm, tple, tpue, targets, starvs, pers, rads, mass, eccs = datavars

        try:
            if len(targets[0]) > 1:
                for i in range(len(targets)):
                    self.write_bf_pred_file((self.G*starvs[i][2]*self.M_sun/(4*math.pi**2))**(1/3), Pk[i], Rk[i], ik[i], ek[i], pers[i], rads[i], mass[i], eccs[i], Pecc[i], targets[i][0], starvs[i], True)

        except:
            self.write_bf_pred_file((self.G*starvs[2]*self.M_sun/(4*math.pi**2))**(1/3), Pk, Rk, ik, ek, pers, rads, mass, eccs, Pecc, targets[0], starvs, True)

        if self.config_parameters["plot"] == "True":
            self.run_plots(datavars, targsys)



    def run_single(self, targsys):
        """Runs DYNAMITE in single mode for the given target name."""

        targets_dict = dynamite_targets_db(self.config_parameters["targets_db"]).get_targets(self.config_parameters["mode"], targsys, self.config_parameters["radmax"], self.config_parameters["removed"])
        limits_dict = dynamite_targets_db(self.config_parameters["targets_db"]).get_limits(self.config_parameters["mode"], targsys, self.config_parameters["radmax"])

        if targets_dict == None or len(targets_dict) == 0:
            print("Error: No targets selected!")
            exit()

        alt_name, R_star, Rse, M_star, Mse, target, target_name, limits = self.set_up(targets_dict, limits_dict, targsys)
        target = [i[:-1] for i in target]
        targ_rem = []

        for i in range(len(target) - 1 if len(self.config_parameters["removed"]) > 0 else len(target)):
            if target[i][0][0] not in self.config_parameters["removed"]:
                targ_rem.append(target[i])

        if len(self.config_parameters["additional"][0]) > 0:
            for i in self.config_parameters["additional"]:
                x = []

                for j in range(len(i) - 1):
                    if (j == 2 and self.config_parameters["use_mass"] == "True") or (j == 1 and i[1][0] == "?"):
                        mv = self.mr_convert(i[2][0], "radius")
                        x.append([mv, self.mr_convert(i[2][0]+i[2][1], "radius") - mv, mv - self.mr_convert(i[2][0]-i[2][2], "radius")])

                    elif j == 2 and i[2][0] == "?":
                        mv = self.mr_convert(i[1][0], "mass")
                        x.append([mv, self.mr_convert(i[1][0]+i[1][1], "mass") - mv, mv - self.mr_convert(i[1][0]-i[1][2], "mass"), "Prediction"])

                    else:
                        x.append(i[j])

                targ_rem.append(x)

        if self.config_parameters["add_unconfirmed"] == "True":
            for i in self.config_parameters["unconfirmed"]:
                x = []

                for j in range(len(i) - 1):
                    #if isinstance(i[j], tuple):
                        #if self.config_parameters["use_mass"] == "True":
                            #x.append(self.mr_convert(i[j][0], "radius"))

                        #else:
                            #x.append(i[j][0])

                    if (j == 2 and self.config_parameters["use_mass"] == "True") or (j == 1 and i[j] == "?"):
                        x.append(self.mr_convert(i[2], "radius"))

                    elif j == 2 and i[j] == "?":
                        x.append(self.mr_convert(i[1], "mass"))

                    else:
                        x.append(i[j])

                targ_rem.append(x)

        target = sorted(targ_rem, key=lambda x: x[0][0])
        data = self.run_monte_carlo(R_star, Rse, M_star, Mse, target, target_name, limits) + ([target[i][0] for i in range(len(target))], [target[i][1] for i in range(len(target))], [target[i][2] for i in range(len(target))], [target[i][4] for i in range(len(target))])

        if self.num_of_nodes == 1:
            np.savez_compressed(self.sdir + "saved_data_" + target_name.replace(" ", "_") + ".npz", *data)

        else:
            np.savez_compressed(self.sdir + "saved_data_" + target_name.replace(" ", "_") + "_" + str(self.node_number) + ".npz", *data)

        return data



    def run_plots(self, data, targsys):
        """Runs the plotting program."""

        Pk, P, PP, Rk, R, PR, ik, il, Pinc, ek, el, Pecc, deltas, tdm, tdle, tdue, tpm, tple, tpue, targets, starvs, pers, rads, mass, eccs = data
        print(datetime.now(), "Creating Plots")
        dynamite_plots(Pk, P, PP, Rk, R, PR, ik, il, Pinc, ek, el, Pecc, deltas, tdm, tdle, tdue, tpm, tple, tpue, targets, starvs, pers, rads, mass, eccs, targsys, self.cfname, (self.ppr if (sys.platform != "darwin" and os.cpu_count() > 16) else None))



    def mt_mc(self, targets_dict, limits_dict, targlist, i):
        """Runs the Monte Carlo code on multiple threads"""

        alt_name, R_star, Rse, M_star, Mse, target, target_name, limits = self.set_up(targets_dict, limits_dict, targlist[i])

        try:
            inds = np.argsort([target[i][0][0] for i in range(len(target))])
            target = [target[inds[i]] for i in range(len(inds))]

        except IndexError:
            print("ERROR IN DICT ON TARGET", target_name)

        data = self.run_monte_carlo(R_star, Rse, M_star, Mse, target, target_name, limits) + ([target[i][0] for i in range(len(target))], [target[i][1] for i in range(len(target))], [target[i][2] for i in range(len(target))], [target[i][4] for i in range(len(target))])

        return data
    


    def process_data(self, data):
        """Processes data for the multithreading component"""

        return tuple([list(z) for z in zip(*itertools.chain(data))])



    def set_up(self, targets_dict, limits_dict, target):
        """Sets up target"""

        t = list(targets_dict[target])
        limits = list(limits_dict[target])
        pnum = len([x[0] for x in t if isinstance(x[0], tuple)])

        if t[0][2] == None:
            t[0][2] = 0.0

        if t[0][4] == None:
            t[0][4] = 0.0

        for x in range(1, pnum + 1):
            if t[x][1][0] == "?":
                if t[x][2][3] == "Msini" and t[x][3][0] != "?" and t[x][3][0] != 0:
                    val = t[x][2][0]/np.sin(t[x][3][0]*np.pi/180)
                    upp = (t[x][2][0] + t[x][2][1])/np.sin((t[x][3][0] + t[x][3][2])*np.pi/180) - val if t[x][2][1] != "?" else "?"
                    low = (t[x][2][0] + t[x][2][2])/np.sin((t[x][3][0] + t[x][3][1])*np.pi/180) - val if t[x][2][1] != "?" else "?"
                    t[x][2] = (val, upp, low, "Msini/sin(i)")

                val = self.mr_convert(t[x][2][0], "radius")
                fupp, flow = self.get_mr_force(t[x][2][0])
                upp = self.mr_convert(t[x][2][0] + t[x][2][1], "radius", fupp) - val if t[x][2][1] != "?" else "?"
                low = self.mr_convert(t[x][2][0] + t[x][2][2], "radius", flow) - val if t[x][2][2] != "?" else "?"
                t[x][1] = (val, upp, low)

            elif t[x][2][0] == "?":
                val = self.mr_convert(t[x][1][0], "mass")
                fupp, flow = self.get_mr_force(val)
                upp = self.mr_convert(t[x][1][0] + t[x][1][1], "mass", fupp) - val if t[x][1][1] != "?" else "?"
                low = self.mr_convert(t[x][1][0] + t[x][1][2], "mass", flow) - val if t[x][1][2] != "?" else "?"
                t[x][2] = (val, upp, low, "Prediction")

        return t[0][0], t[0][1], t[0][2], t[0][3], t[0][4], list([t[i][:-1] for i in range(1, len(t))]), target, limits



    def get_mr_force(self, val):
        """Gets the force value for upper and lower limits"""

        flow = None
        fupp = None

        if self.config_parameters["mass_radius"] == "otegi":
            if self.config_parameters["otegi_rho"] == "volatile":
                if val > 5:
                    flow = "volatile"

                else:
                    fupp = "rocky"

            elif self.config_parameters["otegi_rho"] == "rocky":
                if val > 25:
                    flow = "volatile"

                else:
                    fupp = "rocky"

        return fupp, flow



    def mr_convert(self, meas, pred, force=None):
        """Runs conversion from mass to radius and vice versa"""

        if self.config_parameters["mass_radius"] == "otegi":
            return self.otegi_mr(meas, pred, force)



    def run_monte_carlo(self, R_star, Rse, M_star, Mse, target, target_name, limits):
        """Runs the Monte Carlo analysis."""

        GMfp213 = (self.G*M_star*self.M_sun/(4*math.pi**2))**(1/3)
        per = [target[i][0][0] for i in range(len(target))]
        rad = [target[i][1][0] for i in range(len(target))]
        mas = [target[i][2][0] for i in range(len(target))]
        mts = [target[i][2][3] for i in range(len(target))]
        inc = [target[i][3][0] for i in range(len(target))]
        ecc = [target[i][4][0] for i in range(len(target))]

        if target_name == "Inner Solar System":
            P = np.arange(0.5, 4331.5, 0.1)

        elif target_name == "eps Eri":
            P = np.arange(100, 4331.5, 0.1)

        elif target_name == "Proxima Centauri":
            P = np.arange(0.5, 1412, 0.1)

        else:
            P = np.arange(0.5, 730.5, 0.1)

        print(datetime.now(), "Creating Inclination Prior Distribution for", target_name)

        if self.config_parameters["inclination"] == "rayleigh_iso":
            il, Pinc, cdfi, inew, ib = self.rayleigh_iso_incs(inc, per, GMfp213, R_star)

        elif self.config_parameters["inclination"] == "syssim":
            il, Pinc, cdfi, inew, ib = self.syssim_incs(inc, per, GMfp213, R_star)

        print(datetime.now(), "Creating Eccentricity Prior Distribution for", target_name)

        if self.config_parameters["eccentricity"] == "rayleigh":
            el, Pecc, cdfe = self.rayleigh_eccs()

        elif self.config_parameters["eccentricity"] == "syssim":
            el, Pecc, cdfe = self.syssim_eccs(len(ecc))

        em = np.median(Pecc)

        Pk = []
        Rk = []
        ek = []
        ik = []
        Pk_r = []
        Rk_r = []
        ek_r = []
        ik_r = []
        r_reason = []
        star_values = [R_star, Rse, M_star, Mse]

        per_tup = [target[i][0] for i in range(len(target))]
        rad_tup = [target[i][1] for i in range(len(target))]
        mas_tup = [target[i][2] for i in range(len(target))]
        inc_tup = [target[i][3] for i in range(len(target))]
        ecc_tup = [target[i][4] for i in range(len(target))]
        accept = "NO"
        itr = 0

        while accept != "YES" and itr < 100:
            inct = [self.random_val_from_tup(i) if i[0] != "?" else 90 for i in inc_tup]
            ecct = [self.random_val_from_tup(i) if i[0] != "?" else 0 for i in ecc_tup]
            pert = [self.random_val_from_tup(i) for i in per_tup]
            mast = [self.random_val_from_tup(mas_tup[i]) if mas_tup[i][3] != "Msini" else (self.random_val_from_tup(mas_tup[i])/np.sin(inct[i]*math.pi/180) if inct[i] != 0 and inct[i] != "?" else self.random_val_from_tup(mas_tup[i])*self.M_sun/self.M_earth) for i in range(len(mas_tup))]
            radt = [self.random_val_from_tup(rad_tup[i]) if mas_tup[i][3] != "Msini" else self.mr_convert(mast[i], "radius") for i in range(len(rad_tup))]
            D = np.zeros(len(pert) - 1)

            for i in range(len(D)):
                a1 = GMfp213*(pert[i]*self.seconds_per_day)**(2/3)
                a2 = GMfp213*(pert[i+1]*self.seconds_per_day)**(2/3)
                e1 = ecct[i]
                e2 = ecct[i+1]
                D[i] = 2*(a2*(1 - e2) - a1*(1 + e1))/(a2 + a1) * ((mast[i] + mast[i+1])*self.M_earth/(3*M_star*self.M_sun))**(-1/3)

            accept, _, _ = self.run_stability_analysis(D, M_star, pert, None, mast, inct, ecct, None, None, None, GMfp213, limits)
            itr += 1

        if itr == 100:
            for i in range(len(ecc_tup)):
                if ecc_tup[i][0] != "?":
                    ecc_tup[i] = (0, ecc_tup[i][1], ecc_tup[i][2])

        indq = False

        for i in range(len(inc_tup)):
            if inc_tup[i][0] == "?" and mas_tup[i][3] == "Msini":
                indq = True

        fP = []

        if not indq and np.all([i[0] != "?" for i in rad_tup]):
            print(datetime.now(), "Creating Radius Prior Distribution for", target_name)

            if self.config_parameters["radius"] == "epos":
                R, PR, cdfR = self.epos_rads(min(rad), max(rad))

            elif self.config_parameters["radius"] == "syssim":
                R, PR, cdfR = self.syssim_rads(rad)

            print(datetime.now(), "Creating Period Prior Distribution for", target_name)

            if self.config_parameters["period"] == "epos":
                fP = self.create_fP(P, per)
                PP, deltas, cdfP = self.epos_pers(fP, per, mas, [], ecc, P, R, cdfR, il, Pinc, Pecc, em, cdfe, M_star, GMfp213)

            elif self.config_parameters["period"] == "syssim":
                PP, deltas, cdfP = self.syssim_pers(per, mas, [], ecc, P, R, cdfR, il, Pinc, Pecc, em, cdfe, M_star, GMfp213)

        else:
            Ms = np.zeros((len(il), len(mas)))
            Rs = np.zeros((len(il), len(mas)))
            PR = []

            print(datetime.now(), "Creating Radius Prior Distribution for", target_name)

            for n in range(len(il)):
                for m in range(len(mas)):
                    if rad[m] == "?" or (mts[m] == "Msini" and inc[m] == "?"):
                        Ms[n, m] = min(mas[m]/np.sin(il[n]*math.pi/180) if il[n] != 0 and il[n] != 180 else M_star*self.M_sun/self.M_earth, M_star*self.M_sun/self.M_earth)
                        Rs[n, m] = self.mr_convert(Ms[n, m], "radius")

                    elif rad[m] != "?":
                        Rs[n, m] = rad[m]
                        Ms[n, m] = self.mr_convert(rad[m], "mass")

                    elif mts[m] == "Mass":
                        Ms[n, m] = mas[m]
                        Rs[n, m] = self.mr_convert(mas[m], "radius")
                
                if self.config_parameters["radius"] == "epos":
                    R, fR, _ = self.epos_rads(Rs[n])

                elif self.config_parameters["radius"] == "syssim":
                    R, fR, _ = self.syssim_rads(Rs[n])

                PR.append(fR)

            PR = np.average(PR, axis=0, weights=Pinc)
            cdfR = np.cumsum(PR)
            cdfR = cdfR/cdfR[-1]

            print(datetime.now(), "Creating Period Prior Distribution for", target_name)

            if self.config_parameters["period"] == "epos":
                fP = self.create_fP(P, per)
                PP, deltas, cdfP = self.epos_pers(fP, per, mas, Ms, ecc, P, R, cdfR, il, Pinc, Pecc, em, cdfe, M_star, GMfp213)

            elif self.config_parameters["period"] == "syssim":
                PP, deltas, cdfP = self.syssim_pers(per, mas, Ms, ecc, P, R, cdfR, il, Pinc, Pecc, em, cdfe, M_star, GMfp213)

        print(datetime.now(), "Running Monte Carlo for", target_name)

        if sys.platform != "darwin" and os.cpu_count() > 16:
            PPi, PRi, Ri, deltasi, stable, val, tim, _Pk, _Rk, _ek, _ik = self.ppr.create_processes("mc_test", -int(self.interation_list[self.node_number - 1]), (P, fP, PP, cdfP, per, deltas, R, PR, cdfR, inew, ib, il, cdfi, el, Pecc, em, cdfe, per_tup, rad_tup, mas_tup, inc_tup, ecc_tup, GMfp213, R_star, M_star, limits, indq), self.process_mc_data)

        else:
            res = self.run_new_mp(self.mc_test, np.arange(self.interations), (P, fP, PP, cdfP, per, deltas, R, PR, cdfR, inew, ib, il, cdfi, el, Pecc, em, cdfe, per_tup, rad_tup, mas_tup, inc_tup, ecc_tup, GMfp213, R_star, M_star, limits, indq))
            results = []

            for j in range(len(res[0])):
                res1 = []

                for i in range(len(res)):
                    res1.append(res[i][j])

                results.append(res1)

            PPi, PRi, Ri, deltasi, stable, val, tim, _Pk, _Rk, _ek, _ik = results

        R = Ri[0]
        PP = np.mean(PPi, axis=0)
        PR = np.mean(PRi, axis=0)
        deltas = np.mean(deltasi, axis=0)

        for a in range(len(stable)):
            if stable[a] == "YES":
                Pk.append(_Pk[a])
                Rk.append(_Rk[a])
                ek.append(_ek[a])
                ik.append(_ik[a])

            else:
                r_reason.append(stable[a])
                Pk_r.append(_Pk[a])
                Rk_r.append(_Rk[a])
                ek_r.append(_ek[a])
                ik_r.append(_ik[a])

        print(target_name, "Accepted Injections:", len(Pk))

        if self.config_parameters["stability"] == "spock":
            print("STABILITY THRESHOLD VALUE:", np.mean(val), "+/-", np.std(val), "TIME TO INSTABILITY", np.mean(tim)/per[0], "+/-", np.std(tim)/per[0])

        elif self.config_parameters["stability"] == "specfrac":
            np.savetxt(self.sdir + "unstable_times.txt", tim)

        with open(self.sdir + "rejected_values_" + str(self.node_number) + ".txt", "w") as f:
            for a in range(len(Pk_r)):
                f.write(r_reason[a] + "\t" + str(round(Pk_r[a], 1)) + "\t" + str(round(Rk_r[a], 2)) + "\t" + str(round(ek_r[a], 3)) + "\t" + str(round(ik_r[a], 1)) + "\n")

        tdm, tdle, tdue, tpm, tpue, tple, target_values = self.write_bf_pred_file(GMfp213, Pk, Rk, ik, ek, per, rad, mas, ecc, Pecc, target_name, star_values, (True if self.num_of_nodes == 1 else False))

        return Pk, P, PP, Rk, R, PR, ik, il, Pinc, ek, el, Pecc, deltas, tdm, tdle, tdue, tpm, tpue, tple, target_values, star_values



    def write_bf_pred_file(self, GMfp213, Pk, Rk, ik, ek, per, rad, mas, ecc, Pecc, target_name, star_values, write):
        """Writes out the best-fit prediction file."""

        print(datetime.now(), "Calculating Best Fit Predictions for", target_name)
        R_star, Rse, M_star, Mse = star_values

        per_vals = []
        rad_vals = []
        mas_vals = []
        ecc_vals = []

        for i in range(len(per)):
            per_vals.append(per[i][0] if isinstance(per[i], tuple) else per[i])
            rad_vals.append(rad[i][0] if isinstance(rad[i], tuple) else rad[i])
            mas_vals.append(mas[i][0] if isinstance(mas[i], tuple) else mas[i])
            ecc_vals.append(ecc[i][0] if isinstance(ecc[i], tuple) else ecc[i])

        Pis = np.hstack([np.logspace(-0.3, 0.5, 18), np.logspace(0.51, 2.9, 121)])
        bins = np.append(Pis, 10**2.91)
        widths = [bins[i+1] - bins[i] for i in range(len(bins)-1)]
        PPis, _ = np.histogram(Pk, bins=bins)
        PPi = PPis/widths
        zl = [l for l in range(len(PPi)) if PPi[l] == 0]
        peaks = []

        for i in range(1, len(zl)):
            if zl[i-1] != zl[i] - 1:
                peaks.append([zl[i-1]+1, zl[i]])

        like_sums = [np.sum(PPi[peaks[s][0]:peaks[s][1]]) for s in range(len(peaks))]
        pvs = [max(PPi[i[0]:i[1]]) for i in peaks]
        Pfs = [Pis[peaks[i][0]:peaks[i][1]][np.where(PPi[peaks[i][0]:peaks[i][1]] == pvs[i])[0][0]] for i in range(len(peaks))]
        Pms = []

        for s in range(len(Pfs)):
            ind = np.where(Pis == Pfs[s])[0][0]

            if PPi[ind] > 0.5*np.amax(PPi) and (like_sums[s] > 0.1*sum(like_sums) or like_sums[s] > 0.25*max(like_sums)):
                Pms.append(Pfs[s])
        
        if len(Pms) == 0:
            for s in range(len(Pfs)):
                if PPi[np.where(Pis == Pfs[s])[0][0]] > 0.5*np.amax(PPi):
                    Pms.append(Pfs[s])

        Pmes = np.zeros(len(Pms))
        Ples = np.zeros(len(Pms))
        Pues = np.zeros(len(Pms))
        
        for pm in range(len(Pms)):
            a = []

            try:
                x = np.where((PPi < np.amax(PPi)*0.01) & (Pis < Pms[pm]))[0][-1]
                y = np.where((PPi < np.amax(PPi)*0.01) & (Pis > Pms[pm]))[0][0]
                z = np.where((Pk > Pis[x]) & (Pk < Pis[y]))[0]

                for aa in z:
                    a.append(Pk[aa])

            except:
                pass

            if len(a) > 0:
                Pmes[pm] = np.percentile(a, 50)
                Ples[pm] = Pmes[pm] - np.percentile(a, 16)
                Pues[pm] = np.percentile(a, 84) - Pmes[pm]

            else:
                Pmes[pm] = Pms[pm]
                Ples[pm] = np.where((PPi < np.exp(-1)*PPi[np.where(Pis == Pms[pm])[0][0]]) & (Pis < Pms[pm]))[0][-1]

                try:
                    Pues[pm] = np.where((PPi < np.exp(-1)*PPi[np.where(Pis == Pms[pm])[0][0]]) & (Pis > Pms[pm]))[0][0]

                except:
                    Pues[pm] = max(Pis)

            if pm == len(Pms) - 1:
                a = []

                for aa in range(len(Pk)):
                    tol = 1 if Pk[aa] < 100 else 10

                    try:
                        if Pk[aa] > Pms[pm]/2 and Pk[aa] < Pms[pm]*4 and PPi[np.where(np.isclose(Pis, Pk[aa], atol=tol))[0][0]] > np.amax(PPi)*0.2:
                            a.append(Pk[aa])

                    except IndexError:
                        pass

                if len(a) > 0:
                    Pmes[pm] = np.percentile(a, 50)
                    Ples[pm] = Pmes[pm] - np.percentile(a, 16)
                    Pues[pm] = np.percentile(a, 84) - Pmes[pm]

                else:
                    Pmes[pm] = Pms[pm]
                    Ples[pm] = np.where((PPi < np.exp(-1)*PPi[np.where(Pis == Pms[pm])[0][0]]) & (Pis < Pms[pm]))[0][-1]
                    Pues[pm] = np.where((PPi < np.exp(-1)*PPi[np.where(Pis == Pms[pm])[0][0]]) & (Pis > Pms[pm]))[0][0]

        Pmmax = Pis[np.where(PPi == np.amax(PPi))[0][0]]

        try:
            Plemax = Ples[np.where(Pms == Pmmax)[0][0]]

        except:
            Plemax = Ples[0]

        try:
            Puemax = Pues[np.where(Pms == Pmmax)[0][0]]

        except:
            Puemax = Pues[0]

        Rm = np.percentile(Rk, 50)
        Rle = Rm - np.percentile(Rk, 16)
        Rue = np.percentile(Rk, 84) - Rm
        Mm = self.mr_convert(Rm, "mass")
        fupp, flow = self.get_mr_force(Mm)
        Mle = Mm - self.mr_convert(Rm - Rle, "mass", flow)
        Mue = self.mr_convert(Rm + Rue, "mass", fupp) - Mm
        im = np.percentile(ik, 50)
        ile = im - np.percentile(ik, 16)
        iue = np.percentile(ik, 84) - im
        em = np.percentile(ek, 50)
        ele = em - np.percentile(ek, 16)
        eue = np.percentile(ek, 84) - em
        tdm = (Rm*self.R_earth/(R_star*self.R_sun))**2*1e6
        tdle = 2*(Rm*self.R_earth/(R_star*self.R_sun))**2*math.sqrt((Rle/Rm)**2 + (Rse/R_star)**2)*1e6
        tdue = 2*(Rm*self.R_earth/(R_star*self.R_sun))**2*math.sqrt((Rue/Rm)**2 + (Rse/R_star)**2)*1e6

        print(datetime.now(), "Running Transit Calculations for", target_name)
        tpm = np.zeros(len(Pms))
        tpue = np.zeros(len(Pms))
        tple = np.zeros(len(Pms))
        rvsa = np.zeros(len(Pms))
        rvue = np.zeros(len(Pms))
        rvle = np.zeros(len(Pms))

        for pm in range(len(Pms)):
            ntrans = 0
            ntl = 0
            ntu = 0

            for j in range(len(ik)):
                if math.cos(ik[j]*math.pi/180) < (R_star*self.R_sun + Rm*self.R_earth)/self.K3(Pms[pm], M_star):
                    ntrans += 1

                if math.cos(ik[j]*math.pi/180) < ((R_star - Rse)*self.R_sun + (Rm - Rle)*self.R_earth)/self.K3(Pues[pm], (M_star + Mse)):
                    ntl += 1

                if math.cos(ik[j]*math.pi/180) < ((R_star + Rse)*self.R_sun + (Rm + Rue)*self.R_earth)/self.K3(Ples[pm], (M_star - Mse)):
                    ntu += 1

            tpm[pm] = ntrans/len(ik)
            tple[pm] = max(1e-3, (ntrans - ntl)/len(ik)) if (tpm[pm] != 0 and tpm[pm] != 1) else (ntrans - ntl)/len(ik)
            tpue[pm] = max(1e-3, (ntu - ntrans)/len(ik)) if (tpm[pm] != 0 and tpm[pm] != 1) else (ntu - ntrans)/len(ik)
            rvsa[pm] = 0.6395*Pms[pm]**(-1/3)*Mm*math.sin(im*math.pi/180)*(M_star+Mm*self.MES)**(-2/3)*(1-em**2)**-0.5
            rvle[pm] = rvsa[pm] - (0.6395*(Pms[pm]-Ples[pm])**(-1/3)*(Mm-Mle)*math.sin((im-ile)*math.pi/180)*(M_star+(Mm-Mle)*self.MES)**(-2/3)*(1-(em-ele)**2)**-0.5)
            rvue[pm] = (0.6395*(Pms[pm]+Pues[pm])**(-1/3)*(Mm+Mue)*math.sin((im+iue)*math.pi/180)*(M_star+(Mm+Mue)*self.MES)**(-2/3)*(1-(em+eue)**2)**-0.5) - rvsa[pm]

        print(datetime.now(), "Writing out Best Values for", target_name)

        for pm in range(len(Pms)):
            Pms[pm] = round(Pms[pm], (1 if Pms[pm] > 10 else 3 if Pms[pm] < 1 else 2))
            Pmes[pm] = round(Pmes[pm], (1 if Pmes[pm] > 10 else 3 if Pmes[pm] < 1 else 2))
            Pues[pm] = round(Pues[pm], (1 if Pues[pm] > 10 else 3 if Pues[pm] < 1 else 2))
            Ples[pm] = round(Ples[pm], (1 if Ples[pm] > 10 else 3 if Ples[pm] < 1 else 2))
            tpm[pm] = round(tpm[pm], 3)
            tpue[pm] = round(tpue[pm], 3)
            tple[pm] = round(tple[pm], 3)
            rvsa[pm] = round(rvsa[pm], 2)
            rvue[pm] = round(rvue[pm], 2)
            rvle[pm] = round(rvle[pm], 2)

        Rm = round(Rm, (3 if Rm < 1 else 2))
        Rue = round(Rue, (3 if Rue < 1 else 2))
        Rle = round(Rle, (3 if Rle < 1 else 2))
        Mm = round(Mm, (3 if Mm < 1 else 2))
        Mue = round(Mue, (3 if Mue < 1 else 2))
        Mle = round(Mle, (3 if Mle < 1 else 2))
        im = round(im, 1)
        iue = round(iue, 1)
        ile = round(ile, 1)
        em = round(em, 3)
        eue = round(eue, 3)
        ele = round(ele, 3)
        tdm = int(round(tdm/10, 0))*10 if tdm > 1000 else int(round(tdm, 0))
        tdue = int(round(tdue/100, 0))*100 if tdue > 10000 else int(round(tdue/10, 0))*10 if tdue > 1000 else int(round(tdue, 0))
        tdle = int(round(tdle/10, 0))*10 if tdle > 1000 else int(round(tdle, 0))
        target_values = [target_name, Pmmax, Plemax, Puemax, Rm, Rle, Rue]

        if write:
            with open(self.sdir + "table_" + self.config_parameters["mode"] + "_" + self.config_parameters["period"] + "_" + target_name + ".txt", "w") as f:
                f.write("Name & Period (d) & Planet Radius (R_Earth) & Mass (M_Earth) & Inclination (deg) & Eccentricity & Stellar Radius (R_Sun) & Transit Depth (ppm) & Transit Probability & RV Semi-amplitude (m/s)\\\\\n" + target_name + " & ($")

                for pm in range(len(Pms)):
                    f.write(str(Pms[pm]) + "(" + str(Pmes[pm]) + ")^{" + str(Pues[pm]) + "}_{" + str(Ples[pm]) + ("}$, $" if pm != len(Pms) - 1 else "}$) & $"))

                f.write(str(Rm) + "^{" + str(Rue) + "}_{" + str(Rle) + "}$ & $" + str(Mm) + "^{" + str(Mue) + "}_{" + str(Mle) + "}$ & $" + str(im) + "^{" + str(iue) + "}_{" + str(ile) + "}$ & $" + str(em) + "^{" + str(eue) + "}_{" + str(ele) + "}$ & $" + str(R_star) + "\\pm" + str(Rse) + "$ & $" + str(tdm) + "^{" + str(tdue) + "}_{" + str(tdle) + "}$ & ($")

                for pm in range(len(tpm)):
                    f.write(str(tpm[pm]) + "^{" + str(tpue[pm]) + "}_{" + str(tple[pm]) + ("}$, $" if pm != len(tpm) - 1 else "}$) & ($"))

                for pm in range(len(rvsa)):
                    f.write(str(rvsa[pm]) + "^{" + str(rvue[pm]) + "}_{" + str(rvle[pm]) + ("}$, $" if pm != len(rvsa) - 1 else "}$)\\\\\n"))

                f.write("\nName --- Period --- Radius --- Transit probability --- RV Limit\n" + target_name + " --- ")

                for pm in range(len(Pms)):
                    f.write(str(Pms[pm]) + " (+" + str(Pues[pm]) + "/-" + str(Ples[pm]) + ("), " if pm != len(Pms) - 1 else ") --- "))

                f.write(str(Rm) + "(+" + str(Rue) + "/-" + str(Rle) + ") --- ")

                for pm in range(len(tpm)):
                    f.write(str(tpm[pm]) + " (+" + str(tpue[pm]) + "/-" + str(tple[pm]) + ("), " if pm != len(tpm) - 1 else ") --- "))

                for pm in range(len(rvsa)):
                    f.write(str(rvsa[pm]) + " (+" + str(rvue[pm]) + "/-" + str(rvle[pm]) + ("), " if pm != len(rvsa) - 1 else ")\n"))

            if isinstance(self.config_parameters["system"], list):
                with open(self.sdir + "table_" + self.config_parameters["mode"] + "_" + self.config_parameters["period"] + "_run_" + self.startdatetime + ".txt", "a") as f:
                    f.write(target_name + " --- ")

                    for pm in range(len(Pms)):
                        f.write(str(Pms[pm]) + " (+" + str(Pues[pm]) + "/-" + str(Ples[pm]) + ("), " if pm != len(Pms) - 1 else ")"))

                    f.write(" --- " + str(Rm) + "(+" + str(Rue) + "/-" + str(Rle) + ")")

                    for pm in range(len(tpm)):
                        f.write(" --- " + str(tpm[pm]) + " (+" + str(tpue[pm]) + "/-" + str(tple[pm]) + ("), " if pm != len(tpm) - 1 else ")"))

                    for pm in range(len(rvsa)):
                        f.write(" --- " + str(rvsa[pm]) + " (+" + str(rvue[pm]) + "/-" + str(rvle[pm]) + ("), " if pm != len(rvsa) - 1 else ")"))

                    f.write("\n")

        mtp = max(tpm)
        mtpu = tpue[np.where(tpm == mtp)[0][0]]
        mtpl = tple[np.where(tpm == mtp)[0][0]]

        return tdm, tdle, tdue, mtp, mtpu, mtpl, target_values



    def random_val_from_tup(self, tup):
        """Draws a random value from the given tuple (value, upper_unc, lower_unc)"""

        if tup[1] == "?":
            return tup[0]

        if tup[1] == abs(tup[2]):
            return np.random.normal(tup[0], tup[1])

        else:
            chk = np.random.rand()

            if chk < 0.5:
                x = np.random.normal(tup[0], tup[1])
                return (x if x > tup[0] else 2*tup[0]-x)

            else:
                x = np.random.normal(tup[0], abs(tup[2]))
                return (x if x < tup[0] else 2*tup[0]-x)



    def run_new_mp(self, func, arr, mp_args):
        """Runs new multiprocessing code needed for iOS users."""

        with mp.Pool(processes=mp.cpu_count()) as pool:
            args = [(i, mp_args) for i in arr]
            results = pool.starmap(func, args)

        return results



    #def mc_test(self, P, fP, PP, cdfP, per, deltas, R, PR, cdfR, inew, ib, il, cdfi, el, Pecc, em, cdfe, per_tup, rad_tup, mas_tup, inc_tup, ecc_tup, GMfp213, R_star, M_star, limits, indq, k):



    def mc_test(self, k, args):
        """Runs Monte Carlo in multi-threading."""

        P, fP, PP, cdfP, per, deltas, R, PR, cdfR, inew, ib, il, cdfi, el, Pecc, em, cdfe, per_tup, rad_tup, mas_tup, inc_tup, ecc_tup, GMfp213, R_star, M_star, limits, indq = args
        np.random.seed(self.seed_start + k + np.random.randint(100000))
        accept = "NO"
        itr = 0

        while accept != "YES" and itr < 100:
            inc = []
            ecc = []

            for i in range(len(inc_tup)):
                inc.append(self.random_val_from_tup(inc_tup[i]) if inc_tup[i][0] != "?" else "?")
                ecc.append(self.random_val_from_tup(ecc_tup[i]) if ecc_tup[i][0] != "?" else "?")

            emc = self.cdf_draw(el, np.random.rand(), cdfe)

            for j in range(len(ecc)):
                if ecc[j] == "?":
                    ecc[j] = self.cdf_draw(el, np.random.rand(), cdfe)

            if len(cdfi) == 0:
                fiq = [np.sin(il[j]*math.pi/180) for j in range(len(il))]
                cdfiq = np.cumsum(fiq)*(il[1]-il[0])
                ib = self.cdf_draw(il, np.random.rand(), cdfiq)

                if ib > 90:
                    it = round(180 - ib, 1)

                else:
                    it = round(ib, 1)

                if np.array(inc).any() == "?":
                    while it + 4 > np.arccos(R_star*self.R_sun/(GMfp213*(per[inc.index("?")]*self.seconds_per_day)**(2/3))*180/math.pi):
                        ib = self.cdf_draw(il, np.random.rand(), cdfiq)

                        if ib > 90:
                            it = round(180 - ib, 1)

                        else:
                            it = round(ib, 1)

            if self.config_parameters["inclination"] == "rayleigh_iso":
                for j in range(len(inc)):
                    if inc[j] == "?":
                        inc[j] = self.inc_draw(inc, cdfi, inew, ib, il)

                    if inc[j] < 0:
                        inc[j] = -inc[j]

                imc = self.inc_draw(inc, cdfi, inew, ib, il)

                if imc < 0:
                    imc = -imc

            elif self.config_parameters["inclination"] == "syssim":
                for j in range(len(inc)):
                    if inc[j] == "?":
                        inc[j] = self.cdf_draw(il, np.random.rand(), cdfi)

                    if inc[j] < 0:
                        inc[j] = -inc[j]

                imc = self.cdf_draw(il, np.random.rand(), cdfi)

                if imc < 0:
                    imc = -imc

            per = [self.random_val_from_tup(i) for i in per_tup]
            mas = [self.random_val_from_tup(mas_tup[i]) if mas_tup[i][3] != "Msini" else (self.random_val_from_tup(mas_tup[i])/np.sin(inc[i]*math.pi/180) if inc[i] != 0 else M_star*self.M_sun/self.M_earth) for i in range(len(mas_tup))]
            rad = [self.random_val_from_tup(rad_tup[i]) if mas_tup[i][3] != "Msini" else self.mr_convert(mas[i], "radius") for i in range(len(rad_tup))]
            D = np.zeros(len(per) - 1)

            for i in range(len(D)):
                a1 = GMfp213*(per[i]*self.seconds_per_day)**(2/3)
                a2 = GMfp213*(per[i+1]*self.seconds_per_day)**(2/3)
                e1 = ecc[i]
                e2 = ecc[i+1]
                D[i] = 2*(a2*(1 - e2) - a1*(1 + e1))/(a2 + a1) * ((mas[i] + mas[i+1])*self.M_earth/(3*M_star*self.M_sun))**(-1/3)

            accept, _, _ = self.run_stability_analysis(D, M_star, per, None, mas, inc, ecc, None, None, None, GMfp213, limits)
            itr += 1

        Pmc = self.cdf_draw(P, np.random.rand(), cdfP)
        Rmc = self.cdf_draw(R, np.random.rand(), cdfR)
        Mmc = self.mr_convert(Rmc, "mass")

        if Pmc < per[0]:
            pns = [per[0]]
            rns = [rad[0]]
            ens = [ecc[0]]

        elif Pmc > per[-1]:
            pns = [per[-1]]
            rns = [rad[-1]]
            ens = [ecc[-1]]

        else:
            pind = np.where(per < Pmc)[0][-1]
            pns = [per[pind], per[pind+1]]
            rns = [rad[pind], rad[pind+1]]
            ens = [ecc[pind], ecc[pind+1]]

        Dn = np.zeros(len(pns))

        for i in range(len(Dn)):
            permc = pns[i]
            radmc = rns[i]
            eccmc = ens[i]
            masmc = self.mr_convert(radmc, "mass")
            a1 = GMfp213*((Pmc if Pmc < permc else permc)*self.seconds_per_day)**(2/3)
            a2 = GMfp213*((permc if Pmc < permc else Pmc)*self.seconds_per_day)**(2/3)
            e1 = emc if Pmc < permc else eccmc
            e2 = eccmc if Pmc < permc else emc
            Dn[i] = 2*(a2*(1 - e2) - a1*(1 + e1))/(a2 + a1) * ((Mmc + masmc)*self.M_earth/(3*M_star*self.M_sun))**(-1/3)
        
        accept, val, tim = self.run_stability_analysis(Dn, M_star, per, Pmc, mas, inc, ecc, Mmc, imc, emc, GMfp213, limits)
        return PP, PR, R, deltas, accept, val, tim, Pmc, Rmc, emc, imc



    def run_stability_analysis(self, D, M_star, per, Pmc, mas, inc, ecc, Mmc, imc, emc, GMfp213, limits):
        """Runs through the stability analysis based on the criterion utilized (Mutual Hill Radius, SPOCK, or spectral fraction)"""

        stability = self.config_parameters["stability"]
        accept = "YES"
        val = 0
        tim = 0

        if stability == "hill":
            dlim = spst.lognorm.rvs(spst.norm.rvs(0.40, 0.02, 1), loc=0, scale=np.exp(spst.norm.rvs(1.97, 0.03, 1)), size=1)

            if np.any(D < dlim):
                accept = "HILL RADIUS"

        if stability == "spock" or stability == "specfrac" and np.any(D) < 100:
            try:
                sim = rebound.Simulation()

            except:
                stability = "hill"
                print("WARNING: REBOUND not installed on this machine - not performing N body integration dynamical stability analysis")

            if stability == "spock":
                try:
                    from spock import FeatureClassifier, DeepRegressor

                except:
                    stability = "hill"
                    print("WARNING: SPOCK not installed on this machine - not performing N body integration dynamical stability analysis")

        if stability == "spock" or stability == "specfrac" and np.any(D) < 100:
            sim, pf, mf, incf, eccf, mind = self.create_sim_system(M_star, per, Pmc, mas, inc, ecc, Mmc, imc, emc, GMfp213)
            sim.dt = per[0]/40
            amd = np.zeros((len(per) + 1, 3000))
            Ms = M_star*self.M_sun

            if stability == "spock":
                model = DeepRegressor()
                val = model.predict_stable(sim)
                tim, _, _ = model.predict_instability_time(sim)

                if val < 0.34:
                    accept = "SPOCK"

            else:
                for it in range(3000):
                    sim.integrate(per[0]*5000*(it+1)/3)
                    l = sim.calculate_orbits()

                    for j in range(0, len(l)):
                        try:
                            mp = mas[j]*self.M_earth
                            amd[j, it] = Ms*mp/(Ms + mp)*math.sqrt(self.G*(Ms*mp)*l[j].a*self.au)*(1-math.sqrt(1-l[j].e**2)*math.cos(l[j].inc*math.pi/180))

                        except ValueError:
                            accept = "EJECTED"

                            if j < mind:
                                ps = "PLANET " + str(j+1)

                            elif j > mind:
                                ps = "PLANET " + str(j)

                            else:
                                ps = "INJECTED PLANET"

                            print("SYSTEM ITERATION " + str(self.seed_start + k) + ":", ps, "EJECTED - EXITING INTEGRATION AT TIMESTEP", str(it))
                            tim = it
                            break

                    for j in range(0, len(l) - 1):
                        if l[j].a*(1+l[j].e) >= l[j+1].a*(1-l[j+1].e):
                            accept = "ORBITS CROSSED"

                            if j < mind:
                                ps1 = "PLANET " + str(j+1)
                                ps2 = ("PLANET " + str(j+2)) if j + 1 < mind else ("INJECTED PLANET")

                            elif j > mind:
                                ps1 = "PLANET " + str(j)
                                ps2 = "PLANET " + str(j+1)

                            else:
                                ps1 = "INJECTED PLANET"
                                ps2 = "PLANET " + str(j+1)

                            print("SYSTEM ITERATION " + str(self.seed_start + k) + ":", ps1, "AND", ps2, "CROSSED ORBITS - EXITING INTEGRATION AT TIMESTEP", str(it))
                            tim = it
                            break

                    if accept != "YES":
                        break

            if accept == "YES":
                for j in range(len(amd)):
                    amdps = abs(np.fft.fft([amd[j,it] for it in range(len(amd[j]))]))**2
                    mps = max(amdps)

                    if sum([1 if amdps[it] > 0.01*mps else 0 for it in range(3000)]) > 0.01*len(amdps):
                        accept = "SPECTRAL FRACTION"

                        if j < mind:
                            ps = "PLANET " + str(j+1)

                        elif j > mind:
                            ps = "PLANET " + str(j)

                        else:
                            ps = "INJECTED PLANET"

                        print("SYSTEM ITERATION " + str(k) + ":", ps, "UNSTABLE VIA SPECTRAL FRACTION")

        if accept == "YES" and Pmc != None and len(limits) != 0:
            ttv_lim = []

            for l in limits:
                if l[1] == "rv" and Pmc > l[2] and Pmc < l[3] and 0.6395*Pmc**(-1/3)*Mmc*math.sin(imc*math.pi/180)*(M_star+Mmc*self.MES)**(-2/3)*(1-emc**2)**-0.5 > l[4]:
                    accept = "RV LIMIT"

                elif l[1] == "transit" and Pmc > l[2] and Pmc < l[3] and Rmc > l[4] and imc*math.pi/180 < (math.pi/2 - math.atan((R_star*self.R_sun + Rmc*self.R_earth)/(GMfp213*(Pmc*self.seconds_per_day)**(2/3)))):
                    accept = "TRANSIT LIMIT"
                
                elif l[1] == "ttv":
                    ttv_lim.append([l[2], l[4]])

            if len(ttv_lim) > 0:
                good_mode = False
                es = 0
                ttvmode = self.config_parameters["ttv_mode"]

                while not good_mode:
                    if ttvmode == "rebound":
                        try:
                            sim = rebound.Simulation()
                            good_mode = True

                        except:
                            print("REBOUND is not installed on this machine - trying TTVFaster")
                            es += 1
                            ttvmode == "ttvfaster"

                    if ttvmode == "ttvfaster":
                        try:
                            from ttvfaster import run_ttvfaster
                            good_mode = True

                        except:
                            print("TTVFaster is not installed on this machine - trying REBOUND")
                            es += 1
                            ttvmode == "rebound"

                    if es == 2:
                        print("Neither TTVFaster or REBOUND is installed on this machine - not running TTV analysis")
                        break

                if good_mode and ttvmode == "rebound":
                    sim, pf, mf, incf, eccf, mind = self.create_sim_system(M_star, per, Pmc, mas, inc, ecc, Mmc, imc, emc, GMfp213)
                    sim.dt = per[0]/5
                    sim.move_to_com()
                    sim.save("beginning.bin")
                    ps = sim.particles
                    N = 100
                    p = 0

                    for l in ttv_lim:
                        for pp in range(len(per)):
                            if abs(per[pp] - l[0]) < 0.5:
                                p = pp + 1

                        tts = np.zeros(N)
                        nn = 0

                        while nn < N:
                            yo = ps[p].y - ps[0].y
                            to = sim.t
                            sim.integrate(sim.t + sim.dt)
                            tn = sim.t

                            if yo*(ps[p].y - ps[0].y) < 0 and ps[p].x - ps[0].x > 0:
                                while tn - to > 1e-5:
                                    if yo*(ps[p].y - ps[0].y) < 0:
                                        tn = sim.t

                                    else:
                                        to = sim.t

                                    sim.integrate((tn+to)/2)

                                tts[nn] = sim.t
                                nn += 1
                                sim.integrate(sim.t + sim.dt)

                        A = np.vstack([np.ones(N), range(N)]).T
                        c, m = np.linalg.lstsq(A, tts, rcond=-1)[0]
                        ttvs = tts - (m*np.array(range(N))-c)*3600/(2*math.pi)

                        if np.random.rand() > np.percentile(ttvs, l[1]):
                            accept = "TTV LIMIT"
                            break

                        del sim
                        sim = rebound.Simulation("beginning.bin")

                elif good_mode and ttvmode == "ttvfaster":
                    pf, mf, incf, eccf, _ = self.create_sim_arrays(per, Pmc, mas, inc, ecc, Mmc, imc, emc)
                    ttvaa = []
                    fb = False

                    for it in range(1000):
                        ttva = []
                        lana = [np.random.rand()*2*math.pi - math.pi for j in range(len(pf))]
                        apa = [np.random.rand()*2*math.pi - math.pi for j in range(len(pf))]
                        params = [M_star]

                        for j in range(len(pf)):
                            for i in [mf[j]*self.MES, pf[j], eccf[j]*math.cos(apa[j]), incf[j], lana[j], eccf[j]*math.sin(apa[j]), np.random.rand()*pf[j]]:
                                params.append(i)

                        tts = run_ttvfaster(len(pf), params, 0.0, 1000.0, 6)

                        for p in range(len(tts)):
                            for t in range(len(tts[p])):
                                tts[p][t] -= t*pf[p]

                            ttvs = tts[p] - np.mean(tts[p])
                            ttva.append(np.mean([abs(min(ttvs)), max(ttvs)])*1440)

                        ttvaa.append(ttva)

                    for l in ttv_lim:
                        for p in range(len(pf)):
                            if pf[p] == l[0]:
                                if np.random.rand() > np.percentile([i[p] for i in range(len(ttvaa))], l[1]):
                                    accept = "TTV LIMIT"
                                    fb = True
                                    break

                        if fb:
                            break

        return accept, val, tim



    def process_mc_data(self, data):
        """Processes data for the multithreading component"""

        return tuple([list(z) for z in zip(*itertools.chain(data))])



    def cdf_draw(self, grid, val, cdf):
        """Draws from the cdf of a distribution"""

        if len(np.where(val - cdf < 0)[0]) == 0:
            return grid[-1]

        else:
            return grid[np.where(val - cdf < 0)[0][0]]



    def inc_draw(self, inc, cdfi, inew, ib, il):
        """Does the MC draw for the Rayleigh-isotropic joint distribution."""

        ii = np.random.rand()
        iso = np.random.rand()

        if (len(inc) == 2 and iso > 0.38) or (len(inc) == 3 and iso > 0.19) or (len(inc) >= 4):
            if len(np.where(ii - cdfi < 0)[0]) == 0:
                imc = inew[-1] + ib

            else:
                imc = il[np.where(ii - cdfi < 0)[0][0]]

        else:
            imc = np.arccos(np.random.rand()*2 - 1)*180/math.pi

        return imc



    def create_fP(self, P, per):
        """Creates the initial fP distribution for the EPOS period ratio distribution."""

        fP = np.zeros(len(P))
        ind = 0
        PRgrid = np.logspace(0,1)

        with np.errstate(divide='ignore'):
            Dgrid = np.log(2.*(PRgrid**(2./3.)-1.)/(PRgrid**(2./3.)+1.))

        Dgrid[0] = -4
        pdfPR = spst.norm(-0.9, 0.41).pdf(Dgrid)
        j = 0

        for i in range(len(P)):
            if P[i] < min(per):
                fP[i] = ((P[i]/12)**1.6 if P[i] < 12 else (P[i]/12)**-0.9)*np.interp(min(per)/P[i], PRgrid, pdfPR)

            else:
                ind = i
                break

        for i in range(ind, len(fP)):
            if j < len(per) - 1 and P[i] > per[j+1]:
                j += 1

            fP[i] = (np.interp(P[i]/per[j], PRgrid, pdfPR) if (P[i]/per[j] < 5 or j == len(per) - 1) else 1)

            if j != len(per) - 1:
                fP[i] *= (np.interp(per[j+1]/P[i], PRgrid, pdfPR) if per[j+1]/P[i] < 5 else 1)

        return fP



    def epos_pers(self, fP, per, mas, Ms, ecc, P, R, cdfR, il, Pinc, Pecc, em, cdfe, M_star, GMfp213):
        """Generates probability of periods using dimensionless spacing in period ratios from EPOS (Mulders et al. 2018)"""

        m2 = self.mr_convert(R[np.where(cdfR > 0.5)[0][0]], "mass")
        #dc = 8
        l = 0

        Dv = np.zeros(len(P))
        Ds = np.zeros(len(il))
        pdfP = np.zeros(len(il))

        for i in range(len(P)):
            if l < len(per) - 1 and P[i] > per[l]*math.sqrt(per[l+1]/per[l]) and P[i] < per[l+1]:
                l += 1

            a1 = GMfp213*((P[i] if P[i] < per[l] else per[l])*self.seconds_per_day)**(2/3)
            a2 = GMfp213*((per[l] if P[i] < per[l] else P[i])*self.seconds_per_day)**(2/3)
            e1 = em if P[i] < per[l] else (ecc[l] if ecc[l] != "?" else em)
            e2 = (ecc[l] if ecc[l] != "?" else em) if P[i] < per[l] else em
            aet = 2*(a2*(1 - e2) - a1*(1 + e1))/(a2 + a1)

            if len(Ms) == 0:
                Dv[i] = max(0, aet * ((mas[l] + m2)*(self.M_earth/(3*M_star*self.M_sun)))**(-1/3))
                fP[i] *= spst.lognorm.cdf(Dv[i], spst.norm.rvs(0.4, 0.02, 1), loc=0, scale=np.exp(spst.norm.rvs(1.97, 0.03, 1)))
    
            else:
                for n in range(len(il)):
                    Ds[n] = max(0, aet * ((Ms[n, l] + m2)*(self.M_earth/(3*M_star*self.M_sun)))**(-1/3))
                    pdfP[n] = spst.lognorm.cdf(Ds[n], spst.norm.rvs(0.4, 0.02, 1), loc=0, scale=np.exp(spst.norm.rvs(1.97, 0.03, 1)))

                fP[i] *= np.average(pdfP, weights=Pinc)
                Dv[i] = np.average(Ds, weights=Pinc)

        Du = np.arange(0, 1000)
        fDu = np.zeros(len(Du))

        for i in range(len(Dv)):
            j = int(Dv[i])
            fDu[j] += fP[i]

        cdfP = np.cumsum(fP)

        return fP/cdfP[-1], fDu/cdfP[-1], cdfP/cdfP[-1]



    def syssim_pers(self, per, mas, Ms, ecc, P, R, cdfR, il, Pinc, Pecc, em, cdfE, M_star, GMfp213):
        """Generates probability of periods using clustered periods from He, Ford, and Ragozzine (2019)"""

        sigmap = 0.2
        Np = len(per)
        bools = []

        for i in range(len(per) - 1):
            if per[i+1]/per[i] > spst.lognorm.ppf(0.95, Np*sigmap):
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
                        found_in_key = [k for (k,v) in kt_cluster.items() if kt_planets[i] in v][0]
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
            Pcb = np.round(P[np.argmax([sum(spst.lognorm.pdf(per[k]/P[i], len(per)*sigmap) for k in range(len(per))) for i in range(len(P))])], 1)

        else:
            Pcb = np.round([P[np.argmax([sum(spst.lognorm.pdf(per[j-1]/P, len(Np[i])*sigmap) for j in Np[i])])] for i in range(Nc)], 1)

        Pip = np.linspace(math.sqrt(min(P)/max(P)), math.sqrt(max(P)/min(P)), 1001)
        fPip = []

        for i in range(Nc):
            fPip.append(spst.lognorm.pdf(Pip, (len(per) if Nc == 1 else len(Np[i]))*sigmap))

        fP = np.array([max([np.interp(P[i], (Pcb if Nc == 1 else Pcb[j])*Pip, fPip[j]) for j in range(Nc)]) for i in range(len(P))])
        m2 = self.mr_convert(R[np.where(cdfR > 0.5)[0][0]], "mass")
        l = 0

        Dv = np.zeros(len(P))
        Ds = np.zeros(len(il))
        pdfP = np.zeros(len(il))

        for i in range(len(P)):
            if l < len(per) - 1 and P[i] > per[l]*math.sqrt(per[l+1]/per[l]) and P[i] < per[l+1]:
                l += 1

            a1 = GMfp213*((P[i] if P[i] < per[l] else per[l])*self.seconds_per_day)**(2/3)
            a2 = GMfp213*((per[l] if P[i] < per[l] else P[i])*self.seconds_per_day)**(2/3)
            e1 = em if P[i] < per[l] else (ecc[l] if ecc[l] != "?" else em)
            e2 = (ecc[l] if ecc[l] != "?" else em) if P[i] < per[l] else em
            aet = 2*(a2*(1 - e2) - a1*(1 + e1))/(a2 + a1)

            if len(Ms) == 0:
                Dv[i] = max(0, aet * ((mas[l] + m2)*(self.M_earth/(3*M_star*self.M_sun)))**(-1/3))
                fP[i] *= spst.lognorm.cdf(Dv[i], spst.norm.rvs(0.4, 0.02, 1), loc=0, scale=np.exp(spst.norm.rvs(1.97, 0.03, 1)))
    
            else:
                for n in range(len(il)):
                    Ds[n] = max(0, aet * ((Ms[n, l] + m2)*(self.M_earth/(3*M_star*self.M_sun)))**(-1/3))
                    pdfP[n] = spst.lognorm.cdf(Ds[n], spst.norm.rvs(0.4, 0.02, 1), loc=0, scale=np.exp(spst.norm.rvs(1.97, 0.03, 1)))

                fP[i] *= np.average(pdfP, weights=Pinc)
                Dv[i] = np.average(Ds, weights=Pinc)

        Du = np.arange(0, 1000)
        fDu = np.zeros(len(Du))

        for i in range(len(Dv)):
            j = int(Dv[i])
            fDu[j] += fP[i]

        cdfP = np.cumsum(fP)

        return fP/cdfP[-1], fDu/cdfP[-1], np.cumsum(fP)/cdfP[-1]


    def epos_rads(self, r1, r2):
        """Generates probability of radius from uniform distribution between minimum and maximum planet radii in system (similar to same-R from Mulders et al. 2018)"""

        R = np.arange(r1, r2, 0.025)
        pdfR = spst.uniform(r1, r2-r1).pdf(R)
        cdfR = spst.uniform(r1, r2-r1).cdf(R)
    
        return R, pdfR, cdfR



    def syssim_rads(self, rad):
        """Generates probability of radius using powerlaw or clustered radii from He, Ford, and Ragozzine (2019)"""

        R = np.arange(float(self.config_parameters["radmin"]), float(self.config_parameters["radmax"]) + 0.025, 0.025)

        if self.config_parameters["radtype"] == "powerlaw":
            R1 = -1
            R2 = -5
            Rbreak = 3
            pdfR = [(i/Rbreak)**R1 if i < Rbreak else (i/Rbreak)**R2 for i in R]

        elif self.config_parameters["radtype"] == "clustered":
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

                return (len(reachability_dict), [sorted([p+1 for p in cl]) for cl in reachability_dict.values()])

            Nc, Np = calc_rad_clusters(np.arange(len(rad)), pairwise, allR)

            if Nc == 1:
                Rcb = np.round(R[np.argmax([sum(spst.lognorm.pdf(rad[k], sigmaR, scale=np.exp(R[i])) for k in range(len(rad))) for i in range(len(R))])], 1)

            else:
                Rcb = np.round([R[np.argmax([sum(spst.lognorm.pdf(rad[j-1], sigmaR, scale=np.exp(Ri)) for j in Np[i] for Ri in R)])] for i in range(Nc)], 1)

            
            fRi = [spst.lognorm.pdf(R, sigmaR, scale=np.exp(Rcb) if Nc == 1 else np.exp(Rcb[i])) for i in range(Nc)]
            pdfR = [max([fRi[j][i] for j in range(len(fRi))]) for i in range(len(R))]

        cdfR = np.cumsum(pdfR)

        return R, pdfR/cdfR[-1], cdfR/cdfR[-1]



    def rayleigh_iso_incs(self, inc, per, GMfp213, R_star):
        """Uses the Rayleigh distribution method from Fabrycky (2015) and the isotropic fraction from Mulders et al. (2018) to determine the system inclination and the distribution of mutual inclinations."""

        il = np.linspace(0, 180, 361)
        fi = np.zeros(len(il))
        rylgh = 2
        ibs = []
        fib = []
        incn = []
        incq = [i for i in inc if i != "?"]
        ib = 0

        if len(incq) == 0:
            Pinc = np.sin(il*math.pi/180)
            li = round(np.arccos(R_star*self.R_sun/(GMfp213*(per[0]*self.seconds_per_day)**(2/3)))*180/math.pi, 1)

            for i in range(len(il)):
                if il[i] >= li and il[i] <= 180 - li:
                    Pinc[i] = 0

            cdfi = np.cumsum(Pinc)

            return il, Pinc/cdfi[-1], cdfi/cdfi[-1], np.linspace(0, 180-li, int((180-li)*10)), li

        else:
            for case in [[False] + list(t) for t in list(itertools.product([False,True], repeat=len(incq)-1))]:
                incn.append([180-incq[i] if case[i] else incq[i] for i in range(0, len(incq))])

            fib = self.run_new_mp(self.inc_test, il, (inc, incn, sigmas))

            mv = 0
            ib = 0

            for i in range(len(fib)):
                cmax = max(fib[i])

                if cmax > mv:
                    ib = il[i]
                    mv = cmax
            
            #ibs, fib = self.ppr.create_processes("inc_test_old", (inc, il, incn, sigmas), -len(il), self.process_inc_data_old)
            #ib = ibs[np.where(fib == max(fib))[0][0]]

        if ib > 90:
            ilb = round(ib, 1)
            ib = round(180 - ib, 1)

        else:
            ilb = round(180 - ib, 1)
            ib = round(ib, 1)

        inew = np.linspace(0, ilb, int(ilb*10) + 1)
        rylghi = spst.rayleigh.pdf(inew, 0, rylgh)
        finew = [rylghi[np.where(np.isclose(inew, abs(ib - il[j])))[0][0]] for j in range(len(il) - 1)]
        finew.append(0.0)
        finew = np.array(finew)

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

        for j in range(len(finew)):
            fi[j] += finew[j]

        cdfi = np.cumsum(fi)

        return il, fi/cdfi[-1], cdfi/cdfi[-1], inew, ib



    def process_inc_data(self, data):
        """Processes data for the multithreading component"""
        
        return tuple([list(itertools.chain(*i)) for i in zip(*itertools.chain(data))])



    def inc_test_old(self, inc, il, incn, sigmas, j):
        """Tests the best system inclination."""
        ibs = []
        fib = []

        for k in range(len(incn)):
            test = 0

            for m in range(len(incn[k])):
                if self.config_parameters["inclination"] == "rayleigh_iso":
                    test += spst.rayleigh.pdf(abs(incn[k][m]-il[j]), 2)

                elif self.config_parameters["inclination"] == "syssim":
                    test += spst.lognorm.pdf(abs(incn[k][m]-il[j]), sigmas[len(inc)-1], scale=1.1*((len(inc)+1)/5)**-1.73)

            ibs.append(il[j])
            fib.append(test)
            
        return (ibs, fib)



    def inc_test(self, j, args):
        """Tests the best system inclination."""

        inc, incn, sigmas = args
        fib = []

        for k in range(len(incn)):
            test = 0

            for m in range(len(incn[k])):
                if self.config_parameters["inclination"] == "rayleigh_iso":
                    test += spst.rayleigh.pdf(abs(incn[k][m]-j), 2)

                elif self.config_parameters["inclination"] == "syssim":
                    test += spst.lognorm.pdf(abs(incn[k][m]-j), sigmas[len(inc)-1], scale=1.1*((len(inc)+1)/5)**-1.73)

            fib.append(test)
            
        return fib



    def syssim_incs(self, inc, per, GMfp213, R_star):
        """Uses the power-law distributions from He et al. (2020) to determine the median of the mutual inclination lognormal distribution."""

        il = np.linspace(0, 180, 361)
        fi = np.zeros(len(il))
        ibs = []
        fib = []
        incn = []
        incq = [i for i in inc if i != "?"]
        ib = 0
        sigmas = [0.84, 0.85, 0.86, 0.86, 0.85, 0.84, 0.81, 0.79, 0.77]

        if len(incq) == 0:
            Pinc = np.sin(il*math.pi/180)
            li = round(np.arccos(R_star*self.R_sun/(GMfp213*(per[0]*self.seconds_per_day)**(2/3)))*180/math.pi, 1)

            for i in range(len(il)):
                if il[i] >= li and il[i] <= 180 - li:
                    Pinc[i] = 0

            cdfi = np.cumsum(Pinc)

            return il, Pinc, cdfi/cdfi[-1], np.linspace(0, 180-li, int((180-li)*10)), li

        else:
            for case in [[False] + list(t) for t in list(itertools.product([False,True], repeat=len(incq)-1))]:
                incn.append([180-incq[i] if case[i] else incq[i] for i in range(0, len(incq))])
            
            fib = self.run_new_mp(self.inc_test, il, (inc, incn, sigmas))

            mv = 0
            ib = 0

            for i in range(len(fib)):
                cmax = max(fib[i])

                if cmax > mv:
                    ib = il[i]
                    mv = cmax
            
            #ibs, fib = self.ppr.create_processes("inc_test_old", (inc, il, incn, sigmas), -len(il), self.process_inc_data_old)
            #ib = ibs[np.where(fib == max(fib))[0][0]]

        if ib > 90:
            ilb = ib
            ib = 180 - ib

        else:
            ilb = 180 - ib

        inew = np.linspace(0, ilb, int(ilb*10) + 1)
        lgnrmi = spst.lognorm.pdf(inew, sigmas[len(inc)-1], scale=1.1*((len(inc)+1)/5)**-1.73)
        finew = [lgnrmi[np.where(np.isclose(inew, abs(ib - il[j]), atol=1e-1))[0][0]] for j in range(len(il) - 1)]
        finew.append(0.0)

        for j in range(len(finew)):
            fi[j] += finew[j]

        cdfi = np.cumsum(fi)

        return il, fi/cdfi[-1], cdfi/cdfi[-1], inew, ib



    def rayleigh_eccs(self):
        """Uses a Rayleigh distribution for the eccentricity."""

        el = np.linspace(0, 1, 1001)

        return el, spst.rayleigh.pdf(el, scale=0.02), spst.rayleigh.cdf(el, scale=0.02)



    def syssim_eccs(self, lecc):
        """Uses the power-law distributions from He et al. (2020) to determine the median of the eccentricity lognormal distribution."""

        el = np.linspace(0, 1, 1001)
        sigmas = [0.641, 0.670, 0.689, 0.701, 0.704, 0.685, 0.666, 0.632, 0.612, 0.587]

        return el, spst.lognorm.pdf(el, sigmas[lecc], scale=0.031*((lecc+1)/5)**-1.74), spst.lognorm.cdf(el, sigmas[lecc], scale=0.031*((lecc+1)/5)**-1.74)



    def otegi_mr(self, measurement, predict, force=None):
        """Uses a power-law MR to predict mass in Earth values from radius in Earth values or vice versa"""

        if predict == "radius":
            R_r = 1.03*measurement**0.29
            R_v = 0.7*measurement**0.63
            
            if force == "rocky":
                return R_r

            if measurement > 25 or (measurement > 5 and self.config_parameters["otegi_rho"] == "volatile") or force == "volatile":
                return R_v

            return R_r

        elif predict == "mass":
            M_r = 0.9*measurement**3.45
            M_v = 1.74*measurement**1.58

            if force == "rocky":
                return M_r

            if M_r > 25 or (M_r > 5 and self.config_parameters["otegi_rho"] == "volatile") or force == "volatile":
                return M_v

            return M_r



    def K3(self, P, M):
        """Calculates semi-major axis in cm using period in days and mass in solar masses"""
        
        return (self.G*M*self.M_sun*(P*self.seconds_per_day)**2/(4*math.pi**2))**(1/3)



    def create_sim_system(self, M_star, per, Pmc, mas, inc, ecc, Mmc, imc, emc, GMfp213):
        """Creates the simulated system for REBOUND analysis."""

        sim = rebound.Simulation()
        sim.units = ('day', 'au', 'Msun')
        sim.integrator = "mercurius"
        sim.add(m=M_star)

        pf, mf, incf, eccf, mind = self.create_sim_arrays(per, Pmc, mas, inc, ecc, Mmc, imc, emc)

        for i in range(len(pf)):
            sim.add(m=mf[i]*self.MES, a=GMfp213*(pf[i]*self.seconds_per_day)**(2/3)/self.au, e=eccf[i], inc=incf[i]*math.pi/180)

        return sim, pf, mf, incf, eccf, mind



    def create_sim_arrays(self, per, Pmc, mas, inc, ecc, Mmc, imc, emc):
        """Creates the simulated planetary parameter arrays for REBOUND or TTVFaster analysis."""

        if Pmc != None:
            pf = np.zeros(len(per) + 1)
            mf = np.zeros(len(per) + 1)
            incf = np.zeros(len(per) + 1)
            eccf = np.zeros(len(per) + 1)
            mind = 0
            mc_add = False
        
            for p in range(len(per)):
                if Pmc > per[p]:
                    mind = p + 1
                    pf[p] = per[p]
                    mf[p] = mas[p]
                    incf[p] = inc[p]*math.pi/180
                    eccf[p] = ecc[p]

                else:
                    if not mc_add:
                        pf[p] = Pmc
                        mf[p] = Mmc
                        incf[p] = imc*math.pi/180
                        eccf[p] = emc
                        mc_add = True

                    pf[p+1] = per[p]
                    mf[p+1] = mas[p]
                    incf[p+1] = inc[p]*math.pi/180
                    eccf[p+1] = ecc[p]

        else:
            pf = per
            mf = mas
            incf = inc
            eccf = ecc
            mind = 0

        return pf, mf, incf, eccf, mind



if __name__ == '__main__':
    dynamite()
