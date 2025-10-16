### DYNAmical Multi-planet Injection TEster (DYNAMITE) ###
### Plotting File ###
### Main Author: Jamie Dietrich ###
### Contact: jdietrich@asu.edu ###
### Contributing Authors: Ritvik Basant, Katelyn Ruppert ###
### 2025 July 9 ###
### Version 3.1 ###
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
import numpy as np
import scipy.stats as spst
import matplotlib.pyplot as plt
import astropy.constants as const
import matplotlib.ticker as mticker
import matplotlib.patches as mpatch

from datetime import datetime
from dynamite_targets_db import dynamite_targets_db

# deprecated: from mrexo import predict_from_measurement as pfm

class dynamite_plots:

    def __init__(self, Pk=[], P=[], PP=[], Rk=[], R=[], PR=[], ik=[], inc=[], Pi=[], ek=[], ecc=[], Pe=[], deltas=[], tdm=[], tdle=[], tdue=[], tpm=[], tple=[], tpue=[], targets=[], starvs=[], pers=[], rads=[], mass=[], eccs=[], system="", cfname="dynamite_config.txt"):
        """Sets up plotting routines"""

        self.seconds_per_day = 86400
        self.G = const.G.cgs.value
        self.sigma_sb = const.sigma_sb.cgs.value
        self.au = const.au.cgs.value
        self.M_sun = const.M_sun.cgs.value
        self.R_sun = const.R_sun.cgs.value
        self.L_sun = const.L_sun.cgs.value
        self.M_earth = const.M_earth.cgs.value
        self.R_earth = const.R_earth.cgs.value
        self.config_parameters = {}
        self.sdir = ""

        try:
            config_data = np.loadtxt(cfname, dtype=str, delimiter=':')

        except IOError:
            print("Error, configuration file not found!")
            exit()

        for i in range(len(config_data)):
            self.config_parameters[config_data[i, 0]] = config_data[i, 1] if config_data[i, 1].find("[") == -1 else ast.literal_eval(config_data[i, 1])

        if len(Pk) == 0:
            Pk, P, PP, Rk, R, PR, ik, inc, Pi, ek, ecc, Pe, deltas, tdm, tdle, tdue, tpm, tple, tpue, targets, _starvs, pers, rads, _mass, _eccs = self.read_saved_data()

        if self.config_parameters["mode"] != "single":
            self.sdir = mode + "/"

            if self.sdir != "" and not os.path.exists(self.sdir):
                os.makedirs(self.sdir)
                
            P = P[0]
            R = R[0]
            inc = inc[0]
            ecc = ecc[0]

        functions, args = [], []

        if self.config_parameters["plt_P_R"] == "True" and self.config_parameters["mode"] != "single" and self.config_parameters["mode"] != "all":
            functions.append("plot_P_R"), args.append((Pk, P, Rk, R, targets, pers, rads))

        if self.config_parameters["plt_tdtp"] == "True" and self.config_parameters["mode"] != "single":
            functions.append("plot_td_tp"), args.append((tdm, tdle, tdue, tpm, tple, tpue, targets))

        if self.config_parameters["plt_deltas"] == "True" and self.config_parameters["mode"] != "single":
            functions.append("plot_deltas"), args.append((deltas,))

        if self.config_parameters["plt_ratios"] == "True" and self.config_parameters["mode"] != "single":
            functions.append("plot_ratios"), args.append(([pers[i+1]/pers[i] for i in range(len(pers) - 1)],))

        if self.config_parameters["plt_indpars"] == "True" and self.config_parameters["mode"] == "single":
            self.sdir = self.config_parameters["system"].replace("\"", "") + "/"

            if self.sdir != "" and not os.path.exists(self.sdir):
                os.makedirs(self.sdir)
                
            self.plot_ind_params(Pk, P, PP, Rk, R, PR, ik, inc, Pi, ek, ecc, Pe, system.replace("\"", ""))

        print(datetime.now(), "Finishing Plots")



    def read_saved_data(self):
        """Reads in the saved data"""

        print(datetime.now(), "Reading in Saved Data")

        with np.load("saved_data.npz", allow_pickle=True) as data:
            Pk, P, PP, Rk, R, PR, ik, inc, Pi, ek, ecc, Pe, deltas, tdm, tdle, tdue, tpm, tple, tpue, targets, starvs, pers, rads, mass, eccs = [data["arr_" + str(i)] for i in range(len(data))]

        return Pk, P, PP, Rk, R, PR, ik, inc, Pi, ek, ecc, Pe, deltas, tdm, tdle, tdue, tpm, tple, tpue, targets, starvs, pers, rads, mass, eccs



    def plot_P_R(self, Pk, P, Rk, R, targets, pers, rads):
        """Plots the P, R, and PR distributions"""

        print(datetime.now(), "Creating Period-Planet Radius Figures")

        Pb = list(P)
        Pb.append(730.1)
        Rb = list(R)
        Rb.append(float(self.config_parameters["radmax"]) + 0.01)
        Pf = np.arange(0.1, 730.01, 0.01)
        Rf = np.arange(float(self.config_parameters["radmin"]), float(self.config_parameters["radmax"]) + 0.001, 0.001)

        if self.config_parameters["saved"] == "False":
            n_tot = []
            r_tot = []
            plot_fig_p = []
            plot_fig_r = []

            for i in range(len(Pk)):
                n, _, _ = plt.hist(Pk[i], bins=Pb, weights=np.ones(len(Pk[i]))*10 / len(Pk[i]), label="DYNAMITE Predictions")
                n_tot.append(n)
                plt.close()
                n, _, _ = plt.hist(Rk[i], bins=Rb, weights=np.ones(len(Rk[i]))*10 / len(Rk[i]), label="DYNAMITE Predictions")
                r_tot.append(n)
                plt.close()

            for i in range(len(n_tot)):
                for _ in range(200):
                    plot_fig_p.append(n_tot[i][:1855])
                    plot_fig_r.append(r_tot[i])

            plot_fig_p = np.array(plot_fig_p)
            plot_fig_r = np.array(plot_fig_r)
            np.savetxt("plot_fig_p_" + self.config_parameters["mode"] + ".txt", plot_fig_p)
            np.savetxt("plot_fig_r_" + self.config_parameters["mode"] + ".txt", plot_fig_r)

        elif self.config_parameters["saved"] == "True":
            plot_fig_p = np.loadtxt("plot_fig_p_" + self.config_parameters["mode"] + ".txt")
            plot_fig_r = np.loadtxt("plot_fig_r_" + self.config_parameters["mode"] + ".txt")

        for i in range(len(targets)):
            targets[i][1] = int(float(targets[i][1])*100 - 10)
            targets[i][2] = int(float(targets[i][2])*100)
            targets[i][3] = int(float(targets[i][3])*100)
            targets[i][4] = int(float(targets[i][4])*1000 - 100)
            targets[i][5] = int(float(targets[i][5])*1000)
            targets[i][6] = int(float(targets[i][6])*1000)

        try:
            pfpi = np.zeros((len(plot_fig_p), len(Pf)))

            for i in range(len(plot_fig_p)):
                pfpi[i] = np.interp(Pf, P[:1855], plot_fig_p[i])

            pfri = np.zeros((len(plot_fig_r), len(Rf)))

            for i in range(len(plot_fig_r)):
                pfri[i] = np.interp(Rf, R, plot_fig_r[i])

            prf = np.zeros((len(Rf), len(Pf)))

            for i in range(len(targets)):
                pp = [pfpi[i*200, j] for j in range(len(Pf))]
                pr = [pfri[i*200, j] for j in range(len(Rf))]
                prf += np.outer(pr, pp)

            fig, ax = plt.subplots(1,1, figsize=(20, 12))
            fig.suptitle(r"Period-Planet Radius Relative Likelihood", fontsize=30)
            img = ax.imshow(prf, cmap=plt.cm.viridis, origin="lower", aspect="auto")
            ax.set_xscale("log")
            xlabels = [0.1, 1, 10, 100]
            ylabels = [0, 1, 2, 3, 4, 5]
            ax.set_xticks([9, 90, 990, 9990])
            ax.set_yticks([-100, 900, 1900, 2900, 3900, 4900])
            ax.set_xticklabels(xlabels)
            ax.set_yticklabels(ylabels)
            ax.tick_params(labelsize=14)
            plt.scatter([int(targets[i][1]) for i in range(len(targets))], [int(targets[i][4]) for i in range(len(targets))], color="w", marker="*")

            for i in range(len(targets)):
                ple = int(targets[i][1]) - int(targets[i][2])
                pue = int(targets[i][1]) + int(targets[i][3])
                rle = int(targets[i][4]) - int(targets[i][5])
                rue = int(targets[i][4]) + int(targets[i][6])
                ellipse = mpatch.Ellipse((np.mean([ple, pue]), np.mean([rle, rue])), pue-ple, rue-rle, alpha=0.1, color="w")
                ax.add_patch(ellipse)

            bottom = []
            left = []
            right = []
            tl = []
            tr = []
            bl = []
            br = []
            
            if self.config_parameters["mode"] == "TESS":
                if self.config_parameters["period"] == "epos":
                    bottom = ["TOI 256", "TOI 286", "TOI 270", "TOI 402", "TOI 411", "TOI 712", "TOI 1269", "TOI 1468", "TOI 1726"]
                    left = ["TOI 119", "TOI 487", "TOI 1469"]
                    right = ["TOI 266", "TOI 396", "TOI 703", "TOI 1266", "TOI 1749"]
                    tl = ["TOI 1260", "TOI 1438"]
                    tr = ["TOI 763"]
                    bl = ["TOI 712"]
                    br = ["TOI 713"]

                elif self.config_parameters["period"] == "syssim":
                    bottom = ["TOI 119", "TOI 286", "TOI 487", "TOI 736", "TOI 797", "TOI 1453", "TOI 1720", "TOI 1726", "TOI 1749"]
                    left = ["TOI 256", "TOI 1269"]
                    right = ["TOI 431", "TOI 703", "TOI 1277", "TOI 1730"]
                    tl = ["TOI 125"]
                    tr = ["TOI 266"]
                    bl = ["TOI 836", "TOI 1260"]
                    br = ["TOI 1438"]

            else:
                sep = " " if self.config_parameters["mode"] == "all" else "-"

                for i in range(len(targets)):
                    if targets[i][0] in bottom:
                        plt.annotate(targets[i][0][targets[i][0].find(sep)+1:], (int(targets[i][1]), int(targets[i][4])), color="w", textcoords="offset points", xytext=(0,-20), ha="center", weight='bold', fontsize=16)

                    elif targets[i][0] in right:
                        plt.annotate(targets[i][0][targets[i][0].find(sep)+1:], (int(targets[i][1]), int(targets[i][4])), color="w", textcoords="offset points", xytext=(25,-7), ha="center", weight='bold', fontsize=16)

                    elif targets[i][0] in left:
                        plt.annotate(targets[i][0][targets[i][0].find(sep)+1:], (int(targets[i][1]), int(targets[i][4])), color="w", textcoords="offset points", xytext=(-30,-7), ha="center", weight='bold', fontsize=16)

                    elif targets[i][0] in tl:
                        plt.annotate(targets[i][0][targets[i][0].find(sep)+1:], (int(targets[i][1]), int(targets[i][4])), color="w", textcoords="offset points", xytext=(-25,7), ha="center", weight='bold', fontsize=16)

                    elif targets[i][0] in tr:
                        plt.annotate(targets[i][0][targets[i][0].find(sep)+1:], (int(targets[i][1]), int(targets[i][4])), color="w", textcoords="offset points", xytext=(12,5), ha="center", weight='bold', fontsize=16)

                    elif targets[i][0] in bl:
                        plt.annotate(targets[i][0][targets[i][0].find(sep)+1:], (int(targets[i][1]), int(targets[i][4])), color="w", textcoords="offset points", xytext=(-15,-20), ha="center", weight='bold', fontsize=16)

                    elif targets[i][0] in br:
                        plt.annotate(targets[i][0][targets[i][0].find(sep)+1:], (int(targets[i][1]), int(targets[i][4])), color="w", textcoords="offset points", xytext=(15,-20), ha="center", weight='bold', fontsize=16)

                    else:
                        plt.annotate(targets[i][0][targets[i][0].find(sep)+1:], (int(targets[i][1]), int(targets[i][4])), color="w", textcoords="offset points", xytext=(0,8), ha="center", weight='bold', fontsize=16)

            cb = fig.colorbar(img)
            cb.set_label(label="Probability normalized to 1 injected planet per system", size=20)
            cb.ax.tick_params(labelsize=14)
            plt.xlabel("Log Period (days)", fontsize=20)
            plt.ylabel(r"Radius ($R_{\oplus}$)", fontsize=20)
            plt.ylim(750, 4200)

            if self.config_parameters["period"] == "epos":
                plt.xlim(150, 11000)

            elif self.config_parameters["period"] == "syssim":
                plt.xlim(40, 11000)

            plt.savefig(self.sdir + self.sdir + "PRfig_" + self.config_parameters["mode"] + "_" + self.config_parameters["period"] + "_ell.png", bbox_inches='tight')

            if self.config_parameters["show_plots"] == "True":
                plt.show()

            elif self.config_parameters["show_plots"] == "False":
                plt.close()

            fig, ax = plt.subplots(1,1, figsize=(10,12))
            fig.suptitle("Period Relative Likelihoods for " + self.config_parameters["mode"] + " Systems", fontsize=30)
            img = ax.imshow(pfpi, cmap=plt.cm.Blues, origin="lower", aspect="auto")

            if self.config_parameters["plt_P_scale"] == "linear":
                xlabels = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
                ax.set_xticks([-50, 1500, 3500, 5500, 7500, 9500, 11500, 13500, 15500, 17500])
                ax.set_xticklabels(xlabels)
                ax.set_xlim(-500, 18000)

            elif self.config_parameters["plt_P_scale"] == "log":
                ax.set_xscale("log")
                xlabels = [0.1, 1, 10, 100]
                ax.set_xticks([0, 90, 990, 9990])
                ax.set_xticklabels(xlabels)
                ax.set_xlim(10, 73000)

            ylabels = [targets[i][0] for i in range(len(targets))]
            ax.set_yticks(np.arange(100, len(ylabels)*200 + 100, 200))
            ax.set_yticklabels(ylabels)
            ax.tick_params(labelsize=12)
            cb = fig.colorbar(img)
            cb.set_label(label="Probability normalized to 1 injected planet per system", size=16)
            cb.ax.tick_params(labelsize=12)
            plt.xlabel("Period (days)", fontsize=16)
            plt.ylabel("Systems", fontsize=16)
            pc = 100
            plt.scatter(np.where(np.isclose(Pf, pers[0][0][0], atol=(0.005 if pers[0][0][0] < 1 else 0.05)))[0][0], pc, color="r", s=int(round(rads[0][0][0]*20, 0)), label="Known Planets")

            for i in range(len(pers)):
                for j in range(len(pers[i])):
                    plt.scatter(np.where(np.isclose(Pf, pers[i][j][0], atol=0.006))[0][0], pc, color="r", s=int(round(rads[i][j][0]*20, 0)))

                pc += 200

            plt.savefig(self.sdir + "logPfig_" + self.config_parameters["mode"] + "_" + self.config_parameters["period"] + ".png", bbox_inches='tight')

            if self.config_parameters["show_plots"] == "True":
                plt.show()

            elif self.config_parameters["show_plots"] == "False":
                plt.close()
           
            fig, ax = plt.subplots(1,1, figsize=(10,12))
            fig.suptitle("Planet Radius Relative Likelihoods for " + self.config_parameters["mode"] + " Systems", fontsize=30)
            img = ax.imshow(pfri, cmap=plt.cm.Blues, origin="lower", aspect="auto")
            xlabels = np.arange(int(self.config_parameters["radmax"]) + 1)
            ylabels = [targets[i][0] for i in range(len(targets))]
            ax.set_xticks(np.arange(-100, (int(self.config_parameters["radmax"]) + 1)*1000 - 100, 1000))
            ax.set_yticks(np.arange(100, len(ylabels)*200 + 100, 200))
            ax.set_xticklabels(xlabels)
            ax.set_yticklabels(ylabels)
            ax.tick_params(labelsize=12)
            cb = fig.colorbar(img)
            cb.set_label(label="Probability normalized to 1 injected planet per system", size=16)
            cb.ax.tick_params(labelsize=12)
            plt.xlabel(r"Radius ($R_{\oplus}$)", fontsize=16)
            plt.ylabel("Systems", fontsize=16)
            rc = 100
            plt.scatter(np.where(np.isclose(Rf, rads[0][0][0], atol=0.006))[0][0], rc, color="r", label="Known Planets")

            for i in range(len(rads)):
                for j in range(len(rads[i])):
                    plt.scatter(np.where(np.isclose(Rf, rads[i][j][0], atol=0.006))[0][0], rc, color="r")

                rc += 200

            plt.savefig(self.sdir + "Rfig_" + self.config_parameters["mode"] + ".png", bbox_inches='tight')

            if self.config_parameters["show_plots"] == "True":
                plt.show()

            elif self.config_parameters["show_plots"] == "False":
                plt.close()

        except Exception as e:
            print(e)
            print("P-R plots could not be created! Check DYNAMITE's mode and/or most recent saved data")



    def plot_deltas(self, deltas):
        """Plots the probability as a function of the separation stability criterion in terms of mutual Hill radii"""

        print(datetime.now(), "Creating Stability Criterion Figure")

        if self.config_parameters["saved"] == "False":
            deltas = [sum(x) for x in zip(*deltas)]
            deltas = [deltas[i]/46 for i in range(len(deltas))]
            Du = np.arange(0, len(deltas))
            deltas = np.insert(deltas, 8, np.zeros(99))
            Du = np.insert(Du, 8, np.linspace(7.01, 7.99, 99))
            x = [[Du[i], deltas[i]] for i in range(len(deltas))]
            np.savetxt("Deltas_" + self.config_parameters["mode"] + "_" + self.config_parameters["period"] + ".txt", x, delimiter="\t")

        elif self.config_parameters["saved"] == "True":
            Du, deltas = np.loadtxt("Deltas_" + self.config_parameters["mode"] + "_" + self.config_parameters["period"] + ".txt", delimiter="\t", unpack=True)

        plt.plot(Du, deltas)
        plt.xlabel(r"$\Delta$ (Mutual Hill Radii)")
        plt.ylabel("Frequency")
        plt.title("Planet Pair Separation", fontsize=30)
        plt.savefig(self.sdir + "Deltas_" + self.config_parameters["mode"] + "_" + self.config_parameters["period"] + ".png")

        if self.config_parameters["show_plots"] == "True":
            plt.show()

        elif self.config_parameters["show_plots"] == "False":
            plt.close()



    def plot_ratios(self, ratios):
        """Plots the period ratios for the current subset of data"""

        print(datetime.now(), "Creating Period Ratio Figure")

        ratios = np.array(ratios)
        fig, ax = plt.subplots(figsize=(12,8))
        plt.hist(ratios, bins=np.logspace(0,1.4,15), label="TESS Systems")

        logD = -0.9
        sigma = 0.4
        PRgrid = np.logspace(0,1.4,100)

        with np.errstate(divide='ignore'):
            Dgrid = np.log(2.*(PRgrid**(2./3.)-1.)/(PRgrid**(2./3.)+1.))

        Dgrid[0] = -4
        pdfP = spst.norm(logD, sigma).pdf(Dgrid)
        plt.plot(PRgrid, 17*pdfP, linewidth=4, label=r"$Kepler$ PDF fit")
        plt.xscale("log")
        plt.xlabel(r"Period Ratio $P_i/P_{i-1}$", fontsize=20)
        plt.ylabel("Occurrence", fontsize=20)
        plt.legend(fontsize=16)
        ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
        ax.set_xticklabels([0.1,1,10])
        ax.tick_params(labelsize=14)
        ax.set_yticks([0, 5, 10, 15, 20])
        fig.suptitle(r"Subsample vs $Kepler$ Period Ratios", fontsize=30)
        plt.xlim(1, 10**1.4)
        plt.savefig(self.sdir + self.config_parameters["mode"] + "_Kepler_period_ratios.png")

        if self.config_parameters["show_plots"] == "True":
            plt.show()

        elif self.config_parameters["show_plots"] == "False":
            plt.close()



    def plot_td_tp(self, tdm, tdle, tdue, tpm, tple, tpue, targets):
        """Plots the transit probability vs the transit depth for each system in the current data subset"""

        print(datetime.now(), "Creating Transit Depth-Probability Figure")

        fig, ax = plt.subplots(figsize=(12,8))
        plt.plot(tdm, tpm, "o")

        names = [targets[i][0] for i in range(len(targets))]
        bunch = []
        bottom = []
        left = []
        right = []
        sep = " " if (self.config_parameters["mode"] == "all" or self.config_parameters["mode"] == "tess") else "-"
        for i in range(len(names)):
            tdle1 = tdm[i] - tdle[i]
            tdue1 = tdm[i] + tdue[i]
            tple1 = tpm[i] - tple[i]
            tpue1 = tpm[i] + tpue[i]
            ellipse = mpatch.Ellipse((np.mean([tdle1, tdue1]), np.mean([tple1, tpue1])), tdue1-tdle1, tpue1-tple1, alpha=0.1)
            ax.add_patch(ellipse)

            if names[i] not in bunch:
                if names[i] in bottom:
                    plt.annotate(names[i][names[i].find(sep)+1:], (tdm[i], tpm[i]), textcoords="offset points", xytext=(0,-20), ha='center', fontsize=16)

                elif names[i] in left:
                    plt.annotate(names[i][names[i].find(sep)+1:], (tdm[i], tpm[i]), textcoords="offset points", xytext=(-30,-5), ha='center', fontsize=16)

                elif names[i] in right:
                    plt.annotate(names[i][names[i].find(sep)+1:], (tdm[i], tpm[i]), textcoords="offset points", xytext=(25,-5), ha='center', fontsize=16)

                else:
                    plt.annotate(names[i][names[i].find(sep)+1:], (tdm[i], tpm[i]), textcoords="offset points", xytext=(0,10), ha='center', fontsize=16)

            else:
                if names[i] in bottom:
                    plt.annotate(names[i][names[i].find(sep)+1:], (tdm[i], tpm[i]), textcoords="offset points", xytext=(0,-20), ha='center', fontsize=16)

                elif names[i] in left:
                    plt.annotate(names[i][names[i].find(sep)+1:], (tdm[i], tpm[i]), textcoords="offset points", xytext=(-30,-5), ha='center', fontsize=16)

                elif names[i] in right:
                    plt.annotate(names[i][names[i].find(sep)+1:], (tdm[i], tpm[i]), textcoords="offset points", xytext=(30,-5), ha='center', fontsize=16)

                else:
                    plt.annotate(names[i][names[i].find(sep)+1:], (tdm[i], tpm[i]), textcoords="offset points", xytext=(0,10), ha='center', fontsize=16)

        ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
        plt.xlim(150, max(tdm) + max(tdue))
        ax.tick_params(labelsize=14)
        plt.xlabel("Transit Depth (ppm)", fontsize=20)
        plt.ylabel("Transit Probability", fontsize=20)
        fig.suptitle("Transit Probability vs Transit Depth", fontsize=30)
        plt.savefig(self.sdir + self.config_parameters["mode"] + "_td_tp.png")

        if self.config_parameters["show_plots"] == "True":
            plt.show()

        elif self.config_parameters["show_plots"] == "False":
            plt.close()



    def plot_ind_params(self, Pk, P, PP, Rk, R, PR, ik, inc, Pi, ek, ecc, Pe, system):
        """Plots histograms of each individual distribution"""

        print(datetime.now(), "Creating Histograms for", system)

        color_scheme = ["#" + self.config_parameters["plt_colors"][i] for i in range(len(self.config_parameters["plt_colors"]))]
        targets_dict = dynamite_targets_db(self.config_parameters["targets_db"]).get_targets(self.config_parameters["mode"], system, self.config_parameters["radmax"], self.config_parameters["removed"])
        alt_name, Rs, Ms, Ts, target, names, _ = self.set_up(targets_dict, system)
        target = target[:-1] if len(self.config_parameters["removed"]) > 0 else target

        for i in target:
            if i[3] == "Msini" and i[4] != "?":
                i[2] = i[2]/np.sin(i[4]*np.pi/180)
                i[1] = self.mr_appends("radius", i[2])

        p1, p11, p12, p2, p3, r1, r11, r12, r2, r3, m1, m11, m12, m2, m3, i1, i11, i12, i2, i3, e1, e11, e12, e2, e3, l1, l11, l12, l2, l3 = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
        
        for i in range(len(target)):
            if target[i][0] not in self.config_parameters["removed"]:
                if target[i][6] == 1 or (system != "Inner Solar System" and target[i][4] != "?" and abs(math.cos(target[i][4]*math.pi/180)) * self.K3(target[i][0], Ms) / (Rs*self.R_sun) <= 1) or system.find("TOI") != -1 or (system.find("Kepler") != -1 and target[i][4] != "?" and abs(math.cos(target[i][4]*math.pi/180)) *
                            self.K3(target[i][0], Ms) / (Rs*self.R_sun) <= 1) or (system == "Inner Solar System" or system == "eps Eri" and target[i][1] < 1.5):
                    p1.append(target[i][0])
                    r1.append(target[i][1] if target[i][1] != "?" else self.mr_appends("radius", float(target[i][2])))
                    m1.append(target[i][2] if target[i][2] != "?" else self.mr_appends("mass", float(target[i][1])))
                    i1.append(target[i][4])
                    e1.append(target[i][5])
                    l1.append(names[i])

                else:
                    p12.append(target[i][0])
                    r12.append(target[i][1] if target[i][1] != "?" else self.mr_appends("radius", float(target[i][2])))
                    m12.append(target[i][2] if target[i][2] != "?" else self.mr_appends("mass", float(target[i][1])))
                    i12.append(target[i][4])
                    e12.append(target[i][5])
                    l12.append(names[i])

        if len(self.config_parameters["additional"][0]) > 0:
            for i in self.config_parameters["additional"]:
                p11.append(i[0][0])
                r11.append(i[1][0] if i[1][0] != "?" else self.mr_appends("radius", i[2][0]))
                m11.append(i[2][0] if i[2][0] != "?" else self.mr_appends("mass", i[1][0]))
                i11.append(i[3][0])
                e11.append(i[4][0])
                l11.append(i[5][0])

        if len(self.config_parameters["unconfirmed"][0]) > 0:
            for i in self.config_parameters["unconfirmed"]:
                if i[2] == "?":
                    mv = self.mr_appends("mass", i[1])

                else:
                    if i[2][3] == "Msini" and i[3] != "?":
                        mv = i[2]/math.sin(i[3]*math.pi/180)

                    else:
                        mv = i[2]

                p3.append(i[0])
                r3.append(i[1] if i[1] != "?" else self.mr_appends("radius", mv))
                m3.append(mv)
                i3.append(i[3])
                e3.append(i[4])
                l3.append(i[5])

        if len(self.config_parameters["removed"]) > 0:
            for i in range(len(self.config_parameters["removed"])):
                p2.append([target[j][0] for j in range(len(target)) if target[j][0] == self.config_parameters["removed"][i]][0])
                i2.append([target[j][4] for j in range(len(target)) if target[j][0] == self.config_parameters["removed"][i]][0])
                e2.append([target[j][5] for j in range(len(target)) if target[j][0] == self.config_parameters["removed"][i]][0])
                l2.append([names[j] for j in range(len(target)) if target[j][0] == self.config_parameters["removed"][i]][0])

                for j in target:
                    if j[0] == self.config_parameters["removed"][i]:
                        r2.append(j[1] if j[1] != "?" else self.mr_appends("radius", j[2]))
                        m2.append(j[2] if j[2] != "?" else self.mr_appends("mass", j[1]))

        per = np.hstack([p1, p11, p12, p2, p3])
        rad = np.hstack([r1, r11, r12, r2, r3])
        mas = np.hstack([m1, m11, m12, m2, m3])

        fig, ax = plt.subplots(figsize=(12, 8))
        hands = []
        labs = []
        ul = 0

        ax1 = ax.twinx()
        m, t, ieb, oeb = dynamite_targets_db(self.config_parameters["targets_db"]).get_isochrones(Ms)

        if np.array_equal(ieb[0], ieb[1]):
            ie = ieb[0]

        else:
            ie = [np.interp(Ms, m, [ieb[0][e], ieb[1][e]]) for e in range(len(ieb[0]))]

        if np.array_equal(oeb[0], oeb[1]):
            oe = oeb[0]

        else:
            oe = [np.interp(Ms, m, [oeb[0][e], oeb[1][e]]) for e in range(len(oeb[0]))]

        pi = [(x**1.5/Ms**0.5)*365 for x in ie]
        po = [(x**1.5/Ms**0.5)*365 for x in oe]
        ax1.plot(pi, t, color=color_scheme[6], alpha=0.5)
        ax1.plot(po, t, color=color_scheme[6], alpha=0.5)
        ax1.fill_betweenx(t, pi, po, color=color_scheme[6], alpha=0.25)
        ax1.set_yscale("log")
        ax1.invert_yaxis()
        ax.set_zorder(ax1.get_zorder()-1)
        ax.patch.set_visible(False)
        xlim = (0.5, 1e5 if max(per) > 5e4 else 2*max(per) if max(per) > max(po) else 2*max(po))

        if self.config_parameters["ind_P"] == "linear_zoom":
            n, _, h1 = ax.hist(Pk, bins=np.arange(0.5, round(5*max(per)) + 1, 0.5), weights=np.ones(len(Pk))*2 / (len(Pk)), color=color_scheme[0] if self.config_parameters["stability"] == "hill" else color_scheme[2])
            hands.append(h1[0])
            labs = [(r"$\tau$ Ceti" if system == "tau Cet" else system) +  " Predictions - " + (r"$\bf{Period}$ $\bf{Ratios}$" if self.config_parameters["period"] == "epos" else r"$\bf{Clustered}$ $\bf{Periods}$")]
            if self.config_parameters["plt_PDFs"] == "True":
                h2 = ax.plot(P, PP, color=color_scheme[1], linewidth=3, zorder=ax.get_zorder()+1)
                hands.append(h2[0])
                labs.append("Probability Density Function")

            ul = np.amax(n)*1.5
            plt_savename = system.replace(" ", "") + "_P_linear_zoom_" + self.config_parameters["period"] + "_" + self.config_parameters["mass_radius"] + ("_" + self.config_parameters["otegi_rho"] if self.config_parameters["mass_radius"] == "otegi" else "") + "_" + self.config_parameters["inclination"] + "_" + self.config_parameters["eccentricity"] + "_" + self.config_parameters["stability"] + " .png"

        elif self.config_parameters["ind_P"] == "linear":
            Pb = list(P)
            n, _, h1 = ax.hist(Pk, bins=Pb, weights=np.ones(len(Pk))*10 / len(Pk), color=color_scheme[0] if self.config_parameters["stability"] == "hill" else color_scheme[2])
            hands.append(h1[0])
            labs = [(r"$\tau$ Ceti" if system == "tau Cet" else system) +  " Predictions - " + (r"$\bf{Period}$ $\bf{Ratios}$" if self.config_parameters["period"] == "epos" else r"$\bf{Clustered}$ $\bf{Periods}$")]

            if self.config_parameters["plt_PDFs"] == "True":
                h2 = ax.plot(P, PP, color=color_scheme[1], linewidth=3, zorder=ax.get_zorder()+1)
                hands.append(h2[0])
                labs.append("Probability Density Function")

            ul = np.amax(n)*1.5
            plt_savename = system.replace(" ", "") + "_P_linear_" + self.config_parameters["period"] + "_" + self.config_parameters["mass_radius"] + ("_" + self.config_parameters["otegi_rho"] if self.config_parameters["mass_radius"] == "otegi" else "") + "_" + self.config_parameters["inclination"] + "_" + self.config_parameters["eccentricity"] + "_" + self.config_parameters["stability"] + ".png"

        elif self.config_parameters["ind_P"] == "log":
            if system == "Inner Solar System":
                if self.config_parameters["stability"] == "hill":
                    bins = np.hstack([np.logspace(-0.3, 0.5, 18), np.logspace(0.51, 3.66, 211)])

                else:
                    bins = np.hstack([np.logspace(-0.3, 0.5, 18), np.logspace(0.51, 3.66, 316)])

            elif system == "Proxima Centauri":
                if self.config_parameters["stability"] == "hill":
                    bins = np.hstack([np.logspace(-0.3, 0.5, 18), np.logspace(0.51, 3.15, 177)])

                else:
                    bins = np.hstack([np.logspace(-0.3, 0.5, 18), np.logspace(0.51, 3.15, 265)])

            elif system == "eps Eri":
                if self.config_parameters["stability"] == "hill":
                    bins = np.logspace(2, 3.66, 112)

                else:
                    bins = np.logspace(2, 3.66, 167)

            elif system == "HD 219134":
                bins = np.hstack([np.logspace(-0.3, 0.5, 18), np.logspace(0.51, 2.9, 181)])
                
            else:
                if self.config_parameters["stability"] == "hill":
                    bins = np.hstack([np.logspace(-0.3, 0.5, 18), np.logspace(0.51, 2.9, 121)])

                else:
                    bins = np.hstack([np.logspace(-0.3, 0.5, 18), np.logspace(0.51, 2.9, 241)])

            widths = [bins[i+1] - bins[i] for i in range(len(bins)-1)]
            hist, _ = np.histogram(Pk, bins=bins)
            hist_norm = hist/widths
            h1 = ax.bar(bins[:-1], hist_norm/len(Pk), widths, color=color_scheme[0] if self.config_parameters["stability"] == "hill" else color_scheme[2])
            hands.append(h1)
            labs = [(r"$\tau$ Ceti" if system == "tau Cet" else (system[:system.index("test") - 1] if system.find("test") != -1 else system)) +  " Predictions - " + (r"$\bf{Period}$ $\bf{Ratios}$" if self.config_parameters["period"] == "epos" else r"$\bf{Clustered}$ $\bf{Periods}$")]
            PPl = np.interp(bins, P, PP)

            if self.config_parameters["plt_PDFs"] == "True":
                h2 = ax.plot(bins, PPl, color=color_scheme[1], linewidth=3, zorder=ax.get_zorder()+1)
                hands.append(h2[0])
                labs.append("Probability Density Function")

            ul = np.amax(hist_norm)*1.5/len(Pk)
            ax.set_xscale("log")
            plt_savename = system.replace(" ", "") + "_P_log_" + self.config_parameters["period"] + "_" + self.config_parameters["mass_radius"] + ("_" + self.config_parameters["otegi_rho"] if self.config_parameters["mass_radius"] == "otegi" else "") + "_" + self.config_parameters["inclination"] + "_" + self.config_parameters["eccentricity"] + "_" + self.config_parameters["stability"] + ".png"

        annots = []

        if len(p1) > 0:
            h3 = ax.scatter(p1, np.ones(len(p1))*ul/10, c=color_scheme[4], edgecolors="k", s=[r1[i]*200 for i in range(len(r1))], linewidth=2)
            hands.append(h3)

            if system == "Inner Solar System":
                labs.append("Terrestrial planets" if len(p1) > 1 else "Terrestrial planet")

            else:
                labs.append("Known transiting planets" if len(p1) > 1 else "Known transiting planet")

            annots = self.plot_annots(annots, p1, ul, l1, xlim, ("log" if self.config_parameters["ind_P"] == "log" else "linear"), ax)

        if len(p11) > 0:
            h4 = ax.scatter(p11, np.ones(len(p11))*ul/10, c=color_scheme[3], edgecolors="k", marker="s", s=[r11[i]*200 for i in range(len(r11))], linewidth=2)
            hands.append(h4)
            labs.append("Additional inserted planets" if len(p11) > 1 else "Additional inserted planet")
            annots = self.plot_annots(annots, p11, ul, l11, xlim, ("log" if self.config_parameters["ind_P"] == "log" else "linear"), ax)

        if len(p12) > 0:
            h5 = ax.scatter(p12, np.ones(len(p12))*ul/10, c=color_scheme[4], edgecolors="k", marker="^", s=[r12[i]*200 for i in range(len(r12))], linewidth=2)
            hands.append(h5)

            if system == "Inner Solar System":
                labs.append("Gas giant planets" if len(p12) > 1 else "Gas giant planet")

            else:
                labs.append("Known non-transiting planets" if len(p12) > 1 else "Known non-transiting planet")

            annots = self.plot_annots(annots, p12, ul, l12, xlim, ("log" if self.config_parameters["ind_P"] == "log" else "linear"), ax)

        if len(p2) > 0:
            h6 = ax.scatter(p2, np.ones(len(p2))*ul/10, c="white", edgecolors="k", marker="X", s=[r2[i]*200 for i in range(len(r2))], linewidth=2)
            hands.append(h6)
            labs.append("Known planets removed" if len(p2) > 1 else "Known planet removed")
            annots = self.plot_annots(annots, p2, ul, l2, xlim, ("log" if self.config_parameters["ind_P"] == "log" else "linear"), ax)

        if len(p3) > 0:
            h7 = ax.scatter(p3, np.ones(len(p3))*ul/10, c=color_scheme[5], edgecolors="k", marker='d', s=[r3[i]*200 for i in range(len(r3))], linewidth=2)
            hands.append(h7)
            labs.append("Unconfirmed planet candidates" if len(p3) > 1 else "Unconfirmed planet candidate")
            annots = self.plot_annots(annots, p3, ul, l3, xlim, ("log" if self.config_parameters["ind_P"] == "log" else "linear"), ax)

        hztl = [np.interp(0.003, t, pi), np.interp(0.003, t, po)]
        ax.annotate("Habitable\nZone", (hztl[0]*math.sqrt(hztl[1]/hztl[0]), ul*0.75), textcoords="offset points", xytext=(0, 0), ha="center", weight='bold', fontsize=11)
        ax.set_xlabel("Period (days)", fontsize=24)
        ax.set_ylabel("Relative Likelihood", fontsize=24)
        ax1.set_ylabel("Age (Gyr)", fontsize=24, rotation='vertical')
        lgnd = ax.legend(hands, labs, loc=9, fontsize=14, ncol=2, markerscale=0.75)

        for i in lgnd.legend_handles:
            i._sizes = [200]

        ax.tick_params(labelsize=18)
        ax1.tick_params(axis='y', labelsize=18)
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylim(0, ul)
        tx = ax.yaxis.get_offset_text()
        tx.set_fontsize(20)
        tx.set_position((-0.1,1.05))
        ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
        fig.suptitle((r"$\tau$ Ceti" if system == "tau Cet" else (system[:system.index("test") - 1] if system.find("test") != -1 else system)) + " Period Relative Likelihood", fontsize=30)
        plt.savefig(self.sdir + plt_savename)

        if self.config_parameters["show_plots"] == "True":
            plt.show()

        elif self.config_parameters["show_plots"] == "False":
            plt.close()

        fig, ax = plt.subplots(figsize=(12,8))
        Rk = np.array(Rk)
        Rb = list(R)

        if self.config_parameters["use_mass"] == "True":
            if self.config_parameters["mass_radius"] == "otegi":
                Mk = [self.otegi_mr(Rk[i], 'mass') for i in range(len(Rk))]
                Mb1 = sorted([self.otegi_mr(Rb[i], 'mass') for i in range(len(Rb))])
                Mb = np.arange(Mb1[0], Mb1[-1], 0.25)

            n, _, h1 = plt.hist(Mk, bins=Mb, weights=np.ones(len(Mk))*10 / len(Mk), color=color_scheme[0] if self.config_parameters["stability"] == "hill" else color_scheme[2])
            plt.xlabel(r"Mass ($M_{\oplus}$)", fontsize=30)

        else:
            n, _, h1 = plt.hist(Rk, bins=Rb, weights=np.ones(len(Rk))*10 / len(Rk), color=color_scheme[0] if self.config_parameters["stability"] == "hill" else color_scheme[2])
            plt.xlabel(r"Radius ($R_{\oplus}$)", fontsize=30)

        ul = 1.5*np.amax(n)

        if self.config_parameters["ind_R"] == "linear_zoom":
            plt_savename = system.replace(" ", "") + "_" + ("R" if self.config_parameters["use_mass"] == "False" else "M") + "_linear_zoom" + ("_Rearth" if self.config_parameters["show_Rearth"] == "True" else "") + "_" + self.config_parameters["mass_radius"] + ("_" + self.config_parameters["otegi_rho"] if self.config_parameters["mass_radius"] == "otegi" else "") + "_" + self.config_parameters["inclination"] + "_" + self.config_parameters["eccentricity"] + "_" + self.config_parameters["stability"] + ".png"

            if self.config_parameters["show_Rearth"] == "True":
                if self.config_parameters["use_mass"] == "False":
                    xlim = [0.75 if min(rad) > 1 else min(rad) - 0.25, math.ceil(max(rad))]

                else:
                    xlim = [0.75 if min(mas) > 1 else min(mas) - 0.25, math.ceil(max(mas))]

            else:
                if self.config_parameters["use_mass"] == "False":
                    xlim = [math.floor(min(rad)), math.ceil(max(rad))]

                else:
                    xlim = [math.floor(min(mas)), math.ceil(max(mas))]

        elif self.config_parameters["ind_R"] == "linear":
            plt_savename = system.replace(" ", "") + "_" + ("R" if self.config_parameters["use_mass"] == "False" else "M") + "_linear_" + self.config_parameters["mass_radius"] + ("_" + self.config_parameters["otegi_rho"] if self.config_parameters["mass_radius"] == "otegi" else "") + "_" + self.config_parameters["inclination"] + "_" + self.config_parameters["eccentricity"] + "_" + self.config_parameters["stability"] + ".png"
            ullim = np.where(n > 0.01*ul)[0]
            lims = [Rb[ullim[0]], Rb[ullim[-1]]] if self.config_parameters["use_mass"] == "False" else [Mb[ullim[0]], Mb[ullim[-1]]]

            if self.config_parameters["use_mass"] == "False":
                xlim = (min(lims[0], min(rad)) - 0.25, max(lims[1], max(rad)) + 0.25)

            else:
                xlim = (min(lims[0], min(mas)) - 0.25, max(lims[1], max(mas)) + 0.25)

        hands = [h1[0]]
        labs = [(r"$\tau$ Ceti" if system == "tau Cet" else system) + " Predictions - " + (r"$\bf{NP}$" if self.config_parameters["mass_radius"] == "mrexo" else "\"" + r"$\bf{Rocky}$" + "\"" if self.config_parameters["mass_radius"] == "otegi" and self.config_parameters["otegi_rho"] == "rocky" else "\"" + r"$\bf{Volatile}$" + "\"")]
        annots = []

        if self.config_parameters["plt_PDFs"] == "True":
            if self.config_parameters["use_mass"] == "False":
                h2 = plt.plot(R, PR, color=color_scheme[1], linewidth=3)
                hands.append(h2[0])
                labs.append("Probability Density Function")

        if len(r1) > 0:
            if self.config_parameters["use_mass"] == "False":
                h3 = plt.scatter(r1, np.ones(len(r1))*ul/10, c=color_scheme[4], edgecolors="k", s=200, linewidth=2, zorder=2)
                annots = self.plot_annots(annots, r1, ul, l1, xlim)

            else:
                h3 = plt.scatter(m1, np.ones(len(m1))*ul/10, c=color_scheme[4], edgecolors="k", s=200, linewidth=2, zorder=2)
                annots = self.plot_annots(annots, m1, ul, l1, xlim)

            hands.append(h3)

            if system == "Inner Solar System":
                labs.append("Terrestrial planets" if len(r1) > 1 else "Terrestrial planet")

            else:
                labs.append("Known transiting planets" if len(r1) > 1 else "Known transiting planet")


        if len(r11) > 0:
            if self.config_parameters["use_mass"] == "False":
                h4 = plt.scatter(r11, np.ones(len(r11))*ul/10, c=color_scheme[3], marker="s", edgecolors="k", s=200, linewidth=2, zorder=2)
                annots = self.plot_annots(annots, r11, ul, l11, xlim)

            else:
                h4 = plt.scatter(m11, np.ones(len(m11))*ul/10, c=color_scheme[3], marker="s", edgecolors="k", s=200, linewidth=2, zorder=2)
                annots = self.plot_annots(annots, m11, ul, l11, xlim)

            hands.append(h4)
            labs.append("Additional inserted planets" if len(r11) > 1 else "Additional inserted planet")

        if len(r12) > 0:
            if self.config_parameters["use_mass"] == "False":
                h5 = plt.scatter(r12, np.ones(len(r12))*ul/10, c=color_scheme[4], marker="^", edgecolors="k", s=200, linewidth=2, zorder=2)
                annots = self.plot_annots(annots, r12, ul, l12, xlim)

            else:
                h5 = plt.scatter(m12, np.ones(len(m12))*ul/10, c=color_scheme[4], marker="^", edgecolors="k", s=200, linewidth=2, zorder=2)
                annots = self.plot_annots(annots, m12, ul, l12, xlim)

            hands.append(h5)

            if system == "Inner Solar System":
                labs.append("Gas giant planets" if len(r12) > 1 else "Gas giant planet")

            else:
                labs.append("Known non-transiting planets" if len(r12) > 1 else "Known non-transiting planet")


        if len(r2) > 0:
            if self.config_parameters["use_mass"] == "False":
                h6 = plt.scatter(r2, np.ones(len(r2))*ul/10, c="white", marker="X", s=200, edgecolors="k", linewidth=2, zorder=2)
                annots = self.plot_annots(annots, r2, ul, l2, xlim)

            else:
                h6 = plt.scatter(m2, np.ones(len(m2))*ul/10, c="white", marker="X", s=200, edgecolors="k", linewidth=2, zorder=2)
                annots = self.plot_annots(annots, m2, ul, l2, xlim)

            hands.append(h6)
            labs.append("Known planets removed from system" if len(r2) > 1 else "Known planet removed from system")


        if len(r3) > 0:
            if self.config_parameters["use_mass"] == "False":
                h7 = plt.scatter(r3, np.ones(len(r3))*ul/10, c=color_scheme[5], marker='d', s=200, edgecolors="k", linewidth=2, zorder=2)
                annots = self.plot_annots(annots, r3, ul, l3, xlim)

            else:
                h7 = plt.scatter(m3, np.ones(len(m3))*ul/10, c=color_scheme[5], marker='d', s=200, edgecolors="k", linewidth=2, zorder=2)
                annots = self.plot_annots(annots, m3, ul, l3, xlim)

            hands.append(h7)
            labs.append("Unconfirmed planet candidates" if len(r3) > 1 else "Unconfirmed planet candidate")

        plt.xlim(xlim[0], xlim[1])
        plt.ylim(0,ul)
        plt.ylabel("Relative Likelihood", fontsize=32)
        plt.legend(hands, labs, loc=9, fontsize=14, ncol=2)
        ax.tick_params(labelsize=24)
        fig.suptitle((r"$\tau$ Ceti" if system == "tau Cet" else system) + " Planet " + ("Radius " if self.config_parameters["use_mass"] == "False" else "Mass ") + "Relative Likelihood", fontsize=30)
        plt.savefig(self.sdir + plt_savename)

        if self.config_parameters["show_plots"] == "True":
            plt.show()

        elif self.config_parameters["show_plots"] == "False":
            plt.close()

        fig, ax = plt.subplots(figsize=(12,8))

        if self.config_parameters["ind_i"] == "full":
            ik = np.array(ik)
            plt_savename = system.replace(" ", "") + "_inc_full_" + self.config_parameters["period"] + "_" + self.config_parameters["mass_radius"] + ("_" + self.config_parameters["otegi_rho"] if self.config_parameters["mass_radius"] == "otegi" else "") + "_" + self.config_parameters["inclination"] + "_" + self.config_parameters["eccentricity"] + "_" + self.config_parameters["stability"] + ".png"

        elif self.config_parameters["ind_i"] == "truncated":
            ik = np.array([ik[j] if ik[j] < 90 else 180-ik[j] for j in range(len(ik))])     # toggle for allowing inclinations to be greater than 90 degrees or truncated back to 0-90 range
            plt_savename = system.replace(" ", "") + "_inc_trunc_" + self.config_parameters["period"] + "_" + self.config_parameters["mass_radius"] + ("_" + self.config_parameters["otegi_rho"] if self.config_parameters["mass_radius"] == "otegi" else "") + "_" + self.config_parameters["inclination"] + "_" + self.config_parameters["eccentricity"] + "_" + self.config_parameters["stability"] + ".png"

        n, _, h1 = plt.hist(ik, bins=np.linspace(0, 180.5, 362), weights=np.ones(len(ik)) / (len(ik)), color=color_scheme[0])
        hands = [h1[0]]
        labs = [(r"$\tau$ Ceti" if system == "tau Cet" else system) + " DYNAMITE Predictions"]
        ul = 1.5*np.amax(n)
        annots = []
        j1 = [j for j in i1 if j != "?"]
        j11 = [j for j in i11 if j != "?"]
        j12 = [j for j in i12 if j != "?"]
        j2 = [j for j in i2 if j != "?"]
        j3 = [j for j in i3 if j != "?"]
        xval = [j1]

        if len(j11) > 0:
            xval.append(j11)

        if len(j12) > 0:
            xval.append(j12)

        if len(j2) > 0:
            xval.append(j2)

        if len(j3) > 0:
            xval.append(j3)

        if len(np.hstack(xval)) > 0:
            xlim = [math.floor(min(np.hstack(xval))) - 10, math.ceil(max(np.hstack(xval))) + 10]

        else:
            xlim = [0, 180]

        if self.config_parameters["plt_PDFs"] == "True":
            h2 = plt.plot(inc, Pi/2, color=color_scheme[1], linewidth=3)
            hands.append(h2[0])
            labs.append("Probability Density Function")

        if len(i1) > 0:
            h3 = plt.scatter(i1, np.ones(len(i1))*ul/10, c=color_scheme[4], edgecolors="k", s=200, linewidth=2, zorder=2)
            hands.append(h3)

            if system == "Inner Solar System":
                labs.append("Terrestrial planets" if len(i1) > 1 else "Terrestrial planet")

            else:
                labs.append("Known transiting planets" if len(i1) > 1 else "Known transiting planet")

            annots = self.plot_annots(annots, i1, ul, l1, xlim)

        if len(j11) > 0:
            h4 = plt.scatter(j11, np.ones(len(j11))*ul/10, c=color_scheme[3], marker="s", edgecolors="k", s=200, linewidth=2, zorder=2)
            hands.append(h4)
            labs.append("Additional inserted planets" if len(i11) > 1 else "Additional inserted planet")
            annots = self.plot_annots(annots, j11, ul, [l11[x] for x in range(len(l11)) if i11[x] != "?"], xlim)
        
        if len(j12) > 0:
            h5 = plt.scatter(j12, np.ones(len(j12))*ul/10, c=color_scheme[4], marker="^", edgecolors="k", s=200, linewidth=2, zorder=2)
            hands.append(h5)

            if system == "Inner Solar System":
                labs.append("Gas giant planets" if len(i12) > 1 else "Gas giant planet")

            else:
                labs.append("Known non-transiting planets" if len(i12) > 1 else "Known non-transiting planet")

            annots = self.plot_annots(annots, j12, ul, [l12[x] for x in range(len(l12)) if i12[x] != "?"], xlim)

        if len(j2) > 0:
            h6 = plt.scatter(j2, np.ones(len(j2))*ul/10, c="white", marker="X", s=200, edgecolors="k", linewidth=2, zorder=2)
            hands.append(h6)
            labs.append("Known planets removed from system" if len(i2) > 1 else "Known planet removed from system")
            annots = self.plot_annots(annots, j2, ul, l2, xlim)

        if len(j3) > 0:
            h7 = plt.scatter(j3, np.ones(len(j3))*ul/10, c=color_scheme[5], marker='d', s=200, edgecolors="k", linewidth=2, zorder=2)
            hands.append(h7)
            labs.append("Unconfirmed planet candidates" if len(i3) > 1 else "Unconfirmed planet candidate")
            annots = self.plot_annots(annots, j3, ul, l3, xlim)

        plt.xlabel("Inclination (degrees)", fontsize=24)
        plt.ylabel("Relative Likelihood", fontsize=24)
        plt.xlim(xlim[0], xlim[1])
        plt.ylim(0, np.amax(n)*1.5)
        ax.tick_params(labelsize=18)
        plt.legend(hands, labs, loc=9, fontsize=14, ncol=2)
        fig.suptitle((r"$\tau$ Ceti" if system == "tau Cet" else system) + " Inclination Relative Likelihood", fontsize=30)
        plt.savefig(self.sdir + plt_savename)

        if self.config_parameters["show_plots"] == "True":
            plt.show()

        elif self.config_parameters["show_plots"] == "False":
            plt.close()

        fig, ax = plt.subplots(figsize=(12,8))
        ek = np.array(ek)
        plt_savename = system.replace(" ", "") + "_ecc_" + self.config_parameters["period"] + "_" + self.config_parameters["mass_radius"] + ("_" + self.config_parameters["otegi_rho"] if self.config_parameters["mass_radius"] == "otegi" else "") + "_" + self.config_parameters["inclination"] + "_" + self.config_parameters["eccentricity"] + "_" + self.config_parameters["stability"] + ".png"
        n, _, h1 = plt.hist(ek, bins=np.linspace(0, 1, 201), weights=np.ones(len(ek)) / (len(ek)), color=color_scheme[0])
        hands = [h1[0]]
        labs = [(r"$\tau$ Ceti" if system == "tau Cet" else system) + " DYNAMITE Predictions"]
        ul = 1.5*np.amax(n)
        annots = []
        xlim = [-0.05, max(np.hstack([[j for j in e1 if j != "?"], [j for j in e2 if j != "?"], [j for j in e12 if j != "?"], [j for j in e2 if j != "?"], [j for j in e3 if j != "?"], [np.percentile(ek, 95)]])) + 0.05]

        if self.config_parameters["plt_PDFs"] == "True":
            h2 = plt.plot(ecc, 5*Pe/sum(Pe), color=color_scheme[1], linewidth=3)
            hands.append(h2[0])
            labs.append("Probability Density Function")

        if len([j for j in e1 if j != "?"]) > 0:
            h3 = plt.scatter([j for j in e1 if j != "?"], np.ones(len([j for j in e1 if j != "?"]))*ul/10, c=color_scheme[4], edgecolors="k", s=200, linewidth=2, zorder=2)
            hands.append(h3)

            if system == "Inner Solar System":
                labs.append("Terrestrial planets" if len(e1) > 1 else "Terrestrial planet")

            else:
                labs.append("Known transiting planets" if len(e1) > 1 else "Known transiting planet")

            annots = self.plot_annots(annots, [j for j in e1 if j != "?"], ul, l1, xlim)

        if len([j for j in e11 if j != "?"]) > 0:
            h4 = plt.scatter([j for j in e11 if j != "?"], np.ones(len([j for j in e11 if j != "?"]))*ul/10, c=color_scheme[3], marker="s", edgecolors="k", s=200, linewidth=2, zorder=2)
            hands.append(h4)
            labs.append("Additional inserted planets" if len(e11) > 1 else "Additional inserted planet")
            annots = self.plot_annots(annots, [j for j in e11 if j != "?"], ul, l11, xlim)
        
        if len([j for j in e12 if j != "?"]) > 0:
            h5 = plt.scatter([j for j in e12 if j != "?"], np.ones(len([j for j in e12 if j != "?"]))*ul/10, c=color_scheme[4], marker="^", edgecolors="k", s=200, linewidth=2, zorder=2)
            hands.append(h5)

            if system == "Inner Solar System":
                labs.append("Gas giant planets" if len(e12) > 1 else "Gas giant planet")

            else:
                labs.append("Known non-transiting planets" if len(e12) > 1 else "Known non-transiting planet")

            annots = self.plot_annots(annots, [j for j in e12 if j != "?"], ul, [l12[x] for x in range(len(l12)) if l12[x] != "?"], xlim)

        if len([j for j in e2 if j != "?"]) > 0:
            h6 = plt.scatter([j for j in e2 if j != "?"], np.ones(len([j for j in e2 if j != "?"]))*ul/10, c="white", marker="X", s=200, edgecolors="k", linewidth=2, zorder=2)
            hands.append(h6)
            labs.append("Known planets removed from system" if len(e2) > 1 else "Known planet removed from system")
            annots = self.plot_annots(annots, e2, ul, l2, xlim)

        if len([j for j in e3 if j != "?"]) > 0:
            h7 = plt.scatter([j for j in e3 if j != "?"], np.ones(len([j for j in e3 if j != "?"]))*ul/10, c=color_scheme[5], marker='d', s=200, edgecolors="k", linewidth=2, zorder=2)
            hands.append(h7)
            labs.append("Unconfirmed planet candidates" if len(e3) > 1 else "Unconfirmed planet candidate")
            annots = self.plot_annots(annots, e3, ul, l3, xlim)

        plt.xlabel("Eccentricity", fontsize=24)
        plt.ylabel("Relative Likelihood", fontsize=24)
        plt.xlim(xlim[0], xlim[1])
        plt.ylim(0, np.amax(n)*1.5)
        ax.tick_params(labelsize=18)
        plt.legend(hands, labs, loc=9, fontsize=14, ncol=2)
        fig.suptitle((r"$\tau$ Ceti" if system == "tau Cet" else system) + " Eccentricity Relative Likelihood", fontsize=30)
        plt.savefig(self.sdir + plt_savename)

        if self.config_parameters["show_plots"] == "True":
            plt.show()

        elif self.config_parameters["show_plots"] == "False":
            plt.close()



    def mr_appends(self, pred, val):
        """Finds the mass from the radius or vice versa"""

        if self.config_parameters["mass_radius"] == "otegi":
            ret = self.otegi_mr(val, pred)

        return ret



    def set_up(self, targets_dict_tuple, target):
        """Sets up target"""

        targets_dict = {}

        for k in targets_dict_tuple.keys():
            tx = []
            td = targets_dict_tuple.get(k)

            for i in range(len(td)):
                txx = []

                for j in range(len(td[i])):
                    txx.append(td[i][j][0]) if isinstance(td[i][j], tuple) else txx.append(td[i][j])

                    if i > 0 and j == 2:
                        txx.append(td[i][j][3])

                tx.append(list(txx))

            targets_dict[k] = tuple(tx)
           
        t = list(targets_dict[target])

        if t[0][2] == None:
            t[0][2] = 0.0

        if t[0][4] == None:
            t[0][4] = 0.0

        for x in range(len(t)):
            for y in range(len(t[x])):
                if x > 0 and not isinstance(t[x][y], str) and y < len(t[x]) - 2:
                    t[x][y] = float(t[x][y])

                elif y == 1 and t[x][y] == "?":
                    t[x][y] = float(self.mr_convert(t[x][y+1], "radius"))

                elif y == 2 and t[x][y] == "?":
                    t[x][y] = float(self.mr_convert(t[x][y-1], "mass"))

                elif y == len(t[x]) - 2 and t[x][y] != "?":
                    t[x][y] = int(t[x][y])

        return t[0][0], t[0][1], t[0][3], t[0][5], [t[i][:-1] for i in range(1, len(t))], np.array([t[i][-1] for i in range(1, len(t))]), t[0][5]



    def plot_annots(self, annots, xl, ul, annot_l, xlim, scale="linear", ax=None):
        """Puts planet annotation labels on scatter plots"""

        for i in range(len(xl)):
            xy1 = [xl[i], ul/10 - ul/15]

            if scale == "linear":
                overlap = [annots[j] for j in range(len(annots)) if xy1[0] != "?" and abs(xy1[0] - annots[j][0]) < (xlim[1] - xlim[0])/35]

            else:
                overlap = [annots[j] for j in range(len(annots)) if xy1[0] != "?" and xy1[0] / annots[j][0] < 1.2 and xy1[0] / annots[j][0] > 0.8]

            if len(overlap) == 1:
                xy1[1] += (8*ul/75 if overlap[0][1] == ul/10 - ul/15 else 0)

            elif len(overlap) > 1:
                xy1[1] += 8*ul/75 + (len(overlap) - 1)*ul/20

            if annot_l[i].find("PxP") != -1:
                xy1[1] += 8*ul/75 + (len(overlap))*ul/20

            if ax == None:
                xy = plt.annotate(annot_l[i], (xl[i],ul/10), color="k", xytext=xy1, ha="center", weight='bold', fontsize=24)
            else:
                xy = ax.annotate(annot_l[i], (xl[i], ul/10), color="k", xytext=xy1, ha="center", weight='bold', fontsize=24)

            annots.append(list(xy.xyann))

        return annots




    def otegi_mr(self, measurement, predict, force=None):
        """Uses a power-law MR to predict mass in Earth values from radius in Earth values or vice versa"""

        if predict == "radius":
            R_r = 1.03*measurement**0.29
            R_v = 0.7*measurement**0.63
            R_j = 14.357*(measurement/131.61)**-0.04
            
            if force == "rocky" or measurement < 5:
                return R_r

            if force == "volatile" or (measurement < 131.61 and (measurement > 25 or (measurement > 5 and self.config_parameters["otegi_rho"] == "volatile"))):
                return R_v

            return R_j

        elif predict == "mass":
            M_r = 0.9*measurement**3.45
            M_v = 1.74*measurement**1.58
            M_j = 131.61*(measurement/14.357)**(-25)

            if force == "rocky" or M_r < 5:
                return M_r

            if force == "volatile" or (M_v < 131.61 and (M_r > 25 or (M_r > 5 and self.config_parameters["otegi_rho"] == "volatile"))):
                return M_v

            return M_j


    # DEPRECATED: def mrexo_masses(self, Rk, i):
        """Calculates masses from radius distribution using MRExo"""

        # return pfm(Rk[i], predict='mass', dataset='kepler')[0]



    def K3(self, P, M):
        """Calculates semi-major axis in cm using period in days and mass in solar masses"""
        
        return (self.G*(M*self.M_sun)*(P*self.seconds_per_day)**2/(4*math.pi**2))**(1/3)



    def K3_inv(self, a, M):
        """Calculates period in days using semi-major axis in au and mass in solar masses"""
        
        return (4*math.pi**2*(a*self.au)**3/(self.G*(M*self.M_sun)))**0.5/self.seconds_per_day



    def mr_convert(self, meas, pred):
        """Runs conversion from mass to radius and vice versa"""

        #if self.config_parameters["mass_radius"] == "mrexo":
        #    try:
        #        return pfm(meas, predict=pred, dataset="kepler")[0]

        #    except:
        #        return 1e-10

        if self.config_parameters["mass_radius"] == "otegi":
            return self.otegi_mr(meas, pred)



if __name__ == '__main__':
    if len(sys.argv) > 1:
        dynamite_plots(cfname=sys.argv[1])

    else:
        dynamite_plots()
