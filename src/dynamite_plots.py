### DYNAmical Multi-planet Injection TEster (DYNAMITE) ###
### Plotting File ###
### Jeremy Dietrich ###
### 2020 July 15 ###
### Version 1.2 ###
### Dietrich & Apai (2020), Astronomical Journal in press ###
### http://arxiv.org/pdf/2007.06521.pdf ###

import math
import numpy as np
import scipy.stats as spst
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patches as mpatch
from mrexo import predict_from_measurement as pfm

class dynamite_plots:

    def __init__(self, Pk=[], P=[], PP=[], per=[], Rk=[], R=[], PR=[], ik=[], inc=[], Pi=[], deltas=[], ratios=[], tdm=[], tdle=[], tdue=[], tpm=[], tple=[], tpue=[], targets=[], pers=[], rads=[]):
        """Sets up plotting routines"""

        self.config_parameters = {}

        try:
            config_data = np.loadtxt("dynamite_config.txt", dtype=str, delimiter='::')

        except IOError:
            print("Error, configuration file not found!")
            exit()

        for i in range(len(config_data)):
            self.config_parameters[config_data[i, 0]] = config_data[i, 1]

        if len(Pk) == 0:
            Pk, P, PP, per, Rk, R, PR, ik, inc, Pi, deltas, ratios, tdm, tdle, tdue, tpm, tple, tpue, targets, pers, rads = self.read_saved_data()

        if self.config_parameters["mode"] != "single":
            P = P[0]
            R = R[0]
            inc = inc[0]

        if self.config_parameters["plt_P_R"] == "True":
            self.plot_P_R(Pk, P, Rk, R, targets, pers, rads)

        if self.config_parameters["plt_tdtp"] == "True":
            self.plot_td_tp(tdm, tdle, tdue, tpm, tple, tpue, targets)

        if self.config_parameters["plt_deltas"] == "True":
            self.plot_deltas(deltas)

        if self.config_parameters["plt_ratios"] == "True":
            self.plot_ratios(ratios)

        if self.config_parameters["plt_indpars"] == "True":
            self.plot_ind_params(Pk, P, PP, per, Rk, R, PR, ik, inc, Pi)

        print(datetime.now(), "Finishing Plots")



    def read_saved_data(self):
        """Reads in the saved data"""

        print(datetime.now(), "Reading in Saved Data")

        with np.load("saved_data.npz", allow_pickle=True) as fd:
            Pk, P, PP, per, Rk, R, PR, ik, inc, Pi, deltas, ratios, tdm, tdle, tdue, tpm, tple, tpue, targets, pers, rads = fd["data"]
        
        return Pk, P, PP, per, Rk, R, PR, ik, inc, Pi, deltas, ratios, tdm, tdle, tdue, tpm, tple, tpue, targets, pers, rads



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
            np.savetxt("plot_fig_p.txt", plot_fig_p)
            np.savetxt("plot_fig_r.txt", plot_fig_r)

        elif self.config_parameters["saved"] == "True":
            plot_fig_p = np.loadtxt("plot_fig_p.txt")
            plot_fig_r = np.loadtxt("plot_fig_r.txt")

        for i in range(len(targets)):
            targets[i][1] = int(float(targets[i][1])*100 - 10)
            targets[i][2] = int(float(targets[i][2])*100)
            targets[i][3] = int(float(targets[i][3])*100)
            targets[i][4] = int(float(targets[i][4])*1000 - 100)
            targets[i][5] = int(float(targets[i][5])*1000)
            targets[i][6] = int(float(targets[i][6])*1000)

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
        ax.set_xscale("Log")
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

        if self.config_parameters["mode"] != "Kepler":
            for i in range(len(targets)):
                if targets[i][0] in bottom:
                    plt.annotate(targets[i][0][targets[i][0].find(" "):], (int(targets[i][1]), int(targets[i][4])), color="w", textcoords="offset points", xytext=(0,-20), ha="center", weight='bold', fontsize=16)

                elif targets[i][0] in right:
                    plt.annotate(targets[i][0][targets[i][0].find(" "):], (int(targets[i][1]), int(targets[i][4])), color="w", textcoords="offset points", xytext=(25,-7), ha="center", weight='bold', fontsize=16)

                elif targets[i][0] in left:
                    plt.annotate(targets[i][0][targets[i][0].find(" "):], (int(targets[i][1]), int(targets[i][4])), color="w", textcoords="offset points", xytext=(-30,-7), ha="center", weight='bold', fontsize=16)

                elif targets[i][0] in tl:
                    plt.annotate(targets[i][0][targets[i][0].find(" "):], (int(targets[i][1]), int(targets[i][4])), color="w", textcoords="offset points", xytext=(-25,7), ha="center", weight='bold', fontsize=16)

                elif targets[i][0] in tr:
                    plt.annotate(targets[i][0][targets[i][0].find(" "):], (int(targets[i][1]), int(targets[i][4])), color="w", textcoords="offset points", xytext=(12,5), ha="center", weight='bold', fontsize=16)

                elif targets[i][0] in bl:
                    plt.annotate(targets[i][0][targets[i][0].find(" "):], (int(targets[i][1]), int(targets[i][4])), color="w", textcoords="offset points", xytext=(-15,-20), ha="center", weight='bold', fontsize=16)

                elif targets[i][0] in br:
                    plt.annotate(targets[i][0][targets[i][0].find(" "):], (int(targets[i][1]), int(targets[i][4])), color="w", textcoords="offset points", xytext=(15,-20), ha="center", weight='bold', fontsize=16)

                else:
                    plt.annotate(targets[i][0][targets[i][0].find(" "):], (int(targets[i][1]), int(targets[i][4])), color="w", textcoords="offset points", xytext=(0,8), ha="center", weight='bold', fontsize=16)

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

        plt.savefig("PRfig_" + self.config_parameters["mode"] + "_" + self.config_parameters["period"] + "_ell.png", bbox_inches='tight')

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
            ax.set_xscale("Log")
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
        plt.scatter(np.where(np.isclose(Pf, pers[0][0], atol=(0.005 if pers[0][0] < 1 else 0.05)))[0][0], pc, color="r", s=int(round(rads[0][0]*20, 0)), label="Known Planets")

        for i in range(len(pers)):
            for j in range(len(pers[i])):
                plt.scatter(np.where(np.isclose(Pf, pers[i][j], atol=0.006))[0][0], pc, color="r", s=int(round(rads[i][j]*20, 0)))

            pc += 200

        plt.savefig("logPfig_" + self.config_parameters["mode"] + "_" + self.config_parameters["period"] + ".png", bbox_inches='tight')

        if self.config_parameters["show_plots"] == "True":
            plt.show()

        elif self.config_parameters["show_plots"] == "False":
            plt.close()
       
        fig, ax = plt.subplots(1,1, figsize=(10,12))
        fig.suptitle("Planet Radius Relative Likelihoods for " + self.config_parameters["mode"] + " Systems", fontsize=30)
        img = ax.imshow(pfri, cmap=plt.cm.Blues, origin="lower", aspect="auto")
        xlabels = [0, 1, 2, 3, 4, 5]
        ylabels = [targets[i][0] for i in range(len(targets))]
        ax.set_xticks([-100, 900, 1900, 2900, 3900, 4900]) 
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
        plt.scatter(np.where(np.isclose(Rf, rads[0][0], atol=0.006))[0][0], rc, color="r", label="Known Planets")

        for i in range(len(rads)):
            for j in range(len(rads[i])):
                plt.scatter(np.where(np.isclose(Rf, rads[i][j], atol=0.006))[0][0], rc, color="r")

            rc += 200

        plt.savefig("Rfig_" + self.config_parameters["mode"] + ".png", bbox_inches='tight')

        if self.config_parameters["show_plots"] == "True":
            plt.show()

        elif self.config_parameters["show_plots"] == "False":
            plt.close()



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
        plt.savefig("Deltas_" + self.config_parameters["mode"] + "_" + self.config_parameters["period"] + ".png")

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
        plt.xscale("Log")
        plt.xlabel(r"Period Ratio $P_i/P_{i-1}$", fontsize=20)
        plt.ylabel("Occurrence", fontsize=20)
        plt.legend(fontsize=16)
        ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
        ax.set_xticklabels([0.1,1,10])
        ax.tick_params(labelsize=14)
        ax.set_yticks([0, 5, 10, 15, 20])
        fig.suptitle(r"Subsample vs $Kepler$ Period Ratios", fontsize=30)
        plt.xlim(1, 10**1.4)
        plt.savefig(self.config_parameters["mode"] + "_Kepler_period_ratios.png")

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
        #bunch = ["TOI 261", "TOI 266", "TOI 396", "TOI 411", "TOI 487", "TOI 561", "TOI 703", "TOI 797", "TOI 1238", "TOI 1346", "TOI 1453", "TOI 1469", "TOI 1726"]
        bunch = ["TOI 256", "TOI 714", "TOI 736"]
        #bunch = ["TOI 396", "TOI 411", "TOI 286", "TOI 1469", "TOI 487", "TOI 174", "TOI 261", "TOI 1453", "TOI 713", "TOI 1339", "TOI 431", "TOI 282", "TOI 1346", "TOI 1238", "TOI 266", "TOI 1726", "TOI 797", "TOI 1269", "TOI 703", "TOI 1730", "TOI 696", "TOI 836", "TOI 732", "TOI 1449", "TOI 763", "TOI 1260"]
        bottom = []
        left = []
        right = []

        for i in range(len(names)):
            tdle1 = tdm[i] - tdle[i]
            tdue1 = tdm[i] + tdue[i]
            tple1 = tpm[i] - tple[i]
            tpue1 = tpm[i] + tpue[i]
            ellipse = mpatch.Ellipse((np.mean([tdle1, tdue1]), np.mean([tple1, tpue1])), tdue1-tdle1, tpue1-tple1, alpha=0.1)
            ax.add_patch(ellipse)

            if names[i] not in bunch:
                if names[i] in bottom:
                    plt.annotate(names[i][names[i].find(" "):], (tdm[i], tpm[i]), textcoords="offset points", xytext=(0,-20), ha='center', fontsize=16)

                elif names[i] in left:
                    plt.annotate(names[i][names[i].find(" "):], (tdm[i], tpm[i]), textcoords="offset points", xytext=(-30,-5), ha='center', fontsize=16)

                elif names[i] in right:
                    plt.annotate(names[i][names[i].find(" "):], (tdm[i], tpm[i]), textcoords="offset points", xytext=(25,-5), ha='center', fontsize=16)

                else:
                    plt.annotate(names[i][names[i].find(" "):], (tdm[i], tpm[i]), textcoords="offset points", xytext=(0,10), ha='center', fontsize=16)

            else:
                if names[i] in bottom:
                    plt.annotate(names[i][names[i].find(" "):], (tdm[i], tpm[i]), textcoords="offset points", xytext=(0,-20), ha='center', fontsize=16)

                elif names[i] in left:
                    plt.annotate(names[i][names[i].find(" "):], (tdm[i], tpm[i]), textcoords="offset points", xytext=(-30,-5), ha='center', fontsize=16)

                elif names[i] in right:
                    plt.annotate(names[i][names[i].find(" "):], (tdm[i], tpm[i]), textcoords="offset points", xytext=(30,-5), ha='center', fontsize=16)

                else:
                    plt.annotate(names[i][names[i].find(" "):], (tdm[i], tpm[i]), textcoords="offset points", xytext=(0,10), ha='center', fontsize=16)

        ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
        plt.xlim(150, max(tdm) + max(tdue))
        #plt.xlim(150, 3500)
        #plt.xlim(150, 800)
        #plt.ylim(0.75, 0.85)
        #plt.xlim(150, 1150)
        #plt.ylim(0.8, 0.89)
        ax.tick_params(labelsize=14)
        plt.xlabel("Transit Depth (ppm)", fontsize=20)
        plt.ylabel("Transit Probability", fontsize=20)
        fig.suptitle("Transit Probability vs Transit Depth", fontsize=30)
        plt.savefig(self.config_parameters["mode"] + "_td_tp.png")

        if self.config_parameters["show_plots"] == "True":
            plt.show()

        elif self.config_parameters["show_plots"] == "False":
            plt.close()



    def plot_ind_params(self, Pk, P, PP, per, Rk, R, PR, ik, inc, Pi):
        """Plots histograms of each individual distribution"""

        print(datetime.now(), "Creating Histograms for", self.config_parameters["system"])

        system = self.config_parameters["system"]
        p1 = []     # known planet period
        p11 = []    # inserted planet period for iteration
        p12 = [20, 49.41, 162.87, 636.13]       # known non-transiting planet periods
        p2 = []     # removed known planet period
        p3 = [13.965, 35.362, 94.11]        # unconfirmed planet candidate periods
        r1 = []     # " for planet radius
        r11 = []    # " for planet radius
        r12 = [self.mr_predict(1.75/math.sin(35*math.pi/180), "Radius"), self.mr_predict(1.83/math.sin(35*math.pi/180), "Radius"), self.mr_predict(3.93/math.sin(35*math.pi/180), "Radius"), self.mr_predict(3.93/math.sin(35*math.pi/180), "Radius")]      # " for planet radius
        r2 = []     # " for planet radius
        r3 = [self.mr_predict(2/math.sin(35*math.pi/180), "Radius"), self.mr_predict(3.1/math.sin(35*math.pi/180), "Radius"), self.mr_predict(3.6/math.sin(35*math.pi/180), "Radius")]     # " for planet radius
        i1 = []
        i11 = []
        i12 = [35, 35, 35, 35]
        i2 = []
        i3 = [35, 35, 35]
        l1 = []
        l11 = []
        l12 = ["g", "h", "e", "f"]
        l2 = []
        l3 = ["b", "c", "d"]

        if self.config_parameters["ind_P"] == "linear_zoom":
            fig, ax = plt.subplots(figsize=(12, 8))
            plt.hist(Pk, bins=np.arange(0.5, round(5*max(per)) + 1, 0.5), weights=np.ones(len(Pk))*2 / (len(Pk)), label="DYNAMITE Predictions")

            if self.config_parameters["plt_PDFs"] == "True":
                plt.plot(P, PP, label="PDF")

            plt.xlabel("Period (days)")
            plt.ylabel("Relative Likelihood")
            plt.xlim(0, 5*max(per))
            plt.legend()
            fig.suptitle((r"$\tau$ Ceti" if system == "tau Ceti" else system) + " Period Relative Likelihood", fontsize=30)

            if self.config_parameters["period"] == "epos":
                plt.savefig(system + "_P_linear_zoom_epos.png")

            elif self.config_parameters["period"] == "syssim":
                plt.savefig(system + "_P_linear_zoom_syssim.png")

            if self.config_parameters["show_plots"] == "True":
                plt.show()

            elif self.config_parameters["show_plots"] == "False":
                plt.close()

        elif self.config_parameters["ind_P"] == "linear":
            Pb = list(P)
            Pb.append(730.1)
            fig, ax = plt.subplots(figsize=(12, 8))
            plt.hist(Pk, bins=Pb, weights=np.ones(len(Pk))*10 / len(Pk), label="DYNAMITE Predictions")

            if self.config_parameters["plt_PDFs"] == "True":
                plt.plot(P, PP, label="PDF")

            plt.xlabel("Period (days)")
            plt.ylabel("Relative Likelihood")
            plt.xlim(0, 5*max(per))
            plt.legend()
            fig.suptitle((r"$\tau$ Ceti" if system == "tau Ceti" else system) + (" Period Ratio" if self.config_parameters["period"] == "epos" else " Clustered Periods") + " Relative Likelihood", fontsize=30)

            if self.config_parameters["period"] == "epos":
                plt.savefig(system + "_P_linear_epos.png")

            elif self.config_parameters["period"] == "syssim":
                plt.savefig(system + "_P_linear_syssim.png")

            if self.config_parameters["show_plots"] == "True":
                plt.show()

            elif self.config_parameters["show_plots"] == "False":
                plt.close()

        elif self.config_parameters["ind_P"] == "log":
            Pl = np.logspace(-0.3, 2.86, 317)
            Plb = list(Pl)
            Plb.append(10**2.895)
            bins = np.array(Plb)
            widths = (bins[1:] - bins[:-1])
            hist, _ = np.histogram(Pk, bins=bins)
            hist_norm = hist/widths
            fig, ax = plt.subplots(figsize=(12,8))
            plt.bar(bins[:-1], hist_norm/1e5, widths, label="DYNAMITE Predictions")
            PPl = np.interp(Pl, P, PP)

            if self.config_parameters["plt_PDFs"] == "True":
                plt.plot(Pl, PPl, color='#ff7f0e', label="PDF", linewidth=4)
            
            if len(p1) > 0:
                plt.scatter(p1, np.ones(len(p1))*np.amax(hist_norm)/1e6, c="w", edgecolors="g", s=[r1[i]*200 for i in range(len(r1))], linewidth=2, label=("Known planets in system" if len(p1) > 1 else "Known planet in system"), zorder=2)

                for i in range(len(l1)):
                    plt.annotate(l1[i], (p1[i], np.amax(hist_norm)/1e6), color="k", textcoords="offset points", ha="center", weight='bold', fontsize=16)

            if len(p11) > 0:
                plt.scatter(p11, np.ones(len(p11))*np.amax(hist_norm)/1e6, c="w", edgecolors="purple", marker="s", s=[r11[i]*200 for i in range(len(r11))], linewidth=2, label="Highest relative likelihood for inserted planet", zorder=2)

                for i in range(len(l11)):
                    plt.annotate(l11[i], (p11[i], np.amax(hist_norm)/1e6), color="k", textcoords="offset points", ha="center", weight='bold', fontsize=16)

            if len(p12) > 0:
                plt.scatter(p12, np.ones(len(p12))*np.amax(hist_norm)/1e6, c="w", edgecolors="g", marker="^", s=[r12[i]*200 for i in range(len(r12))], linewidth=2, label=("Known non-transiting planets" if len(p12) > 1 else "Known non-transiting planet"), zorder=2)

                for i in range(len(l12)):
                    plt.annotate(l12[i], (p12[i], np.amax(hist_norm)/1e6), color="k", textcoords="offset points", ha="center", weight='bold', fontsize=16)

            if len(p2) > 0:
                plt.scatter(p2, np.ones(len(p2))*np.amax(hist_norm)/1e6, c="w", edgecolors="r", marker="X", s=[r2[i]*200 for i in range(len(r2))], linewidth=2, label=("Known planets removed from system" if len(p2) > 1 else "Known planet removed from system"), zorder=2)

                for i in range(len(l2)):
                    plt.annotate(l2[i], (p2[i], np.amax(hist_norm)/1e6), color="k", textcoords="offset points", ha="center", weight='bold', fontsize=16)

            if len(p3) > 0:
                plt.scatter(p3, np.ones(len(p3))*np.amax(hist_norm)/1e6, c="w", edgecolors="y", marker='$?$', s=[r3[i]*200 for i in range(len(r3))], linewidth=2, label=("Unconfirmed planet candidates" if len(p3) > 1 else "Unconfirmed planet candidate"), zorder=2)

                for i in range(len(l3)):
                    plt.annotate(l3[i], (p3[i], np.amax(hist_norm)/1e6), color="k", textcoords="offset points", ha="center", weight='bold', fontsize=16)

            plt.xlabel("Log Period (days)", fontsize=20)
            plt.ylabel("Relative Likelihood", fontsize=20)
            plt.xscale("Log")
            plt.xlim(0.5, 730)
            plt.ylim(0, np.amax(hist_norm)*1.2/1e5)
            plt.legend(fontsize=16)
            ax.tick_params(labelsize=14)
            ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
            fig.suptitle((r"$\tau$ Ceti" if system == "tau Ceti" else system) + " Period Relative Likelihood", fontsize=30)

            if self.config_parameters["period"] == "epos":
                plt.savefig(system + "_P_log_epos.png")

            elif self.config_parameters["period"] == "syssim":
                plt.savefig(system + "_P_log_syssim.png")

            if self.config_parameters["show_plots"] == "True":
                plt.show()

            elif self.config_parameters["show_plots"] == "False":
                plt.close()

        if self.config_parameters["ind_R"] == "linear_zoom":
            Rk = np.array(Rk)
            fig, ax = plt.subplots(figsize=(12,8))
            plt.hist(Rk, bins=np.arange(R[0], R[-1] + 0.01, 0.005), weights=np.ones(len(Rk))*20 / len(Rk), label="DYNAMITE Predictions")

            if self.config_parameters["plt_PDFs"] == "True":
                plt.plot(R, PR, label="PDF")

            plt.xlabel(r"Radius ($R_{\oplus}$)")
            plt.ylabel("Relative Likelihood")
            plt.legend()
            fig.suptitle((r"$\tau$ Ceti" if system == "tau Ceti" else system) + " Planet Radius Relative Likelihood")
            plt.savefig(system + "_R_linear_zoom.png")

            if self.config_parameters["show_plots"] == "True":
                plt.show()

            elif self.config_parameters["show_plots"] == "False":
                plt.close()

        elif self.config_parameters["ind_R"] == "linear":
            Rb = list(R)
            Rb.append(float(self.config_parameters["radmax"]) + 0.01)
            fig, ax = plt.subplots(figsize=(12,8))
            n, _, _ = plt.hist(Rk, bins=Rb, weights=np.ones(len(Rk))*100 / len(Rk), label="DYNAMITE Predictions")

            if self.config_parameters["plt_PDFs"] == "True":
                plt.plot(R, PR, label="PDF")

            if len(r1) > 0:
                plt.scatter(r1, np.ones(len(r1))*np.amax(n)/10, c="g", edgecolors="w", s=200, linewidth=2, label=("Known planets in system" if len(p1) > 1 else "Known planet in system"), zorder=2)

            if len(r11) > 0:
                plt.scatter(r11, np.ones(len(r11))*np.amax(n)/10, c="purple", marker="s", edgecolors="w", s=200, linewidth=2, label="Highest relative likelihood for inserted planet", zorder=2)

            if len(r12) > 0:
                plt.scatter(r12, np.ones(len(r12))*np.amax(n)/10, c="g", marker="^", edgecolors="k", s=200, linewidth=2, label=("Known non-transiting planets" if len(r12) > 1 else "Known non-transiting planet"), zorder=2)

            if len(r2) > 0:
                plt.scatter(r2, np.ones(len(r2))*np.amax(n)/10, c="r", marker="X", s=200, edgecolors="w", linewidth=2, label=("Known planets removed from system" if len(r2) > 1 else "Known planet removed from system"), zorder=2)

            if len(r3) > 0:
                plt.scatter(r3, np.ones(len(r3))*np.amax(n)/10, c="y", marker='$?$', s=200, label=("Unconfirmed planet candidates" if len(r3) > 1 else "Unconfirmed planet candidate"), zorder=2)

            plt.ylim(0,1.2*np.amax(n))
            plt.xlabel(r"Radius ($R_{\oplus}$)", fontsize=20)
            plt.ylabel("Relative Likelihood", fontsize=20)
            plt.legend(loc=1, fontsize=16)
            ax.tick_params(labelsize=14)
            fig.suptitle((r"$\tau$ Ceti" if system == "tau Ceti" else system) + " Planet Radius Relative Likelihood", fontsize=30)
            plt.savefig(system + "_R_linear.png")

            if self.config_parameters["show_plots"] == "True":
                plt.show()

            elif self.config_parameters["show_plots"] == "False":
                plt.close()

        if self.config_parameters["ind_i"] == "full":
            ik = np.array(ik)
            fname = system + "_inc_full.png"

        elif self.config_paramters["ind_i"] == "truncated":
            ik = np.array([ik[j] if ik[j] < 90 else 180-ik[j] for j in range(len(ik))])     # toggle for allowing inclinations to be greater than 90 degrees or truncated back to 0-90 range
            fname = system + "_inc_trunc.png"

        n, _, _ = plt.hist(ik, bins=np.linspace(0, 180.5, 362), weights=np.ones(len(ik)) / (len(ik)), label="DYNAMITE Predictions")

        if self.config_parameters["plt_PDFs"] == "True":
            plt.plot(inc, Pi, label="PDF")

        if len(i1) > 0:
            plt.scatter(i1, np.ones(len(i1))*np.amax(n)/10, c="g", edgecolors="w", s=200, linewidth=2, label=("Known planets in system" if len(i1) > 1 else "Known planet in system"), zorder=2)

        if len(i11) > 0:
            plt.scatter(i11, np.ones(len(i11))*np.amax(n)/10, c="purple", marker="s", edgecolors="w", s=200, linewidth=2, label="Highest relative likelihood for inserted planet", zorder=2)

        if len(i12) > 0:
            plt.scatter(i12, np.ones(len(i12))*np.amax(n)/10, c="g", marker="^", edgecolors="k", s=200, linewidth=2, label=("Known non-transiting planets" if len(i12) > 1 else "Known non-transiting planet"), zorder=2)

        if len(i2) > 0:
            plt.scatter(i2, np.ones(len(i2))*np.amax(n)/10, c="r", marker="X", s=200, edgecolors="w", linewidth=2, label=("Known planets removed from system" if len(i2) > 1 else "Known planet removed from system"), zorder=2)

        if len(i3) > 0:
            plt.scatter(i3, np.ones(len(i3))*np.amax(n)/10, c="y", marker='$?$', s=200, label=("Unconfirmed planet candidates" if len(i3) > 1 else "Unconfirmed planet candidate"), zorder=2)

        plt.xlabel("Inclination (degrees)")
        plt.ylabel("Relative Likelihood")
        #plt.xlim(np.amin([i1, i11, i12, i2, i3]) - 10, np.amax([i1, i11, i12, i2, i3]) + 10)
        plt.ylim(0, np.amax(n)*1.2)
        plt.legend()
        fig.suptitle(r"$\tau$ Ceti" if system == "tau Ceti" else system + " Inclination Relative Likelihood", fontsize=30)
        plt.savefig(fname)

        if self.config_parameters["show_plots"] == "True":
            plt.show()

        elif self.config_parameters["show_plots"] == "False":
            plt.close()



    def mr_predict(self, meas, value):
        """Calls the mass-radius relationship prediction code"""

        return pfm(measurement=meas, predict=value, dataset="kepler")[0]



if __name__ == '__main__':
    dynamite_plots()
