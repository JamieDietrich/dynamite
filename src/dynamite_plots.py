import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patches as mpatch
from mrexo import predict_from_measurement as pfm

class dynamite_plots:

    def __init__(self, Pk=None, P=None, PP=None, per=None, Rk=None, R=None, PR=None, ik=None, inc=None, Pi=None, deltas=None, ratios=None, targets=None, pers=None, rads=None):
        """Sets up plotting routines"""

        self.config_parameters = {}

        try:
            config_data = np.loadtxt("dynamite_config.txt", dtype=str, delimiter='::')

        except IOError:
            print("Error, configuration file not found!")
            exit()

        for i in range(len(config_data)):
            self.config_parameters[config_data[i, 0]] = config_data[i, 1]

        n_tot = []
        r_tot = []

        if Pk == None or self.config_parameters["saved"] == "True":
            Pk, P, PP, per, Rk, R, PR, ik, inc, Pi, deltas, ratios = self.read_saved_data()

        Pb = list(P)
        Pb.append(730.1)
        n, b, _ = plt.hist(Pk, bins=Pb, weights=np.ones(len(Pk))*10 / len(Pk), label="DYNAMITE Predictions")
        n_tot.append(n)
        plt.close()
        Rb = list(R)
        Rb.append(float(self.config_parameters["radmax"]) + 0.01)
        n, b, _ = plt.hist(Rk, bins=Rb, weights=np.ones(len(Rk))*10 / len(Rk), label="DYNAMITE Predictions")
        r_tot.append(n)
        plt.close()

        if self.config_parameters["plt_P_R"] == "True":
            self.set_up_P_R(targets, pers, rads, n_tot, r_tot)

        if self.config_parameters["plt_tdtp"] == "True":
            self.plot_td_tp()

        if self.config_parameters["plt_deltas"] == "True":
            self.plot_deltas(deltas)

        if self.config_parameters["plt_ratios"] == "True":
            self.plot_ratios(ratios)

        if self.config_parameters["plt_indpars"] == "True":
            self.plot_ind_params(Pk, P, PP, per, Rk, R, PR, ik, inc, Pi)



    def read_saved_data(self, pers, rads):
        """Reads in the saved data"""

        Pk, P, PP, per, Rk, R, PR, ik, inc, Pi, deltas, ratios = np.loadtxt("saved_full_data.txt", delimiter='\t', unpack=True)
        
        return Pk, P, PP, per, Rk, R, PR, ik, inc, Pi, deltas, ratios


    def set_up_P_R(self, targets, pers, rads, n_tot, r_tot):
        """Sets up the P and R distributions for plotting"""

        Pf = np.arange(0.5, 730.01, 0.01)
        Rf = np.arange(float(self.config_parameters["radmin"]), float(self.config_parameters["radmax"]) + 0.001, 0.001)

        if self.config_parameters["saved"] == "False":
            P = np.arange(0.5, 730.1, 0.1)
            R = np.arange(float(self.config_parameters["radmin"]), float(self.config_parameters["radmax"]) + 0.01, 0.01)
            plot_fig_p = []
            plot_fig_r = []

            for i in range(len(n_tot)):
                for j in range(200):
                    plot_fig_p.append(n_tot[i][:1855])
                    plot_fig_r.append(r_tot[i])

            plot_fig_p = np.array(plot_fig_p)
            plot_fig_r = np.array(plot_fig_r)
            np.savetxt("plot_fig_p.txt", plot_fig_p)
            np.savetxt("plot_fig_r.txt", plot_fig_r)

        elif self.config_parameters["saved"] == "True":
            plot_fig_p = np.loadtxt("plot_fig_p.txt")
            plot_fig_r = np.loadtxt("plot_fig_r.txt")
            targets = np.loadtxt("targets.txt", dtype='str', delimiter='\t')

        self.plot_figs(plot_fig_p, plot_fig_r, P, R, Pf, Rf, pers, rads, targets)



    def plot_P_R(self, plot_fig_p, plot_fig_r, P, R, Pf, Rf, pers, rads, targets):
        """Plots P, R, and PR figures"""

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
        ax.set_xticks([0, 90, 990, 9990])
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
                bottom = ["TOI 119", "TOI 286", "TOI 487", "TOI 736", "TOI 797", "TOI 1453", "TOI 1720", "TOI 1726", "TOI 1749"]
                left = ["TOI 256", "TOI 1269"]
                right = ["TOI 174", "TOI 431", "TOI 703", "TOI 1277", "TOI 1730"]
                tl = ["TOI 125"]
                tr = ["TOI 266"]
                bl = ["TOI 836", "TOI 1260"]
                br = ["TOI 1438"]

            elif self.config_parameters["period"] == "syssim":
                bottom = ["TOI 119", "TOI 286", "TOI 487", "TOI 736", "TOI 797", "TOI 1453", "TOI 1720", "TOI 1726", "TOI 1749"]
                left = ["TOI 256", "TOI 1269"]
                right = ["TOI 174", "TOI 431", "TOI 703", "TOI 1277", "TOI 1730"]
                tl = ["TOI 125"]
                tr = ["TOI 266"]
                bl = ["TOI 836", "TOI 1260"]
                br = ["TOI 1438"]

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
            plt.savefig("PRfig_epos_ell.png", bbox_inches='tight')

        elif self.config_parameters["period"] == "syssim":
            plt.xlim(40, 11000)    
            plt.savefig("PRfig_syssim_ell.png", bbox_inches='tight')

        if self.config_parameters["show_plots"] == "True":
            plt.show()

        elif self.config_parameters["show_plots"] == "False":
            plt.close()

        fig, ax = plt.subplots(1,1, figsize=(10,12))
        fig.suptitle(r"Period Relative Likelihoods for $TESS$ Systems", fontsize=30)
        img = ax.imshow(pfpi, cmap=plt.cm.Blues, origin="lower", aspect="auto")

        if self.config_parameters["plt_P_scale"] == "linear":
            xlabels = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
            ax.set_xticks([-50, 1500, 3500, 5500, 7500, 9500, 11500, 13500, 15500, 17500])
            plt.xlim(-500, 18000)

        elif self.config_parameters["plt_P_scale"] == "log":
            ax.set_xscale("Log")
            xlabels = [0.1, 1, 10, 100]
            ax.set_xticks([0, 90, 990, 9990])
            plt.xlim(10, 73000)

        ylabels = [tn for tn in targets_dict.keys() if tn.find("TOI") != -1]
        ax.set_yticks(np.arange(100, len(ylabels)*200 + 100, 200))
        ax.set_xticklabels(xlabels)
        ax.set_yticklabels(ylabels)
        ax.tick_params(labelsize=12)
        cb = fig.colorbar(img)
        cb.set_label(label="Probability normalized to 1 injected planet per system", size=16)
        cb.ax.tick_params(labelsize=12)
        plt.xlabel("Period (days)", fontsize=16)
        plt.ylabel("Systems", fontsize=16)
        pc = 100
        plt.scatter(np.where(np.isclose(Pf, pers[0][0], atol=0.06))[0][0], pc, color="r", s=int(round(rads[0][0]*20, 0)), label="Known Planets")

        for i in range(len(pers)):
            for j in range(len(pers[i])):
                plt.scatter(np.where(np.isclose(Pf, pers[i][j], atol=0.06))[0][0], pc, color="r", s=int(round(rads[i][j]*20, 0)))

            pc += 200

        if self.config_parameters["period"] == "epos":
            plt.savefig("logPfig_epos.png", bbox_inches='tight')

        elif self.config_parameters["period"] == "syssim":
            plt.savefig("logPfig_syssim.png", bbox_inches='tight')

        if self.config_parameters["show_plots"] == "True":
            plt.show()

        elif self.config_parameters["show_plots"] == "False":
            plt.close()
       
        fig, ax = plt.subplots(1,1, figsize=(10,12))
        fig.suptitle(r"Planet Radius Relative Likelihoods for $TESS$ Systems", fontsize=30)
        img = ax.imshow(pfri, cmap=plt.cm.Blues, origin="lower", aspect="auto")
        xlabels = [0, 1, 2, 3, 4, 5]
        ylabels = [tn for tn in targets_dict.keys() if tn.find("TOI") != -1]
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

        plt.savefig("Rfig.png", bbox_inches='tight')

        if self.config_parameters["show_plots"] == "True":
            plt.show()

        elif self.config_parameters["show_plots"] == "False":
            plt.close()



    def plot_deltas(self, deltas):
        """Plots the probability as a function of the separation stability criterion in terms of mutual Hill radii"""

        if self.config_parameters["saved"] == "False":
            deltas = [sum(x) for x in zip(*deltas)]
            deltas = [deltas[i]/46 for i in range(len(deltas))]
            Du = np.arange(0, len(deltas))
            deltas = np.insert(deltas, 8, np.zeros(99))
            Du = np.insert(Du, 8, np.arange(7.01, 8, 0.01))
            
            if self.config_parameters["period"] == "epos":
                np.savetxt("Deltas_epos.txt", np.transpose([Du, deltas]), delimiter="\t")

            elif self.config_parameters["period"] == "syssim":
                np.savetxt("Deltas_syssim.txt", np.transpose([Du, deltas]), delimiter="\t")

        elif self.config_parameters["saved"] == "True":
            if self.config_parameters["period"] == "epos":
                Du, deltas = np.loadtxt("Deltas_epos.txt", delimiter="\t", unpack=True)

            elif self.config_parameters["period"] == "syssim":
                Du, deltas = np.loadtxt("Deltas_syssim.txt", delimiter="\t", unpack=True)

        plt.plot(Du, deltas)

        if self.config_parameters["period"] == "epos":
            plt.savefig("Deltas_epos.png")

        elif self.config_parameters["period"] == "syssim":
            plt.savefig("Deltas_syssim.png")

        if self.config_parameters["show_plots"] == "True":
            plt.show()

        elif self.config_parameters["show_plots"] == "False":
            plt.close()



    def plot_ratios(self, ratios):
        """Plots the period ratios for the current subset of data"""

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
        cdfP = spst.norm(logD, sigma).cdf(Dgrid)
        plt.plot(PRgrid, 17*pdfP, linewidth=4, label=r"$Kepler$ PDF fit")
        plt.xscale("Log")
        plt.xlabel(r"Period Ratio $P_i/P_{i-1}$", fontsize=20)
        plt.ylabel("Occurrence", fontsize=20)
        plt.legend(fontsize=16)
        ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
        ax.set_xticklabels([0.1,1,10])
        ax.tick_params(labelsize=14)
        ax.set_yticks([0, 5, 10, 15, 20])
        plt.suptitle(r"Subsample vs $Kepler$ Period Ratios", fontsize=30)
        plt.xlim(1, 10**1.4)
        plt.savefig(self.config_parameters["mode"] + "_Kepler_period_ratios.png")

        if self.config_parameters["show_plots"] == "True":
            plt.show()

        elif self.config_parameters["show_plots"] == "False":
            plt.close()



    def plot_td_tp(self):
        """Plots the transit probability vs the transit depth for each system in the current data subset"""

        x = []

        for i in range(len(self.tdm)):
            x.append([self.tdm[i], self.tdle[i], self.tdue[i], self.tpm[i], self.tple[i], self.tpue[i]])

        np.savetxt("td_tp.txt", x)
        fig, ax = plt.subplots(figsize=(12,8))
        plt.plot(self.tdm, self.tpm, "o")
        
        names = [set_up(tn)[5] for tn in targets_dict.keys() if tn.find("TOI") != -1]
        #bunch = ["TOI 261", "TOI 266", "TOI 396", "TOI 411", "TOI 487", "TOI 561", "TOI 703", "TOI 797", "TOI 1238", "TOI 1346", "TOI 1453", "TOI 1469", "TOI 1726"]
        bunch = ["TOI 256", "TOI 714", "TOI 736"]
        #bunch = ["TOI 396", "TOI 411", "TOI 286", "TOI 1469", "TOI 487", "TOI 174", "TOI 261", "TOI 1453", "TOI 713", "TOI 1339", "TOI 431", "TOI 282", "TOI 1346", "TOI 1238", "TOI 266", "TOI 1726", "TOI 797", "TOI 1269", "TOI 703", "TOI 1730", "TOI 696", "TOI 836", "TOI 732", "TOI 1449", "TOI 763", "TOI 1260"]
        bottom = []
        left = []
        right = []

        for i in range(len(names)):
            tdle = self.tdm[i] - self.tdle[i]
            tdue = self.tdm[i] + self.tdue[i]
            tple = self.tpm[i] - self.tple[i]
            tpue = self.tpm[i] + self.tpue[i]
            ellipse = mpatch.Ellipse((np.mean([tdle, tdue]), np.mean([tple, tpue])), tdue-tdle, tpue-tple, alpha=0.1)
            ax.add_patch(ellipse)

            if names[i] not in bunch:
                if names[i] in bottom:
                    plt.annotate(names[i][names[i].find(" "):], (self.tdm[i], self.tpm[i]), textcoords="offset points", xytext=(0,-20), ha='center', fontsize=16)

                elif names[i] in left:
                    plt.annotate(names[i][names[i].find(" "):], (self.tdm[i], self.tpm[i]), textcoords="offset points", xytext=(-30,-5), ha='center', fontsize=16)

                elif names[i] in right:
                    plt.annotate(names[i][names[i].find(" "):], (self.tdm[i], self.tpm[i]), textcoords="offset points", xytext=(25,-5), ha='center', fontsize=16)

                else:
                    plt.annotate(names[i][names[i].find(" "):], (self.tdm[i], self.tpm[i]), textcoords="offset points", xytext=(0,10), ha='center', fontsize=16)

            else:
                if names[i] in bottom:
                    plt.annotate(names[i][names[i].find(" "):], (self.tdm[i], self.tpm[i]), textcoords="offset points", xytext=(0,-20), ha='center', fontsize=16)

                elif names[i] in left:
                    plt.annotate(names[i][names[i].find(" "):], (self.tdm[i], self.tpm[i]), textcoords="offset points", xytext=(-30,-5), ha='center', fontsize=16)

                elif names[i] in right:
                    plt.annotate(names[i][names[i].find(" "):], (self.tdm[i], self.tpm[i]), textcoords="offset points", xytext=(30,-5), ha='center', fontsize=16)

                else:
                    plt.annotate(names[i][names[i].find(" "):], (self.tdm[i], self.tpm[i]), textcoords="offset points", xytext=(0,10), ha='center', fontsize=16)

        ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
        plt.xlim(150, max(self.tdm) + max(self.tdue))
        #plt.xlim(150, 3500)
        #plt.xlim(150, 800)
        #plt.ylim(0.75, 0.85)
        #plt.xlim(150, 1150)
        #plt.ylim(0.8, 0.89)
        ax.tick_params(labelsize=14)
        plt.xlabel("Transit Depth (ppm)", fontsize=20)
        plt.ylabel("Transit Probability", fontsize=20)
        plt.show()


    def plot_ind_params(self, Pk, P, PP, per, Rk, R, PR, ik, inc, Pi):
        """Plots histograms of each individual distribution"""

        system = self.config_parameters["system"]
        p1 = []     # known planet period
        p11 = []    # inserted planet period for iteration
        p12 = [20, 49.41, 162.87, 636.13]       # known non-transiting planet periods
        p2 = []     # removed known planet period
        p3 = [13.965, 35.362, 94.11]        # unconfirmed planet candidate periods
        r1 = []     # " for planet radius
        r11 = []    # " for planet radius
        r12 = [pfm(measurement=1.75/math.sin(35*math.pi/180), predict="Radius", dataset="kepler")[0], pfm(measurement=1.83/math.sin(35*math.pi/180), predict="Radius", dataset="kepler")[0], pfm(measurement=3.93/math.sin(35*math.pi/180), predict="Radius", dataset="kepler")[0], pfm(measurement=3.93/math.sin(35*math.pi/180), predict="Radius", dataset="kepler")[0]]      # " for planet radius
        r2 = []     # " for planet radius
        r3 = [pfm(measurement=2/math.sin(35*math.pi/180), predict="Radius", dataset="kepler")[0], pfm(measurement=3.1/math.sin(35*math.pi/180), predict="Radius", dataset="kepler")[0], pfm(measurement=3.6/math.sin(35*math.pi/180), predict="Radius", dataset="kepler")[0]]     # " for planet radius
        i1 = []
        i11 = []
        i12 = [35, 35, 35, 35]
        i2 = []
        i3 = [35, 35, 35]
        ul = 0.01       # upper limit of plot

        if self.config_parameters["ind_P"] == "linear_zoom":
            n, b, _ = plt.hist(Pk, bins=np.arange(0.5, round(5*max(per)) + 1, 0.5), weights=np.ones(len(Pk))*2 / (len(Pk)), label="DYNAMITE Predictions")
            plt.plot(P, PP, label="PDF")
            plt.xlabel("Period (days)")
            plt.ylabel("Relative Likelihood")
            plt.xlim(0, 5*max(per))
            plt.legend()
            fig.suptitle(system + " Period Relative Likelihood", fontsize=30)
            plt.savefig(system + "_P_linear_zoom.png")

            if self.config_parameters["show_plots"] == "True":
                plt.show()

            elif self.config_parameters["show_plots"] == "False":
                plt.close()

        elif self.config_parameters["ind_P"] == "linear":
            Pb = list(P)
            Pb.append(730.1)
            n, b, _ = plt.hist(Pk, bins=Pb, weights=np.ones(len(Pk))*10 / len(Pk), label="DYNAMITE Predictions")
            plt.plot(P, PP, label="PDF")
            plt.xlabel("Period (days)")
            plt.ylabel("Relative Likelihood")
            plt.xlim(0, 5*max(per))
            plt.legend()
            fig.suptitle(system + " Period Relative Likelihood", fontsize=30)
            plt.savefig(system + "_P_linear.png")

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
            hist = np.histogram(Pk, bins=bins)
            hist_norm = hist[0]/widths
            fig, ax = plt.subplots(figsize=(12,8))
            plt.bar(bins[:-1], hist_norm/100000, widths, label="DYNAMITE Predictions")
            PPl = np.interp(Pl, P, PP)
            plt.plot(Pl, PPl, color='#ff7f0e', label="PDF", linewidth=4)
            
            if len(p1) > 0:
                plt.scatter(p1, np.ones(len(p1))*ul/20, c="g", s=r1*100, label=("Known planets in system" if len(p1) > 1 else "Known planet in system"), zorder=2)

            if len(p11) > 0:
                plt.scatter(p11, np.ones(len(p11))*ul/20, c="purple", marker="s", s=r11*100, label="Highest relative likelihood for inserted planet", zorder=2)

            if len(p12) > 0:
                plt.scatter(p12, np.ones(len(p12))*ul/20, c="g", marker="^", s=r12*100, label=("Known non-transiting planets" if len(p12) > 1 else "Known non-transiting planet"), zorder=2)

            if len(p2) > 0:
                plt.scatter(p2, np.ones(len(p2))*ul/20, c="w", marker="X", edgecolors="k", s=r2*100, linewidth=2, label=("Known planets removed from system" if len(p2) > 1 else "Known planet removed from system"), zorder=2)

            if len(p3) > 0:
                plt.scatter(p3, np.ones(len(p3))*ul/20, c="y", marker='$?$', s=r3*100, label=("Unconfirmed planet candidates" if len(p3) > 1 else "Unconfirmed planet candidate"), zorder=2)

            plt.xlabel("Log Period (days)", fontsize=20)
            plt.ylabel("Relative Likelihood", fontsize=20)
            plt.xscale("Log")
            plt.xlim(0.5, 730)
            plt.ylim(0, ul)
            plt.legend(fontsize=16)
            ax.tick_params(labelsize=14)
            ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
            fig.suptitle(system + " Period Relative Likelihood", fontsize=30)
            plt.savefig(system + "_P_log.png")

            if self.config_parameters["show_plots"] == "True":
                plt.show()

            elif self.config_parameters["show_plots"] == "False":
                plt.close()

        if self.config_parameters["ind_R"] == "linear_zoom":
            Rk = np.array(Rk)
            n, b, _ = plt.hist(Rk, bins=np.arange(R[0], R[-1] + 0.01, 0.005), weights=np.ones(len(Rk))*20 / len(Rk), label="DYNAMITE Predictions")
            plt.plot(R, PR, label="PDF")
            plt.xlabel(r"Radius ($R_{\oplus}$)")
            plt.ylabel("Relative Likelihood")
            plt.legend()
            plt.savefig(system + "_R_linear_zoom.png")

            if self.config_parameters["show_plots"] == "True":
                plt.show()

            elif self.config_parameters["show_plots"] == "False":
                plt.close()

        elif self.config_parameters["ind_R"] == "linear":
            Rb = list(R)
            Rb.append(float(self.config_parameters["radmax"]) + 0.01)
            n, b, _ = plt.hist(Rk, bins=Rb, weights=np.ones(len(Rk))*10 / len(Rk), label="DYNAMITE Predictions")
            plt.plot(R, PR, label="PDF")

            if len(r1) > 0:
                plt.scatter(r1, np.ones(len(r1)*ul/20), c="g", label=("Known planets in system" if len(p1) > 1 else "Known planet in system"), zorder=2)

            if len(r11) > 0:
                plt.scatter(r11, np.ones(len(r11))*ul/20, c="purple", marker="s", label="Highest relative likelihood for inserted planet", zorder=2)

            if len(r12) > 0:
                plt.scatter(r12, np.ones(len(r12)*ul/20, c="g", marker="^", label=("Known non-transiting planets" if len(r12) > 1 else "Known non-transiting planet"), zorder=2)

            if len(r2) > 0:
                plt.scatter(r2, np.ones(len(r2))*ul/20, c="w", marker="X", edgecolors="k", linewidth=2, label=("Known planets removed from system" if len(r2) > 1 else "Known planet removed from system"), zorder=2)

            if len(r3) > 0
                plt.scatter(r3, np.ones(len(r3)*ul/20, c="y", marker='$?$', label=("Unconfirmed planet candidates" if len(r3) > 1 else "Unconfirmed planet candidate"), zorder=2)

            plt.xlabel(r"Radius ($R_{\oplus}$)")
            plt.ylabel("Relative Likelihood")
            plt.legend()
            fig.suptitle(system + " Planet Radius Relative Likelihood", fontsize=30)
            plt.savefig(system + "_R_linear.png")

            if self.config_parameters["show_plots"] == "True":
                plt.show()

            elif self.config_parameters["show_plots"] == "False":
                plt.close()

        if self.config_parameters["ind_i"] == "full":
            ik = np.array(ik)

        elif self.config_paramters["ind_i"] == "truncated":
            ik = np.array([ik[j] if ik[j] < 90 else 180-ik[j] for j in range(len(ik))])     # toggle for allowing inclinations to be greater than 90 degrees or truncated back to 0-90 range

        plt.hist(ik, bins=np.linspace(0, 180.5, 362), weights=np.ones(len(ik)) / (len(ik)), label="DYNAMITE Predictions")
        plt.plot(inc, Pi, label="PDF")
        plt.xlabel("Inclination (degrees)")
        plt.ylabel("Relative Likelihood")
        plt.legend()
        fig.suptitle(system + " Inclination Relative Likelihood", fontsize=30)

        if self.config_parameters["show_plots"] == "True":
            plt.show()

        elif self.config_parameters["show_plots"] == "False":
            plt.close()


if __name__ == '__main__':
    dynamite_plots()
