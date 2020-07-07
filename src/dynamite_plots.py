import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patches as mpatch
from mrexo import predict_from_measurement as pfm

class dynamite_plots:

    def __init__(self, Pk=None, P=None, PP=None, per=None, Rk=None, R=None, PR=None, ik=None, i=None, Pi=None):
        """Sets up plotting routines"""

        self.config_parameters = {}

        try:
            config_data = np.loadtxt("dynamite_config.txt", dtype=str, delimiter='::')

        except IOError:
            print("Error, configuration file not found!")
            exit()

        for i in range(len(config_data)):
            self.config_parameters[config_data[i, 0]] = config_data[i, 1]

        self.n_tot = []
        self.r_tot = []

        if Pk == None:
            self.read_saved_data()

        if self.config_parameters["plt_P_R"] == "True":
            self.set_up_P_R()

        if self.config_parameters["plt_tdtp"] == "True":
            self.plot_td_tp()

        if self.config_parameters["plt_deltas"] == "True":
            self.plot_deltas()

        if self.config_parameters["plt_ratios"] == "True":
            self.plot_ratios()

        if self.config_parameters["plt_indpars"] == "True":
            self.plot_ind_params(Pk, P, PP, per, Rk, R, PR, ik, i, Pi)



    def read_saved_data(self):
        """Reads in the saved data"""

        return None



    def set_up_P_R(self):
        """Sets up the P and R distributions for plotting"""

        P = np.arange(0.5, 730.1, 0.1)
        Pf = np.arange(0.1, 730.01, 0.01)
        R = np.arange(0.1, float(self.config_parameters["radcut"]) + 0.01, 0.01)
        Rf = np.arange(0.1, float(self.config_parameters["radcut"]) + 0.001, 0.001)
        plot_fig_p = []
        plot_fig_r = []

        for i in range(len(self.n_tot)):
            for j in range(200):
                plot_fig_p.append(self.n_tot[i][:1855])
                plot_fig_r.append(self.r_tot[i])

        plot_fig_p = np.array(plot_fig_p)
        plot_fig_r = np.array(plot_fig_r)
        np.savetxt("plot_fig_p.txt", plot_fig_p)
        np.savetxt("plot_fig_r.txt", plot_fig_r)
        np.savetxt("targets.txt", self.targets, fmt='%s', delimiter='\t')
        self.plot_figs(plot_fig_p, plot_fig_r, pers, rads, P, R, Pf, Rf, self.targets)



    def plot_P_R(self, plot_fig_p, plot_fig_r, pers, rads, P, R, Pf, Rf, targets):
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

        bottom = ["TOI 119", "TOI 286", "TOI 487", "TOI 736", "TOI 797", "TOI 1453", "TOI 1720", "TOI 1726", "TOI 1749"]        # text placements - edit each time for optimal annotating
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
        #plt.xlim(150, 11000)
        plt.xlim(40, 11000)
        plt.ylim(750, 4200)
        #plt.savefig("PRfig_epos_ell.png", bbox_inches='tight')
        plt.savefig("PRfig_syssim_ell.png", bbox_inches='tight')
        #plt.show()
        plt.close()

        fig, ax = plt.subplots(1,1, figsize=(10,12))
        fig.suptitle(r"Period Relative Likelihoods for $TESS$ Systems", fontsize=30)
        img = ax.imshow(pfpi, cmap=plt.cm.Blues, origin="lower", aspect="auto")
        ax.set_xscale("Log")
        #xlabels = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
        xlabels = [0.1, 1, 10, 100]
        ylabels = [tn for tn in targets_dict.keys() if tn.find("TOI") != -1]
        #ax.set_xticks([-50, 1500, 3500, 5500, 7500, 9500, 11500, 13500, 15500, 17500]) 
        ax.set_xticks([0, 90, 990, 9990])
        ax.set_yticks(np.arange(100, len(ylabels)*200 + 100, 200))
        ax.set_xticklabels(xlabels)
        ax.set_yticklabels(ylabels)
        ax.tick_params(labelsize=12)
        cb = fig.colorbar(img)
        cb.set_label(label="Probability normalized to 1 injected planet per system", size=16)
        cb.ax.tick_params(labelsize=12)
        plt.xlabel("Period (days)", fontsize=16)
        plt.ylabel("Systems", fontsize=16)
        plt.xlim(10, 73000)
        #plt.xlim(-500, 18000)
        pc = 100
        plt.scatter(np.where(np.isclose(Pf, pers[0][0], atol=0.06))[0][0], pc, color="r", s=int(round(rads[0][0]*20, 0)), label="Known Planets")

        for i in range(len(pers)):
            for j in range(len(pers[i])):
                plt.scatter(np.where(np.isclose(Pf, pers[i][j], atol=0.06))[0][0], pc, color="r", s=int(round(rads[i][j]*20, 0)))

            pc += 200

        #plt.savefig("logPfig_epos.png", bbox_inches='tight')
        plt.savefig("logPfig_syssim.png", bbox_inches='tight')
        #plt.show()
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
        #plt.show()
        plt.close()



    def plot_deltas(self):
        """Plots the probability as a function of the separation stability criterion in terms of mutual Hill radii"""

        PDs = [sum(x) for x in zip(*self.PD)]
        PDs = [PDs[i]/46 for i in range(len(PDs))]
        Du = np.arange(0, len(PDs), self.DD)
        PDs = np.insert(PDs, 8, np.zeros(99))
        Du = np.insert(Du, 8, np.arange(7.01, 8, 0.01))
        
        f = open("Deltas_syssim.txt", "w")

        for i in range(len(Du)):
            f.write(str(Du[i]) + "\t" + str(PDs[i]) + "\n")

        f.close()
        plt.plot(Du, PDs)
        plt.show()



    def plot_ratios():
        """Plots the period ratios for the current subset of data"""

        self.PRs = np.array(self.PRs)
        fig, ax = plt.subplots(figsize=(12,8))
        plt.hist(self.PRs, bins=np.logspace(0,1.4,15), label="TESS Systems")

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
        plt.suptitle(r"$Kepler$ and $TESS$ Period Ratios", fontsize=30)
        plt.xlim(1, 10**1.4)
        plt.show()


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


    def plot_ind_params(self, Pk, P, PP, per, Rk, R, PR, ik, i, Pi):
        """Plots histograms of each individual distribution"""

        p1 = []     # known planet period
        p11 = []    # inserted planet period for iteration
        p12 = [20, 49.41, 162.87, 636.13]       # known non-transiting planet periods
        p2 = []     # removed known planet period
        p3 = [13.965, 35.362, 94.11]        # unconfirmed planet candidate periods
        r1 = []     # " for planet radius
        r11 = []    # " for planet radius
        r12 = [round(pfm(measurement=1.75/math.sin(35*math.pi/180), predict="Radius", dataset="kepler")[0]*100, 0), round(pfm(measurement=1.83/math.sin(35*math.pi/180), predict="Radius", dataset="kepler")[0]*100, 0), round(pfm(measurement=3.93/math.sin(35*math.pi/180), predict="Radius", dataset="kepler")[0]*100, 0), round(pfm(measurement=3.93/math.sin(35*math.pi/180), predict="Radius", dataset="kepler")[0]*100, 0)]      # " for planet radius
        r2 = []     # " for planet radius
        r3 = [round(pfm(measurement=2/math.sin(35*math.pi/180), predict="Radius", dataset="kepler")[0]*100, 0), round(pfm(measurement=3.1/math.sin(35*math.pi/180), predict="Radius", dataset="kepler")[0]*100, 0), round(pfm(measurement=3.6/math.sin(35*math.pi/180), predict="Radius", dataset="kepler")[0]*100, 0)]     # " for planet radius
        ul = 0.01

        n, b, _ = plt.hist(Pk, bins=np.arange(0.5, round(5*max(per)) + 1, 0.5), weights=np.ones(len(Pk))*2 / (len(Pk)), label="DYNAMITE Predictions")
        plt.plot(P, PP, label="PDF")
        plt.xlabel("Period (days)")
        plt.ylabel("Relative Likelihood")
        plt.xlim(0, 5*max(per))
        plt.legend()
        #plt.show()
        Pb = list(P)
        Pb.append(730.1)
        n, b, _ = plt.hist(Pk, bins=Pb, weights=np.ones(len(Pk))*10 / len(Pk), label="DYNAMITE Predictions")
        self.n_tot.append(n)
        plt.plot(P, PP, label="PDF")
        plt.xlabel("Period (days)")
        plt.ylabel("Relative Likelihood")
        plt.xlim(0, 5*max(per))
        plt.legend()
        #plt.show()
        plt.close()
        Pl = np.logspace(-0.3, 2.86, 317)
        Plb = list(Pl)
        Plb.append(10**2.895)
        #nl, _, _ = plt.hist(Pk, bins=Plb, label="DYNAMITE Predictions")
        bins = np.array(Plb)
        widths = (bins[1:] - bins[:-1])
        hist = np.histogram(Pk, bins=bins)
        hist_norm = hist[0]/widths
        fig, ax = plt.subplots(figsize=(12,8))
        plt.bar(bins[:-1], hist_norm/100000, widths, label="DYNAMITE Predictions")
        PPl = np.interp(Pl, P, PP)
        plt.plot(Pl, PPl, color='#ff7f0e', label="PDF", linewidth=4)
        #plt.scatter(p1, np.ones(len(p1))*ul/20, c="g", s=r1, label="Known planets still in system", zorder=2)
        #plt.scatter(p11, np.ones(len(p11))*ul/20, c="purple", marker="s", s=r11, label="Highest relative likelihood for inserted planet", zorder=2)
        plt.scatter(p12, np.ones(len(p12))*ul/20, c="g", marker="^", s=r12, label=("Known non-transiting planets") if len(p12) > 1 else ("Known non-transiting planet"), zorder=2)
        #plt.scatter(p2, np.ones(len(p2))*ul/20, c="w", marker="X", edgecolors="k", s=r2, linewidth=2, label=("Known planets removed from system") if len(p2) > 1 else ("Known planet removed from system"), zorder=2)
        plt.scatter(p3, np.ones(len(p3))*ul/20, c="y", marker='$?$', s=r3, label=("Unconfirmed planet candidates") if len(p3) > 1 else ("Unconfirmed planet candidate"), zorder=2)
        #self.n_tot_log.append(nl)
        plt.xlabel("Log Period (days)", fontsize=20)
        plt.ylabel("Relative Likelihood", fontsize=20)
        plt.xscale("Log")
        plt.xlim(0.5, 730)
        plt.ylim(0, ul)
        plt.legend(fontsize=16)
        ax.tick_params(labelsize=14)
        ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
        fig.suptitle("tau Ceti", fontsize=30)
        plt.show()
        exit()
        plt.close()
        Rk = np.array(Rk)
        n, b, _ = plt.hist(Rk, bins=np.arange(R[0], R[-1] + 0.01, 0.005), weights=np.ones(len(Rk))*20 / len(Rk), label="DYNAMITE Predictions")
        plt.plot(R, PR, label="PDF")
        plt.scatter([1.78, 4.12, 2.26, 2.95], np.ones(4)*0.1, c="g", label="Known planets still in system", zorder=2)
        plt.scatter(1.5, 0.1, c="r", marker="x", label="Known planet removed from system", zorder=2)
        plt.xlabel(r"Radius ($R_{\oplus}$)")
        plt.ylabel("Relative Likelihood")
        plt.legend()
        plt.show()
        Rb = list(R)
        Rb.append(15.01)
        n, b, _ = plt.hist(Rk, bins=Rb, weights=np.ones(len(Rk))*10 / len(Rk), label="DYNAMITE Predictions")
        self.r_tot.append(n)
        plt.plot(R, PR, label="PDF")
        plt.xlabel(r"Radius ($R_{\oplus}$)")
        plt.ylabel("Relative Likelihood")
        plt.legend()
        #plt.show()
        plt.close()
        ik = np.array(ik)
        ik = np.array([ik[i] if ik[i] < 90 else 180-ik[i] for i in range(len(ik))])     # toggle for allowing inclinations to be greater than 90 degrees or truncated back to 0-90 range
        plt.hist(ik, bins=np.linspace(0, 180.5, 362), weights=np.ones(len(ik)) / (len(ik)), label="DYNAMITE Predictions")
        plt.plot(i, Pi, label="PDF")
        plt.xlabel("Inclination (degrees)")
        plt.ylabel("Relative Likelihood")
        plt.legend()
        #plt.show()
        plt.close()



    def plot_R_I(self):
        """Plots 2D histogram of R and i values along with planet data if known"""

        fig, ax = plt.subplots()
        n, xedges, yedges, _ = ax.hist2d(Rk, ik, bins=[50,27], cmap=plt.cm.Blues)
        plt.plot(1.5, 88.6, "r*", label="Planet Data")
        #plt.plot(np.mean(Rk), np.mean(ik), "bo", label="DYNAMITE Predictions")
        #plt.plot(np.percentile(Rk, 50), np.percentile(ik, 50), "bo", label="DYNAMITE Predictions")
        plt.plot(xedges[16]-(xedges[16]-xedges[15])/2, yedges[26]-(yedges[26]-yedges[25])/2, "bo", label="DYNAMITE Predictions")
        ellipse = mpatch.Ellipse((1.565, 88.6), 0.55, 0.6, alpha=0.3, color="r")
        ax.add_patch(ellipse)
        ellipse = mpatch.Ellipse((1.59, 88.6), 1.1, 1.2, alpha=0.2, color="r")
        ax.add_patch(ellipse)
        ellipse = mpatch.Ellipse((1.635, 88.6), 1.65, 1.8, alpha=0.1, color="r")
        ax.add_patch(ellipse)
        plt.xlim(0.5,5)
        plt.ylim(84,90)
        plt.legend(loc=4, fontsize=16)
        plt.xlabel(r"Planet Radius ($R_\oplus$)", fontsize=18)
        plt.ylabel("Inclination (degrees)", fontsize=18)
        ax.tick_params(labelsize=12)
        plt.title("Predicted and Observed Radius and Inclination for Kepler-154 f", fontsize=20)
        plt.show()

if __name__ == '__main__':
    dynamite_plots()
