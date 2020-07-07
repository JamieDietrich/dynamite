import math
import itertools
import numpy as np
import dynamite_plots
import scipy.stats as spst
import scipy.optimize as spopt
import astropy.constants as const
import scipy.interpolate as spinter
from dynamite_plots import dynamite_plots
from dynamite_targets import dynamite_targets
from mrexo import predict_from_measurement as pfm

class dynamite:

    def __init__(self):
        """Runs the script"""

        np.random.seed(1)
        self.config_parameters = {}

        try:
            config_data = np.loadtxt("dynamite_config.txt", dtype=str, delimiter='::')

        except IOError:
            print("Error, configuration file not found!")
            exit()

        for i in range(len(config_data)):
            self.config_parameters[config_data[i, 0]] = config_data[i, 1]

        targets_dict = dynamite_targets().get_targets(self.config_parameters["mode"], self.config_parameters["system"], self.config_parameters["radmax"])

        if len(targets_dict) == 0:
            print("Error: No targets selected!")
            exit()
       
        pers = []
        rads = []
        self.n_tot_log = []
        self.PRs = []
        self.PD = []
        self.DD = 1
        self.tdm = []
        self.tdue = []
        self.tdle = []
        self.tpm = []
        self.tpue = []
        self.tple = []
        self.targets = []

        def set_up(target):
            """Sets up target"""

            def get_arccos(star_pars, planet_pars):
                return round(np.arccos(planet_pars[0]/(self.K3(planet_pars[1], star_pars[2])/(star_pars[0]*const.R_sun.cgs.value)))*180/math.pi, 3)
               
            t = list(targets_dict[target])
            for x in range(len(t)):
                for y in range(len(t[x])):
                    if isinstance(t[x][y], tuple):
                        t[x][y] = locals()[t[x][y][0]](t[0],t[x][y][1])

            return t[0][0], t[0][1], t[0][2], t[0][3], np.array(t[1:]), target


        if self.config_parameters["saved"] == "False":
            if self.config_parameters["mode"] == "all":
                for tn in targets_dict.keys():
                    R_star, Rse, M_star, Mse, target, target_name = set_up(self.config_parameters["system"])
                    target = target[target[:, 2].argsort()]
                    self.run_monte_carlo(R_star, Rse, M_star, Mse, target, target_name)

            elif self.config_parameters["mode"] == "TESS":
                for tn in targets_dict.keys():
                    if tn.find("TOI") != -1:
                        R_star, Rse, M_star, Mse, target, target_name = set_up(tn)
                        target = target[target[:, 2].argsort()]
                        self.run_monte_carlo(R_star, Rse, M_star, Mse, target, target_name)

            elif self.config_parameters["mode"] == "Kepler":
                for tn in targets_dict.keys():
                    if tn.find("Kepler") != -1 or tn.find("K2") != -1 or tn.find("KOI") != -1:
                        R_star, Rse, M_star, Mse, target, target_name = set_up(tn)
                        target = target[target[:, 2].argsort()]
                        self.run_monte_carlo(R_star, Rse, M_star, Mse, target, target_name)

            elif self.config_parameters["mode"] == "single":
                R_star, Rse, M_star, Mse, target, target_name = set_up(self.config_parameters["system"])
                target = target[target[:, 2].argsort()]
                self.run_monte_carlo(R_star, Rse, M_star, Mse, target, target_name)



    def run_monte_carlo(self, R_star, Rse, M_star, Mse, target, target_name):
        """Runs the Monte Carlo analysis."""

        inc = [target[i][0] for i in range(len(target))]
        rad = [target[i][1] for i in range(len(target))]
        per = [target[i][2] for i in range(len(target))]

        for i in range(1, len(per)):
            self.PRs.append(per[i]/per[i-1])

        p0 = min(per)
        r1 = min(rad)
        r2 = max(rad)
        P = np.arange(0.5, 730.1, 0.1)

        if self.config_parameters["period"] == "epos":
            PP, PD = self.epos_pers(p0, per, rad, P, M_star)

        elif self.config_parameters["period"] == "syssim":
            PP, PD = self.syssim_pers(per, rad, P, M_star)

        self.PD.append(PD)

        if self.config_parameters["radius"] == "epos":
            R, PR, cdfR = self.epos_rads(r1, r2)

        elif self.config_parameters["radius"] == "syssim":
            R, PR, cdfR = self.syssim_rads(self.config_parameters["radtype"], rad)

        i = np.linspace(0, 180.1, 1802)
        fi = np.zeros(len(i))
        rylgh = 2
        ibs = []
        fib = []
        incn = []

        for case in [[False] + list(t) for t in list(itertools.product([False,True], repeat=len(inc)-1))]:
            incn.append([180-inc[i] if case[i] else inc[i] for i in range(0, len(inc))])

        for j in range(len(i)):
            for k in range(len(incn)):
                test = 0

                for m in range(len(incn[k])):
                    test += spst.rayleigh.pdf(abs(incn[k][m]-i[j]), rylgh)

                ibs.append(i[j])
                fib.append(test)

        ib = ibs[np.where(fib == max(fib))[0][0]]

        if ib > 90:
            ib = 180 - ib

        inew = np.linspace(0, 10, 101)
        finew = spst.rayleigh.pdf(inew + ib, ib, rylgh)

        if len(inc) == 2:
            finew = finew*0.62

            for j in range(len(i)):
                fi[j] = np.sin(i[j]*math.pi/180)*76/300

            fi = fi*76/(300*np.trapz(fi, i))


        elif len(inc) == 3:
            finew = finew*0.81

            for j in range(len(i)):
                fi[j] = np.sin(i[j]*math.pi/180)*38/300

            fi = fi*38/(300*np.trapz(fi, i))

        i_ib = np.where(np.isclose(i,ib))[0][0]

        for j in range(len(inew)):
            fi[i_ib + j] += finew[j]

        cdfi = np.array([1 - math.exp(-(inew[j])**2/(2*(rylgh)**2)) for j in range(len(inew))])
        Pi = fi/2
        Pk = []
        Rk = []
        ik = []
        Nk = []

        for k in range(10000):
            N = 0

            for j in range(len(PP)):
                if np.random.rand() < PP[j]:
                    N += 1
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

            Nk.append(N)

        Pi = np.arange(0.5, 146.001, 0.001)

        if self.config_parameters["period"] == "epos":
            PPi, _ = self.epos_pers(p0, per, rad, Pi, M_star)

            if self.config_parameters["mode"] == "TESS":
                low_gap_P = ["TOI 561", "TOI 431", "TOI 1238", "TOI 732", "TOI 696", "TOI 175", "TOI 663", "TOI 1469", "TOI 1260", "TOI 270", "TOI 396", "TOI 836", "TOI 411", "TOI 1269", "TOI 1453", "TOI 714", "TOI 1749", "TOI 125", "TOI 1438", "TOI 119", "TOI 763", "TOI 1136", "TOI 1064", "TOI 266", "TOI 178", "TOI 776", "TOI 1339", "TOI 214", "TOI 700", "TOI 1266", "TOI 553", "TOI 699", "TOI 1277"]
                interior = ["TOI 282"]
        
                if target_name in low_gap_P or target_name in interior:
                    PPz = PPi

                else:
                    Pz = np.where(PPi == 0)[0]
                    PPz = PPi[Pz[1]:Pz[-1]]

            else:
                PPz = PPi

            Pm = Pi[np.where(PPi == np.amax(PPz))[0][0]]
            PPm = np.amax(PPz)

        elif self.config_parameters["period"] == "syssim":
            PPi, _ = self.syssim_pers(per, rad, Pi, M_star)
            Pm = Pi[np.where(PPi == np.amax(PPi))[0][0]]
            PPm = np.amax(PPi)

        Ple = Pm - Pi[np.where((Pi < Pm) & (PPi < 0.606*PPm))][-1]
        Pue = Pi[np.where((Pi > Pm) & (PPi < 0.606*PPm))][0] - Pm
        Rm = np.percentile(Rk, 50)
        Rle = Rm - np.percentile(Rk, 16)
        Rue = np.percentile(Rk, 84) - Rm
        tdm = (Rm*const.R_earth.cgs.value/(R_star*const.R_sun.cgs.value))**2*1e6
        tdle = 2*(Rm*const.R_earth.cgs.value/(R_star*const.R_sun.cgs.value))**2*math.sqrt((Rle/Rm)**2 + (Rse/R_star)**2)*1e6
        tdue = 2*(Rm*const.R_earth.cgs.value/(R_star*const.R_sun.cgs.value))**2*math.sqrt((Rue/Rm)**2 + (Rse/R_star)**2)*1e6
        ntrans = 0
        ntl = 0
        ntu = 0

        for j in range(len(ik)):
            if math.cos(ik[j]*math.pi/180) < (R_star*const.R_sun.cgs.value + Rm*const.R_earth.cgs.value)/self.K3(Pm, M_star):
                ntrans += 1

            if math.cos(ik[j]*math.pi/180) < ((R_star - Rse)*const.R_sun.cgs.value + (Rm - Rle)*const.R_earth.cgs.value)/self.K3(Pm, (M_star + Mse)):
                ntl += 1

            if math.cos(ik[j]*math.pi/180) < ((R_star + Rse)*const.R_sun.cgs.value + (Rm + Rue)*const.R_earth.cgs.value)/self.K3(Pm, (M_star - Mse)):
                ntu += 1

        tpm = ntrans/len(ik)
        tple = max(1e-3, ntrans/len(ik) - ntl/len(ik))
        tpue = max(1e-3, ntu/len(ik) - ntrans/len(ik))
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

        print("\t\t" + target_name + " & $" + str(Pm) + "^{" + str(Pue) + "}_{" + str(Ple) + "}$ & $" + str(Rm) + "^{" + str(Rue) + "}_{" + str(Rle) + "}$ & $" + str(R_star) + "\pm" + str(Rse) + "$ & $" + str(tdm) + "^{" + str(tdue) + "}_{" + str(tdle) + "}$ & $" + str(tpm) + "^{" + str(tpue) + "}_{" + str(tple) + "}$ \\\\")

        self.tdm.append(tdm)
        self.tdle.append(tdle)
        self.tdue.append(tdue)
        self.tpm.append(tpm)
        self.tpue.append(tpue)
        self.tple.append(tple)
        self.targets.append([target_name, Pm, Ple, Pue, Rm, Rle, Rue])

        if self.config_parameters["plot"] == "True":
            plots = dynamite_plots(Pk, P, PP, per, Rk, R, PR, ik, i, Pi)



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
            m[k] = pfm(measurement=rad[k], predict='Mass', dataset='kepler')[0]

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

        Du = np.arange(0, max(fD) + self.DD, self.DD)
        fDu = np.zeros(len(Du))

        for i in range(len(fD)):
            j = int(fD[i] // self.DD)
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
            m[k] = pfm(measurement=rad[k], predict='Mass', dataset='kepler')[0]

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

        Du = np.arange(0, max(fD) + self.DD, self.DD)
        fDu = np.zeros(len(Du))

        for i in range(len(fD)):
            j = int(fD[i] // self.DD)
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



    def K3(self, P, M):
        """Calculates semi-major axis in cm using period in days and mass in solar masses"""

        seconds_per_day = 86400
        
        return (const.G.cgs.value*(M*const.M_sun.cgs.value)*(P*seconds_per_day)**2/(4*math.pi**2))**(1/3)

dynamite()