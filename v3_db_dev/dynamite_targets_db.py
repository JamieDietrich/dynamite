### DYNAmical Multi-planet Injection TEster (DYNAMITE) ###
### Targets List w/ SQLite DB ###
### Jeremy Dietrich ###
### jdietrich@asu.edu ###
### 2023 August 22 ###
### Version 3.0 ###
### Dietrich & Apai (2020), AJ, 160, 107D ###
### https://iopscience.iop.org/article/10.3847/1538-3881/aba61d ###
### Dietrich & Apai (2021), AJ, 161, 17D ###
### Dietrich, Apai, & Malhotra (2022), AJ, 163, 88D ###
### Basant, Dietrich, & Apai (2022), AJ, 164, 12B ###
### Basant, Dietrich, & Apai (2022), RNAAS, 6, 213 ###

import sqlite3

class dynamite_targets_db:

    def __init__(self, dbname):
        """Creates the database connection for the target dictionary"""

        self.DBNAME = "dynamite_targets.db" if dbname == None else dbname



    def get_targets(self, mode, system, radmax, removed=[]):
        """Gets the specific set of targets based on config parameters from the database"""

        def create_tuple(val, upp, low):
            """Tests parameters for value vs. None and creates tuples"""

            if isinstance(low, str):
                print(val, upp, low)

            if upp != None and upp != 0 and val != None and val/upp < 10:
                upp = val/10

            if low != None and low != 0 and val != None and val/low > -10:
                low = val/-10

            if val != None and upp != None and low != None:
                return (val, upp, low)

            elif val != None and ((upp != None and low == None) or (upp == None and low != None) or (upp == None and low == None)):
                return (val, "?", "?")

            elif val == None:
                if upp != None and low == None:
                    return ("?", upp, "?")

                elif upp == None and low != None:
                    return ("?", "?", low)

                return ("?", "?", "?")

        try:
            dbconn = sqlite3.connect(self.DBNAME)
            dbcursor = dbconn.cursor()

            targets = {}
            dbcursor.execute("select p1.tname, min(period) as minperiod from planet p1 where not exists (select * from planet p where p.tname = p1.tname and p.pradius > " + str(radmax) + ") group by p1.tname order by minperiod")
            target_names = [i[0] for i in dbcursor.fetchall()]

            for n in target_names:
                dbcursor.execute("select sradius, sradius_unc_upper, smass, smass_unc_upper, stemp, rv_lim from target where tname = '" + n + "'")
                result = dbcursor.fetchone()
                target_params = [[i for i in result]]
                #dbcursor.execute("select period, pradius, pmass, pinc, pecc, pname from planet where tname = '" + n + "' order by period")
                dbcursor.execute("select * from planet where tname = '" + n + "' order by period")
                planets = dbcursor.fetchall()

                for p in planets:
                    #period, radius, mass, inc, ecc, name = p
                    _, name, period, pul, pll, radius, rul, rll, mass, mul, mll, inc, iul, ill, ecc, eul, ell = p

                    if isinstance(pll, str) or isinstance(rll, str) or isinstance(mll, str) or isinstance(ill, str) or isinstance(ell, str):
                        print(p)
                    #target_params.append([period, radius if radius != None else "?", mass if mass != None else "?", inc if inc != None else "?", ecc if ecc != None else "?", name])
                    target_params.append([create_tuple(period,pul,pll), create_tuple(radius,rul,rll), create_tuple(mass,mul,mll), create_tuple(inc,iul,ill), create_tuple(ecc,eul,ell), name])

                targets[n] = tuple(target_params)

            sub_targets = {}

            if mode == "single":
                for key in targets.keys():
                    if key.find(system) != -1 and len(key) == len(system):
                        if len(removed) > 0:
                            sub_targets[key] = tuple(list(targets[key]) + [removed])

                        else:
                            sub_targets[key] = targets[key]

                        break

                return sub_targets

            for key in targets.keys():
                if mode == "tess" and key.find("TOI") != -1 and key.find("test") == -1:
                    sub_targets[key] = targets[key]

                elif mode == "kepler" and key.find("Kepler") != -1 and key.find("test") == -1:
                    sub_targets[key] = targets[key]

                elif mode == "k2" and key.find("K2") != -1 and key.find("test") == -1:
                    sub_targets[key] = targets[key]

                elif mode == "all" and key.find("test") == -1:
                    sub_targets[key] = targets[key]

                elif mode == "test" and key.find("test 3") != -1:
                    sub_targets[key] = targets[key]

            print(sub_targets)
            return sub_targets

        except Exception as e:
            print(e)

        finally:
            dbconn.close()



    def get_isochrones(self, mass):
        """Gets the isochrone data from the database."""

        try:
            dbconn = sqlite3.connect(self.DBNAME)
            dbcursor = dbconn.cursor()
            dbcursor.execute("select smass, age, inner_hz, outer_hz from isochrones where smass = (select smass from isochrones where smass >= " + str(mass) + " order by smass ASC limit 1) order by age")
            iso_upp = dbcursor.fetchall()
            dbcursor.execute("select smass, age, inner_hz, outer_hz from isochrones where smass = (select smass from isochrones where smass <= " + str(mass) + " order by smass DESC limit 1) order by age")
            iso_low = dbcursor.fetchall()

            if iso_upp == []:
                m = [iso_low[0][0], iso_low[0][0]]
                t = [i[1] for i in iso_low]
                ie = [[iso_low[i][2] for i in range(len(iso_low))], [iso_low[i][2] for i in range(len(iso_low))]]
                oe = [[iso_low[i][3] for i in range(len(iso_low))], [iso_low[i][3] for i in range(len(iso_low))]]

            else:
                m = [iso_low[0][0], iso_upp[0][0]]
                t = [i[1] for i in iso_upp]
                ie = [[iso_low[i][2] for i in range(len(iso_upp))], [iso_upp[i][2] for i in range(len(iso_upp))]]
                oe = [[iso_low[i][3] for i in range(len(iso_upp))], [iso_upp[i][3] for i in range(len(iso_upp))]]

            return m, t, ie, oe

        except Exception as e:
            print(e)

        finally:
            dbconn.close()
