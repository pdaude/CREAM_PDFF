import numpy as np

class Custom_FWspectrum:
    def __init__(self, **kwargs):
        self.initFV()
        for key, value in kwargs.items():
            if key in self.FVs:
                self.FVs[key] = value
        self.FVs["realAmps"] = np.array(self.setFattyPeaks())
    def getFatFreq(self):
        return self.FVs["fatCS"]

    def getAlpha(self):
        alpha = np.zeros([2 + self.FVs["UFV"], len(self.FVs["fatCS"]) + 1], dtype=np.float32)
        alpha[0, 0] = 1.  # Water component
        alpha[1:, 1:] = self.FVs["realAmps"]
        return alpha

    def initFV(self):
        self.FVs = {"UFV": 0,
                    "fatCS": [],
                    "realAmps": []}

    def setFattyPeaks(self):
        return self.FVs["realAmps"]

    def getUFV(self):
        return self.FVs["UFV"]

class ISMRMChallenge_FWspectrum(Custom_FWspectrum):
    # Hamilton G, et al. NMR Biomed. 24(7):784-90, 2011. PMID: 21834002
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def initFV(self):
        self.FVs = {"UFV": 0,
                    "fatCS": [0.90, 1.30, 2.1, 2.76, 4.31, 5.3] ,
                    "realAmps": np.array([87, 693, 128, 4, 39, 48]) / 1000}

class Berglund2012_FWspectrum(Custom_FWspectrum):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def initFV(self):
        self.FVs = {"UFV": 0,
                    "CL": 17.4,
                    "P2U": 0.2,
                    "UD": 2.6,
                    "fatCS": [0.90, 1.30, 1.59, 2.03, 2.25, 2.77, 4.1, 4.3, 5.21, 5.31],
                    "realAmps": []}

    def setFattyPeaks(self):
        CL = self.FVs["CL"]
        UD = self.FVs["UD"]
        P2U = self.FVs["P2U"]
        amps = 0
        if self.FVs["UFV"] == 0:
            amps = np.array([9, (CL - 4) * 6 + UD * (2 * P2U - 8), 6, 4 * (UD * (1 - P2U)), 6, 2 * UD * P2U, 2, 2, 1, 2 * UD])
            amps = amps / np.sum(amps)

        elif self.FVs["UFV"] == 1:
            # F1 = 9A+6(CL-4)B+6C+6E+2G+2H+I
            # F2 = (2P2U-8)B+4(1-P2U)D+2P2UF+2J
            amps = np.array([[9, 6 * (CL - 4), 6, 0, 6, 0, 2, 2, 1, 0],
                    [0, 2 * P2U - 8, 0, 4 * (1 - P2U), 0, 2 * P2U, 0, 0, 0, 2]])
        elif self.FVs["UFV"] == 2:
            # F1 = 9A+6(CL-4)B+6C+6E+2G+2H+I
            # F2 = -8B+4D+2J
            # F3 = 2B-4D+2F
            amps = np.array([[9, 6 * (CL - 4), 6, 0, 6, 0, 2, 2, 1, 0],
                   [0, -8, 0, 4, 0, 0, 0, 0, 0, 2],
                    [0, 2, 0, -4, 0, 2, 0, 0, 0, 0]])
        elif self.FVs["UFV"] == 3:
            # F1 = 9A-24B+6C+6E+2G+2H+I
            # F2 = -8B+4D+2J
            # F3 = 2B-4D+2F
            # F4 = 6B
            amps = np.array([[9, -24, 6, 0, 6, 0, 2, 2, 1, 0],
                    [0, -8, 0, 4, 0, 0, 0, 0, 0, 2],
                    [0, 2, 0, -4, 0, 2, 0, 0, 0, 0],
                    [0, 6, 0, 0, 0, 0, 0, 0, 0, 0]])
        else: raise Exception("Ukwnown Fat variable should be between 0 and 3")

        return amps


class HamiltonLiver2011_FWspectrum(Custom_FWspectrum):
    # Hamilton G, et al. NMR Biomed. 24(7):784-90, 2011. PMID: 21834002
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def initFV(self):
        self.FVs = {"UFV": 0,
                    "fatCS": [0.90, 1.30, 1.60, 2.02, 2.24, 2.75, 4.20, 5.19, 5.29],
                    "realAmps": np.array([88, 642, 58, 62, 58, 6, 39, 10, 37]) / 1000}


class Hamiltonndb2011_FWspectrum(Custom_FWspectrum):
    # Hamilton G, et al. NMR Biomed. 24(7):784-90, 2011. PMID: 21834002
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def initFV(self):
        self.FVs = {"UFV": 0,
                    "NDB": 0,
                    "NMIDB": 0,
                    "CL": 0,
                    "fatCS": [0.90, 1.30, 1.60, 2.02, 2.24, 2.75, 4.20, 5.19, 5.29],
                    "realAmps": []}

    def setFattyPeaks(self):
        CL = self.FVs["CL"]
        nmidb = self.FVs["NMIDB"]
        ndb = self.FVs["NDB"]
        if self.FVs["UFV"] == 0:
            amps = np.array([9, ((CL - 4) * 6) - (ndb * 8) + nmidb * 2, 6, (ndb - nmidb) * 4, 6, nmidb * 2, 4, 1, 2 * ndb])
            amps = amps / np.sum(amps)

        elif self.FVs["UFV"] == 1:
            #Unkwown Fat variable = ndb
            # F1 = 9A+(6(CL-4)+2nmidb)B+6C-4(nmidb)D+6E+2(nmidb)F+4(G+H)+I
            # F2 = -8B+4D+2J
            amps = np.array([[9, 6 * (CL - 4)+2*nmidb, 6, nmidb*-4, 6, nmidb * 2, 4, 1, 0],
                             [0,- 8, 0,4, 0,0, 0, 0, 2]])
        elif self.FVs["UFV"] == 2:
            # Unkwown Fat variable = ndb and nmidb
            # F1 = 9A+(6(CL-4)B+6C+6E+4(G+H)+I
            # F2 = -8B+4D+2J
            # F3 = 2B -4D+2F
            amps = np.array([[9, 6 * (CL - 4), 6,0, 6,0, 4, 1, 0],
                             [0, - 8, 0, 4, 0, 0, 0, 0, 2],
                             [0, 2, 0, -4, 0, 2, 0, 0, 0]])

        elif self.FVs["UFV"] == 3:
            # Unkwown Fat variable = ndb, nmidb and CL
            # F1 = 9A-24B+6C+6E+4(G+H)+I
            # F2 = -8B+4D+2J
            # F3 = 2B -4D+2F
            # F4 = 6B
            amps = np.array([[9, -24, 6, 0, 6, 0, 4, 1, 0],
                             [0, - 8, 0, 4, 0, 0, 0, 0, 2],
                             [0, 2, 0, -4, 0, 2, 0, 0, 0],
                             [0, 6, 0, 0, 0, 0, 0, 0, 0]])
        else:
            raise Exception("Ukwnown Fat variable should be between 0 and 3")
        return amps


class Bydder2011_FWspectrum(Hamiltonndb2011_FWspectrum):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.FVs["NMIDB"] = 0.093 * self.FVs["NDB"] ** 2
        self.FVs["CL"] = 16.8 + 0.25 * self.FVs["NDB"]

    def initFV(self):
        self.FVs = {"UFV": 0,
                    "NDB": 2.5,
                    "NMIDB": 0,
                    "CL": 0,
                    "fatCS": [0.90, 1.30, 1.60, 2.02, 2.24, 2.75, 4.20, 5.19, 5.29],
                    "realAmps": []}


class Hodson2008_FWspectrum(Hamiltonndb2011_FWspectrum):
    # Hodson L, et al. Prog Lipid Res. 2008;47:348-380. PMID: 18435934 DOI: 10.1016/j.plipres.2008.03.003
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def initFV(self):
        self.FVs = {"UFV": 0,
                    "NDB": 2.69,
                    "NMIDB": 0.58,
                    "CL": 17.29,
                    "fatCS": [0.90, 1.30, 1.60, 2.02, 2.24, 2.75, 4.20, 5.19, 5.29],
                    "realAmps": []}


class Leporq2014SAT_FWspectrum(Hamiltonndb2011_FWspectrum):
    # Leporq B, et al. NMR Biomed. 2014; 27: 1211–1221 DOI:10.1002/nbm.3175
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def initFV(self):
        self.FVs = {"UFV": 0,
                    "NDB": 2.53,
                    "NMIDB": 0.94,
                    "CL": 17.47,
                    "fatCS": [0.90, 1.30, 1.60, 2.02, 2.24, 2.75, 4.20, 5.19, 5.29],
                    "realAmps": []
                    }


class Leporq2014VAT_FWspectrum(Hamiltonndb2011_FWspectrum):
    # Leporq B, et al. NMR Biomed. 2014; 27: 1211–1221 DOI:10.1002/nbm.3175
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def initFV(self):
        self.FVs = {"UFV": 0,
                    "NDB": 2.72,
                    "NMIDB": 0.84,
                    "CL": 17.43,
                    "fatCS": [0.90, 1.30, 1.60, 2.02, 2.24, 2.75, 4.20, 5.19, 5.29],
                    "realAmps": []}


