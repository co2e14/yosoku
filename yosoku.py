import cctbx
from cctbx import miller
from mmtbx.scaling import matthews
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from math import isnan
import os
import time

path = os.getcwd()
os.system("unset LIBTBX_BUILD")


class Phasepred(object):
    """
    I23 S-SAD phasing prediction
    """

    def __init__(self, sg_in, uc_in, d_min_in, s_in, outfile=None):
        os.system("unset LIBTBX_BUILD")
        self.sg_in = sg_in
        self.uc_in = uc_in
        self.d_min_in = d_min_in
        self.d_min_static = d_min_in
        self.s_in = s_in
        self.s, self.asu_pred = self.matthewsrupp()
        self.outfile = outfile or os.path.join(path, time.strftime("%Y%m%dT%H%M%S.png"))

    def asu_control(self, asu_in):
        self.asu_mol = int(asu_in)
        self.predict()
        self.ref_per_static = self.ref_per_s
        self.fit_eq = self.resrange()
        self.makegraph()

    def predict(self):
        ms = miller.build_set(
            crystal_symmetry=cctbx.crystal.symmetry(
                space_group_symbol=self.sg_in, unit_cell=(self.uc_in)
            ),
            anomalous_flag=True,
            d_min=float(self.d_min_in),
        )
        refl = int(ms.size())
        self.ref_per_s = refl / (self.s * self.asu_mol)

    def matthewsrupp(self):
        if self.s_in.isdigit():
            s = int(self.s_in)
            return s, None
        elif type(self.s_in) == str:
            s = int(self.s_in.count("C") + self.s_in.count("M"))
            self.s_in = self.s_in.replace(" ", "")
            s_num = len(self.s_in)
            asu_pred = matthews.matthews_rupp(
                crystal_symmetry=cctbx.crystal.symmetry(
                    space_group_symbol=self.sg_in, unit_cell=(self.uc_in),
                ),
                n_residues=s_num,
            )
            return s, asu_pred
        else:
            return None, None

    def objective_exp(self, x, a, b, c):
        return a * np.exp(-b * x) + c

    def objective_log_find_x(self, y, a, b, c):
        return np.log((y - c) / a) / -b

    def resrange(self):
        self.res_v_refl = []
        for high_lim in [x / 10.0 for x in range(14, 46, 1)]:
            self.d_min_in = high_lim
            self.predict()  # make local (return)
            self.res_v_refl += [(high_lim, self.ref_per_s)]
        xpred, ypred = zip(*self.res_v_refl)
        fit_eq, _ = curve_fit(self.objective_exp, xpred, ypred)
        return fit_eq

    def makegraph(self):
        plt.clf()
        a, b, c = self.fit_eq
        plt.xlabel("d (Å)")
        plt.ylabel("# reflections / anomalous scatterer")
        plt.scatter(x=float(self.d_min_static), y=self.ref_per_static, c="b")
        plt.plot(*zip(*self.res_v_refl), label="res-ref (predict)")
        plt.annotate(
            "current crystal situation",
            xy=(float(self.d_min_static), self.ref_per_static),
            xytext=(10, 10),
            textcoords="offset pixels",
        )
        redline = 800
        find_redline = self.objective_log_find_x(redline, *self.fit_eq)
        if isnan(find_redline) == False:
            redline = plt.axhline(
                redline,
                c="r",
                linestyle="--",
                label="borderline = " + str(round(find_redline, 1)) + "Å",
            )
        yellowline = 1100
        find_yellowline = self.objective_log_find_x(yellowline, a, b, c)
        if isnan(find_yellowline) == False:
            yellowline = plt.axhline(
                yellowline,
                c="y",
                linestyle="--",
                label="acceptable = " + str(round(find_yellowline, 1)) + "Å",
            )
        greenline = 2000
        find_greenline = self.objective_log_find_x(greenline, a, b, c)
        greenline = plt.axhline(
            greenline,
            c="g",
            linestyle="--",
            label="ideal = " + str(round(find_greenline, 1)) + "Å",
        )
        plt.legend(loc="upper right")
        plt.savefig(self.outfile)


if __name__ == "__main__":
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"
    ENDC = "\033[0m"
    email = "\nTo receive sample mounts and arrange beamtime (via BAG or rapid access), please email armin.wagner@diamond.ac.uk\n"

    sg_in = input("Space Group (# or name): ")
    uc_in = input("Unit Cell (a, b, c, al, be, ga): ")
    d_min_in = input("High res: ")
    s_in = input("Input sequence or number of scatterers: ")

    predicted = Phasepred(sg_in, uc_in, d_min_in, s_in)
    if predicted.asu_pred is None:
        pass
    else:
        print(predicted.asu_pred.table)
        print(
            "Based on the table above, you likely have",
            predicted.asu_pred.n_copies,
            "molecule(s) in the asu.\n",
        )
    asu_in = input("Number of molecules in the asu: ")
    predicted.asu_control(asu_in)

    if predicted.ref_per_s == 0:
        print(f"{FAIL}\nSomething went wrong...{ENDC}")
    if 0 < predicted.ref_per_s < 500:
        print(f"{FAIL}\nPhasing is highly unlikely to succeed with this crystal{ENDC}")
        print(
            "\nIf you would like to discuss phasing alternatives, please email armin.wagner@diamond.ac.uk"
        )
    if 500 <= predicted.ref_per_s < 800:
        print(f"{FAIL}\nPhasing is unlikely to succeed with this crystal{ENDC}")
        print(
            "\nIf you would like to discuss this project, please email armin.wagner@diamond.ac.uk"
        )
    if 800 <= predicted.ref_per_s < 1100:
        print(
            f"{OKGREEN}\nPhasing is possible with this crystal, though will be a borderline case{ENDC}"
        )
        print(email)
    if 1100 <= predicted.ref_per_s < 2000:
        print(f"{OKGREEN}\nPhasing is likely to succeed with this crystal{ENDC}")
        print(email)
    if 2000 <= predicted.ref_per_s < 10000:
        print(f"{OKGREEN}\nPhasing is highly to succeed with this crystal{ENDC}")
        print(email)
    if 10000 <= predicted.ref_per_s:
        print(
            f"{OKGREEN}\nPhasing is essentially guaranteed to succeed with this crystal{ENDC}"
        )
        print(email)
