from re import VERBOSE
import cctbx
from cctbx import miller
import ast
from mmtbx.scaling import matthews
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from math import isnan


class Phasepred(object):
    """
    docstring here please
    """
    def __init__(self,sg_in,uc_in,s_in,d_min_in=3.0):
        self.sg_in = sg_in
        self.uc_in = uc_in
        self.s_in = s_in
        self.d_min_in = d_min_in


    def colours(self):
        OKGREEN = "\033[92m"
        WARNING = "\033[93m"
        FAIL = "\033[91m"
        BOLD = "\033[1m"
        UNDERLINE = "\033[4m"
        ENDC = "\033[0m"


    def fxn():
        warnings.warn("runtime", RuntimeWarning)


    def predict(self,sg_in,uc_in,asu_mol,d_min_in,s):
        ms = miller.build_set(
            crystal_symmetry=cctbx.crystal.symmetry(
                space_group_symbol=sg_in, unit_cell=(uc_in)
            ),
            anomalous_flag=True,
            d_min=d_min_in,
        )
        refl = int(ms.size())
        ref_per_s = refl / (s * asu_mol)
        return ref_per_s


    def matthewsrupp(self,s_in):
        if s_in.isdigit():
            s = int(s_in)
            return s
        elif type(s_in) == str:
            s = s_in.replace(" ", "")
            s = int(s_in.count("C") + s_in.count("M"))
            s_num = len(s_in)
            asu_pred = matthews.matthews_rupp(
                crystal_symmetry=cctbx.crystal.symmetry(
                    space_group_symbol=sg_in, unit_cell=(uc_in),
                ),
                n_residues=s_num,
            )
            return s_num, asu_pred
        else:
            pass


    def objective_poly(self,x, a, b, c):
        return a * x + b * x ** 2 + c


    def objective_exp(self,x, a, b, c):
        return a * np.exp(-b * x) + c


    def objective_log_find_x(self,y, a, b, c):
        return np.log((y - c) / a) / -b

    def resrange(self):
        res_v_refl = []
        for high_lim in [x / 10.0 for x in range(14, 46, 1)]:
            ref_per_s_theory = predict(sg_in, uc_in, asu_mol, high_lim, s)
            res_v_refl += [(high_lim, ref_per_s_theory)]
        xpred, ypred = zip(*res_v_refl)
        fit_eq, _ = curve_fit(objective_exp, xpred, ypred)
        a, b, c = fit_eq


    def makegraph(self):
        plt.xlabel("d (Å)")
        plt.ylabel("# reflections / anomalous scatterer")
        predictline = plt.plot(*zip(*res_v_refl), label="res-ref (predict)")
        inputblob = plt.scatter(x=d_min_in, y=ref_per_s, c="b")
        plt.annotate(
            "current crystal situation",
            xy=(d_min_in, ref_per_s),
            xytext=(10, 10),
            textcoords="offset pixels",
        )
        redline = 800
        find_redline = objective_log_find_x(redline, a, b, c)
        if isnan(find_redline) == False:
            redline = plt.axhline(
                redline,
                c="r",
                linestyle="--",
                label="borderline = " + str(round(find_redline, 1)) + "Å",
            )
        yellowline = 1100
        find_yellowline = objective_log_find_x(yellowline, a, b, c)
        if isnan(find_yellowline) == False:
            yellowline = plt.axhline(
                yellowline,
                c="y",
                linestyle="--",
                label="acceptable = " + str(round(find_yellowline, 1)) + "Å",
            )
        greenline = 2000
        find_greenline = objective_log_find_x(greenline, a, b, c)
        greenline = plt.axhline(
            greenline,
            c="g",
            linestyle="--",
            label="ideal = " + str(round(find_greenline, 1)) + "Å",
        )
        plt.legend(loc="upper right")
        plt.show()


if __name__ == "__main__":
    sg_in = input("Space Group (# or name): ")
    uc_in = input("Unit Cell (a, b, c, al, be, ga): ")
    d_min_in = float(input("High res: "))
    s_in = input("Input sequence or number of scatterers: ")
    asu_mol = int(input("Number of molecules in the ASU: "))


    predict = Phasepred(sg_in,uc_in,d_min_in,s_in)

    if ref_per_s == 0:
        print(f"{Phasepred.colours.FAIL}\nSomething went wrong...{Phasepred.colours.ENDC}")
    if 0 < ref_per_s < 500:
        print(
            f"{Phasepred.colours.FAIL}\nPhasing is highly unlikely to succeed with this crystal{Phasepred.colours.ENDC}"
        )
        print(
            "\nIf you would like to discuss phasing alternatives, please email armin.wagner@diamond.ac.uk"
        )
    if 500 <= ref_per_s < 800:
        print(
            f"{Phasepred.colours.FAIL}\nPhasing is unlikely to succeed with this crystal{Phasepred.colours.ENDC}"
        )
        print(
            "\nIf you would like to discuss this project, please email armin.wagner@diamond.ac.uk"
        )
    if 800 <= ref_per_s < 1100:
        print(
            f"{Phasepred.colours.OKGREEN}\nPhasing is possible with this crystal, though will be a borderline case{Phasepred.colours.ENDC}"
        )
        print(email)
    if 1100 <= ref_per_s < 2000:
        print(
            f"{Phasepred.colours.OKGREEN}\nPhasing is likely to succeed with this crystal{Phasepred.colours.ENDC}"
        )
        print(email)
    if 2000 <= ref_per_s < 10000:
        print(
            f"{Phasepred.colours.OKGREEN}\nPhasing is highly to succeed with this crystal{Phasepred.colours.ENDC}"
        )
        print(email)
    if 10000 <= ref_per_s:
        print(
            f"{Phasepred.colours.OKGREEN}\nPhasing is essentially guaranteed to succeed with this crystal{Phasepred.colours.ENDC}"
        )
        print(email)