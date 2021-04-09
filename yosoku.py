import anvil.server
import anvil.mpl_util
import matplotlib
matplotlib.use('Agg')
import cctbx
from cctbx import miller
from mmtbx.scaling import matthews
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from math import isnan
import os

anvil.server.connect("IYNVTM7ISDIQZIIMXRBDTE5B-VLHZXVIKFTEZ6F35")

path = os.getcwd()
os.system("unset LIBTBX_BUILD")


@anvil.server.callable
def matthewsrupp(sg_in, uc_in, s_in):
    if s_in.isdigit():
        s = int(s_in)
        return s, None
    elif type(s_in) == str:
        s_in = s_in.replace(" ", "")
        s_in = s_in.upper()
        s = int(s_in.count("C") + s_in.count("M"))
        s_num = len(s_in)
        asu_pred = matthews.matthews_rupp(
            crystal_symmetry=cctbx.crystal.symmetry(
                space_group_symbol=sg_in, unit_cell=(uc_in),
            ),
            n_residues=s_num,
        )
        return s, str(asu_pred.table), int(asu_pred.n_copies)
    else:
        return None, None


@anvil.server.callable
def predict(sg_in, uc_in, d_min_in):
    ms = miller.build_set(
        crystal_symmetry=cctbx.crystal.symmetry(
            space_group_symbol=sg_in, unit_cell=(uc_in)
        ),
        anomalous_flag=True,
        d_min=float(d_min_in),
    )
    refl = int(ms.size())
    return refl


@anvil.server.callable
def resrange(sg_in, uc_in):
    res_v_refl = []
    for high_lim in [x / 10.0 for x in range(14, 46, 1)]:
        d_min_in = high_lim
        ref_per_s = predict(sg_in, uc_in, d_min_in)
        res_v_refl += [(high_lim, ref_per_s)]
    xpred, ypred = zip(*res_v_refl)
    fit_eq, _ = curve_fit(objective_exp, xpred, ypred, maxfev=500000)
    return fit_eq, res_v_refl


def objective_exp(x, a, b, c):
    return a * np.exp(-b * x) + c


def objective_log_find_x(y, a, b, c):
    return np.log((y - c) / a) / -b


@anvil.server.callable
def makegraph(fit_eq, d_min_static, ref_per_static, res_v_refl):
    plt.clf()
    a, b, c = fit_eq
    plt.xlabel("d (Å)")
    plt.ylabel("# reflections / anomalous scatterer")
    plt.scatter(x=float(d_min_static), y=ref_per_static, c="b")
    plt.plot(*zip(*res_v_refl), label="res-ref (predict)")
    plt.annotate(
        "current crystal situation",
        xy=(float(d_min_static), ref_per_static),
        xytext=(10, 10),
        textcoords="offset pixels",
    )
    redline = 800
    find_redline = objective_log_find_x(redline, *fit_eq)
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
    print()
    return anvil.mpl_util.plot_image()


if __name__ == "__main__":
    print("yosoku web ver 0.1")
    anvil.server.wait_forever()
    # can use 'supervisor'
