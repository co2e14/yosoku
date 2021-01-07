import cctbx
from cctbx import miller
from iotbx import pdb
from iotbx.pdb import hierarchy
import matplotlib.pyplot as plt
import numpy as np
from numpy.core.numeric import NaN
from scipy.optimize import curve_fit
from scipy.special import expit
import warnings
from math import isnan


# import os, libtbx.env_config
# import libtbx.load_env
# print(os.getenv("LIBTBX_BUILD"))
# print("Env:", repr(os.getenv("LIBTBX_BUILD")))
# print(libtbx.env_config.get_installed_path())
# if above gives error like '', run "unset LIBTBX_BUILD" in terminal


class colours:
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"
    ENDC = "\033[0m"


def fxn():
    warnings.warn("runtime", RuntimeWarning)


def predict(sg_in, uc_in, asu_mol, d_min, s, s_type):
    ms = miller.build_set(
        crystal_symmetry=cctbx.crystal.symmetry(
            space_group_symbol=sg_in, unit_cell=(uc_in)
        ),
        anomalous_flag=True,
        d_min=d_min,
    )
    refl = int(ms.size())
    if s_type == "p":
        ref_per_s = refl / s
    else:
        ref_per_s = refl / (s * asu_mol)
    return ref_per_s


def objective_poly(x, a, b, c):
    return a * x + b * x ** 2 + c


def objective_exp(x, a, b, c):
    return a * np.exp(-b * x) + c


def objective_log_find_x(y, a, b, c):
    return np.log((y - c) / a) / -b


with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()

print("\n***** Prediction of SAD Phasing on I23 *****\n")

# sg_in = "P321"
# uc_in = "150 150 45 90 90 120"
# asu_mol = 1
# d_min_in = 2.6
# s = 14

sg_in = input("Space Group (# or name): ")
uc_in = input("Unit Cell (a, b, c, al, be, ga): ")
asu_mol = int(input("Number of molecules in the ASU: "))
d_min_in = float(input("High res: "))
s_type = input("Supply (p)db file, (s)equence, or (n)umber of scatterers: ")

if s_type == "p":
    pdb_file = str(input("Path to PDB file: "))
    s = 0
    pdb_in = hierarchy.input(file_name="6fax.pdb")
    for chain in pdb_in.hierarchy.only_model().chains():
        for residue_group in chain.residue_groups():
            for atom_group in residue_group.atom_groups():
                for atom in atom_group.atoms():
                    if atom.element.strip().upper() == "S":
                        s += 1
if s_type == "n":
    s = int(input("Number of S atoms: "))
if s_type == "s":
    s_in = str((input("Sequence (letters only): ")).replace(" ", ""))
    s = s_in.count("C") + s_in.count("M")

ref_per_s = predict(sg_in, uc_in, asu_mol, d_min_in, s, s_type)

print("\n***** RESULT *****")
print(
    "\nThere are", s, "sulphur atoms in your protein and", (s * asu_mol), "in the ASU"
)
print("The number of reflections per sulphur atom is", ref_per_s)

email = "\nTo receive sample mounts and arrange beamtime (via BAG or rapid access), please email armin.wagner@diamond.ac.uk\n"

if ref_per_s == 0:
    print(f"{colours.FAIL}\nSomething went wrong...{colours.ENDC}")
if 0 < ref_per_s < 500:
    print(
        f"{colours.FAIL}\nPhasing is highly unlikely to succeed with this crystal{colours.ENDC}"
    )
    print(
        "\nIf you would like to discuss phasing alternatives, please email armin.wagner@diamond.ac.uk"
    )
if 500 <= ref_per_s < 800:
    print(
        f"{colours.FAIL}\nPhasing is unlikely to succeed with this crystal{colours.ENDC}"
    )
    print(
        "\nIf you would like to discuss this project, please email armin.wagner@diamond.ac.uk"
    )
if 800 <= ref_per_s < 1100:
    print(
        f"{colours.OKGREEN}\nPhasing is possible with this crystal, though will be a borderline case{colours.ENDC}"
    )
    print(email)
if 1100 <= ref_per_s < 2000:
    print(
        f"{colours.OKGREEN}\nPhasing is likely to succeed with this crystal{colours.ENDC}"
    )
    print(email)
if 2000 <= ref_per_s < 10000:
    print(
        f"{colours.OKGREEN}\nPhasing is highly to succeed with this crystal{colours.ENDC}"
    )
    print(email)
if 10000 <= ref_per_s:
    print(
        f"{colours.OKGREEN}\nPhasing is essentially guaranteed to succeed with this crystal{colours.ENDC}"
    )
    print(email)

res_v_refl = []
for high_lim in [x / 10.0 for x in range(14, 46, 1)]:
    ref_per_s_theory = predict(sg_in, uc_in, asu_mol, high_lim, s, s_type)
    res_v_refl += [(high_lim, ref_per_s_theory)]

xpred, ypred = zip(*res_v_refl)
fit_eq, _ = curve_fit(objective_exp, xpred, ypred)
a, b, c = fit_eq
print("\nThe equation for this crystal is: y = %.5f e (-%.5fx) + %.5f" % (a, b, c))


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
    print(
        "\nPhasing is essentially impossible if your crystal does not diffract to at least %.2fÅ"
        % (find_redline)
    )
    redline = plt.axhline(
        redline,
        c="r",
        linestyle="--",
        label="borderline = " + str(round(find_redline, 1)) + "Å",
    )
yellowline = 1100
find_yellowline = objective_log_find_x(yellowline, a, b, c)
if isnan(find_yellowline) == False:
    print(
        "For a chance at solving, you need a crystal which diffracts to %.2fÅ"
        % (find_yellowline)
    )
    yellowline = plt.axhline(
        yellowline,
        c="y",
        linestyle="--",
        label="acceptable = " + str(round(find_yellowline, 1)) + "Å",
    )
greenline = 2000
find_greenline = objective_log_find_x(greenline, a, b, c)
print(
    "For a very good chance at solving, you need a crystal which diffracts to %.2fÅ"
    % (find_greenline)
)
greenline = plt.axhline(
    greenline,
    c="g",
    linestyle="--",
    label="ideal = " + str(round(find_greenline, 1)) + "Å",
)
plt.legend(loc="upper right")
plt.show()
