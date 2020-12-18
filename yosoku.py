import cctbx
from cctbx import miller
from iotbx.pdb import hierarchy


class colours:
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"
    ENDC = "\033[0m"


print("\n***** Prediction of SAD Phasing on I23 *****\n")

sg_in = input("Space Group (X123): ")
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

ms = miller.build_set(
    crystal_symmetry=cctbx.crystal.symmetry(
        space_group_symbol=sg_in, unit_cell=(uc_in)
    ),
    anomalous_flag=True,
    d_min=d_min_in,
)
refl = int(ms.size())
ref_per_s = refl / (s * asu_mol)

print("\n***** RESULT *****")
print(
    "\nThere are", s, "sulphur atoms in your protein and", (s * asu_mol), "in the ASU"
)
print("The number of reflections per sulphur atom is", ref_per_s)

email = "\nTo receive sample mounts and arrange beamtime (via BAG or rapid access), please email armin.wagner@diamond.ac.uk"

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
    print("\nIf you would like to discuss this project, please email armin.wagner@diamond.ac.uk")
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
