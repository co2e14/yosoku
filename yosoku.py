import cctbx
from cctbx import miller


# Aim is to make a web interface to put crystal details in and give a likely outcome of a project on I23.

ms = miller.build_set(
    crystal_symmetry=cctbx.crystal.symmetry(
        space_group_symbol="P321", unit_cell=(150, 150, 45, 90, 90, 120)
    ),
    anomalous_flag=True,
    d_min=1.9,
)
print(ms.size())
#print(list(ms.indices()))
print(ms.d_max_min())
print(ms.crystal_symmetry())
print(ms.minimum_wavelength_based_on_d_min())
print(ms.generate_bivoet_mates())
#print(ms.show_completeness(2))
print(ms.minimum_wavelength_based_on_d_min())
print(ms.minimum_wavelength_based_on_d_min())
print(dir(ms))
ms.miller_indices_as_pdb_file(file_name='mil.pdb')