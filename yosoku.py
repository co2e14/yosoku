import cctbx
from cctbx import miller

# Aim is to make a web interface to put crystal details in and give a likely outcome of a project on I23.

ms = miller.build_set(
    crystal_symmetry=crystal.symmetry(
        space_group_symbol="P212121", unit_cell=(6, 6, 6, 90, 90, 90)
    ),
    anomalous_flag=False,
    d_min=3.0,
)
print(ms.size())

print(list(ms.indices()))
print(ms.d_max_min())