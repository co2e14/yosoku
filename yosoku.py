import cctbx
from cctbx import miller
from iotbx.pdb import hierarchy

ms = miller.build_set(
    crystal_symmetry=cctbx.crystal.symmetry(
        space_group_symbol="P321", unit_cell=(150, 150, 45, 90, 90, 120)
    ),
    anomalous_flag=True,
    d_min=1.9,
)

refl = int(ms.size())
print(ms.size())
#print(list(ms.indices()))
print(ms.d_max_min())
print(ms.crystal_symmetry())
print(ms.generate_bivoet_mates())
#print(dir(ms))
#ms.miller_indices_as_pdb_file(file_name='mil.pdb')
#msa = miller.array(ms)
sgs = 0
pdb_in = hierarchy.input(file_name="6fax.pdb")
for chain in pdb_in.hierarchy.only_model().chains():
  for residue_group in chain.residue_groups():
    for atom_group in residue_group.atom_groups():
      for atom in atom_group.atoms():
        if (atom.element.strip().upper() == "S"):
            sgs += 1
print(sgs)
print(refl/sgs)

#pdb_atoms = pdb_in.hierarchy.atoms()
#sel_cache = pdb_in.hierarchy.atom_selection_cache()
#c_alpha_sel = sel_cache.selection("name sg") # XXX not case sensitive!
#c_alpha_atoms = pdb_atoms.select(c_alpha_sel)
#print(len(c_alpha_atoms))