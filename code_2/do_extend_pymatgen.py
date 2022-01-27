import pymatgen.core as mg

#inputfile name (poscar file)
inputfile="POSCAR_STO_unit_cell_relaxed"

#supercell size (axbxc)
extend_x=2
extend_y=2
extend_z=2
extend_supercell=[extend_x,extend_y,extend_z]

#outputfile name 
outputfile="POSCAR_STO_2x2x2"


str_poscar=mg.Structure.from_str(open(inputfile).read(),fmt="poscar")

str_poscar.make_supercell(extend_supercell)

str_poscar.to(filename=outputfile)


