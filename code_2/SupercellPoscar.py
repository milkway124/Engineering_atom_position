#! /usr/bin/python
#############################################################
#Writer Jaewook Lee, (2021-07-07)
#   Make the Superceel POSCAR file
#           2021.07.07
#############################################################
import sys
import os

if __name__=="__main__":

    if len(sys.argv) != 6 :
        print("This code make supercell using POSCAR file!! (need pymatgen!!)")
        print("[1] : input POSCAR file, [2] output POSCAR file, [3] scaling x-factor, [4] scaling y-factor, [5] scaling z-factor")
        sys.exit()

    if not (os.path.isfile(sys.argv[1])) :
        print("There is no POSCAR file!!")
        print("[1] : ",sys.argv[1])
        sys.exit() 

    inputfile=sys.argv[1]
    outputfile=sys.argv[2]

    extend_supercell=[sys.argv[3],sys.argv[4],sys.argv[5]]

    try:
        import pymatgen.core as mg
    except:
        print("you need the 'pymatgen' module!!")
        print("Now install pymatgen!!")
        print("\n'    conda install pymatgen  '\n")
        sys.exit()

    #read poscar file
    str_poscar=mg.Structure.from_str(open(inputfile).read(),fmt="poscar")
    
    #make supercell
    str_poscar.make_supercell(extend_supercell)

    """
     scaling_matrix: A scaling matrix for transforming the lattice
         vectors. Has to be all integers. Several options are possible:
         a. A full 3x3 scaling matrix defining the linear combination
            the old lattice vectors. E.g., [[2,1,0],[0,3,0],[0,0,
            1]] generates a new structure with lattice vectors a' =
            2a + b, b' = 3b, c' = c where a, b, and c are the lattice
            vectors of the original structure.
         b. An sequence of three scaling factors. E.g., [2, 1, 1]
            specifies that the supercell should have dimensions 2a x b x
            c.
         c. A number, which simply scales all lattice vectors by the
            same factor.
    """
    #make new poscar file
    try:
        str_poscar.to(filename=outputfile)
    except:
        print("The format of Outputfile file is weird!!")
        print("ex)output = POSCAR~~~")
        sys.exit()
    print("Done make supercell :",outputfile)

