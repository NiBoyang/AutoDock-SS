import os
import numpy as np

def create_mapfld(folder, filename, spacing, npts, center):

    with open(f"{folder}{filename}.maps.fld", "w") as mapfld:

        mapfld.write(f"""
# AVS field file
#
# AutoDock Atomic Affinity and Electrostatic Grids
#
# Created by ./../autogrid4.
#
#SPACING {spacing}
#NELEMENTS {npts[0]} {npts[1]} {npts[-1]}
#CENTER {center[0]:.3f}  {center[1]:.3f} {center[-1]:.3f}
#MACROMOLECULE ZINC000000002212_receptor.pdbqt
#GRID_PARAMETER_FILE D:\AutoDock_Workplace\script test\..\scripttest.gpf
#
ndim=3			# number of dimensions in the field
dim1={npts[0]+1}			# number of x-elements
dim2={npts[1]+1}			# number of y-elements
dim3={npts[-1]+1}			# number of z-elements
nspace=3		# number of physical coordinates per point
veclen=9		# number of affinity values at each point
data=float		# data type (byte, integer, float, double)
field=uniform		# field type (uniform, rectilinear, irregular)
coord 1 file=ligand.maps.xyz filetype=ascii offset=0
coord 2 file=ligand.maps.xyz filetype=ascii offset=2
coord 3 file=ligand.maps.xyz filetype=ascii offset=4
label=A-affinity	# component label for variable 1
label=C-affinity	# component label for variable 2
label=HD-affinity	# component label for variable 3
label=N-affinity	# component label for variable 4
label=NA-affinity	# component label for variable 5
label=OA-affinity	# component label for variable 6
label=F-affinity
label=P-affinity
label=SA-affinity
label=S-affinity
label=Cl-affinity
label=Br-affinity
label=I-affinity
label=Electrostatics
label=Desolvation
#
# location of affinity grid files and how to read them
#
variable 1 file={filename}.A.map filetype=ascii skip=6
variable 2 file={filename}.C.map filetype=ascii skip=6
variable 3 file={filename}.HD.map filetype=ascii skip=6
variable 4 file={filename}.N.map filetype=ascii skip=6
variable 5 file={filename}.NA.map filetype=ascii skip=6
variable 6 file={filename}.OA.map filetype=ascii skip=6
variable 7 file={filename}.F.map filetype=ascii skip=6
variable 8 file={filename}.P.map filetype=ascii skip=6
variable 9 file={filename}.SA.map filetype=ascii skip=6
variable 10 file={filename}.S.map filetype=ascii skip=6
variable 11 file={filename}.Cl.map filetype=ascii skip=6
variable 12 file={filename}.Br.map filetype=ascii skip=6
variable 13 file={filename}.I.map filetype=ascii skip=6
variable 14 file={filename}.e.map filetype=ascii skip=6
variable 15 file={filename}.d.map filetype=ascii skip=6
""")

def affinity_map(ligand, folder, spacing, npts, gridcenter, gridcoord, ligandcoord, filename, at):
    import numpy
    # halogens = ["F", "Cl", "Br", "I"]
    # aromatic_atoms = ["A", "NA", "OA", "SA"]
    # if at in halogens:
    #     compensation = set(halogens)-set(at)
    # elif at in aromatic_atoms:
    #     compensation = set(aromatic_atoms)-set(at)
    # else:
    #     compensation = [None]
    # extract data of an atom type
    # extract the coordiante and convert from string to float
    at_pdbqt = []
    # compensation_pdbqt = []

    for i in range(0, len(ligand)):
        if ligand[i][-1] == at:
            at_pdbqt.append(ligandcoord[i])
        # elif ligand[i][-1] in compensation:
        #     compensation_pdbqt.append(ligandcoord[i])

    # generate map file of carbon
    with open(f"{folder}{filename}.{at}.map", "w") as gridmap:

        # write the document information into each map file
        gridmap.write("GRID_PARAMETER_FILE " + folder + "\ligand.gpf\n")
        gridmap.write(f"GRID_DATA_FILE {filename}.maps.fld\n")
        gridmap.write(f"MACROMOLECULE {filename}.pdbqt\n")
        gridmap.write("SPACING " + str(spacing) + "\n")
        gridmap.write("NELEMENTS " + str(npts[0]) + " " + str(npts[1]) + " " + str(npts[2]) + "\n")
        gridmap.write("CENTER " + str(gridcenter[0]) + " " + str(gridcenter[1]) + " " + str(gridcenter[2]) + "\n")

        # if the carbon atom is within 5 angstrom of a grid point
        # energy will be counted into that grid point based on the distance
        for i in range(0, len(gridcoord)):
            energy = 0
            dist_lis = []
            for k in range(0, len(at_pdbqt)):
            #     if (gridcoord[i][0] - 1.54) <= at_pdbqt[k][0] <= (gridcoord[i][0] + 1.54)\
            #         and (gridcoord[i][1] - 1.54) <= at_pdbqt[k][1] <= (gridcoord[i][1] + 1.54)\
            #             and (gridcoord[i][2] - 1.54) <= at_pdbqt[k][2] <= (gridcoord[i][2] + 1.54):
                if (gridcoord[i][0]-at_pdbqt[k][0])**2+ \
                    (gridcoord[i][1]-at_pdbqt[k][1])**2+ \
                    (gridcoord[i][2]-at_pdbqt[k][2])**2 <= (1*1.54)**2:

                    dist = np.linalg.norm(np.array(gridcoord[i]) - np.array(at_pdbqt[k]))
                    # energy = 0.3125*(energy - 10) / (dist**0.5)
                    
                    dist_lis.append(dist)

            if len(dist_lis) != 0:
                dist_lis = np.array(dist_lis)
                sorted_ind = np.argsort(dist_lis)
                sorted_dist = dist_lis[sorted_ind[::-1]]
                for d in sorted_dist:
                    energy = 0.225*(energy - 10) / (d**0.5)
    
            # add positive value to the grid point with 0 energy
            # if energy == 0:
            #     dist1 = [np.linalg.norm(np.array(gridcoord[i]) - np.array(at_pdbqt[j])) for j in range(0, len(at_pdbqt))]
            #     energy = min(dist1)

            # if 0 <= energy <= 50:
            #     energy = 0.029 * energy ** 2.5 + 1

            gridmap.write(f"{energy:.3f}\n")

def general_map(at, target_ligand, folder, spacing, npts, gridcenter, gridmap_coords, target_coords, filename):
    with open(os.path.join(folder, f"{filename}.{at}.map"), "w") as gridmap:
        gridmap.write(f"GRID_PARAMETER_FILE\n\n")
        gridmap.write(f"GRID_CENTER {gridcenter[0]:.3f} {gridcenter[1]:.3f} {gridcenter[2]:.3f}\n")
        gridmap.write(f"GRID_SIZE {npts[0]} {npts[1]} {npts[2]}\n")
        gridmap.write(f"GRID_SPACING {spacing}\n\n")
        gridmap.write(f"MACROMOLECULE {filename}.pdbqt\n\n")
        gridmap.write(f"GRID_DATA_FILE {filename}.{at}.map\n\n")
        gridmap.write(f"ATOM_TYPE {at}\n\n")

        for i in range(len(gridmap_coords)):
            gridmap.write("0.00\n")

def elec_map(folder, spacing, npts, gridcenter, gridmap_coords, ligand_coords, atomic_partial_charge, filename):
    with open(os.path.join(folder, f"{filename}.e.map"), "w") as e_gridmap:
        e_gridmap.write("GRID_PARAMETER_FILE " + folder + "\ligand.gpf\n")
        e_gridmap.write(f"GRID_DATA_FILE {filename}.maps.fld\n")
        e_gridmap.write(f"MACROMOLECULE {filename}.pdbqt\n")
        e_gridmap.write("SPACING " + str(spacing) + "\n")
        e_gridmap.write("NELEMENTS " + str(npts[0]) + " " + str(npts[1]) + " " + str(npts[2]) + "\n")
        e_gridmap.write("CENTER " + str(gridcenter[0]) + " " + str(gridcenter[1]) + " " + str(gridcenter[2]) + "\n")

        gridmap_coords_array = np.array(gridmap_coords)
        ligand_coords_array = np.array(ligand_coords)
        atomic_partial_charge_array = np.array(atomic_partial_charge)

        for i in gridmap_coords_array:
            r = np.linalg.norm(i - ligand_coords_array, axis=1)
            loop_potential = atomic_partial_charge_array / (-0.1465 * r)
            potential = loop_potential.sum()
            e_gridmap.write(f"{potential:.3f}\n")
