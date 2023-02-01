def affinity_mapping(ligand, folder, spacing, npts, gridcenter, gridcoord, ligandcoord, filename, at):
    import numpy as np
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



def general_map(ligand, folder, spacing, npts, gridcenter, gridcoord, ligandcoord, filename, at):
    with open(f"{folder}{filename}.{at}.map", "w") as gridmap:
        # write the document information into each map file
        gridmap.write("GRID_PARAMETER_FILE " + folder + "\ligand.gpf\n")
        gridmap.write(f"GRID_DATA_FILE {filename}.maps.fld\n")
        gridmap.write(f"MACROMOLECULE {filename}.pdbqt\n")
        gridmap.write("SPACING " + str(spacing) + "\n")
        gridmap.write("NELEMENTS " + str(npts[0]) + " " + str(npts[1]) + " " + str(npts[2]) + "\n")
        gridmap.write("CENTER " + str(gridcenter[0]) + " " + str(gridcenter[1]) + " " + str(gridcenter[2]) + "\n")
        for i in range(0, len(gridcoord)):
            gridmap.write("0" + "\n")



def elec_map(folder, spacing, npts, gridcenter, gridmap_coords, ligand_coords, atomic_partial_charge, filename):
    import numpy as np
    with open(f"{folder}{filename}.e.map", "w") as e_gridmap:
        e_gridmap.write("GRID_PARAMETER_FILE " + folder + "\ligand.gpf\n")
        e_gridmap.write(f"GRID_DATA_FILE {filename}.maps.fld\n")
        e_gridmap.write(f"MACROMOLECULE {filename}.pdbqt\n")
        e_gridmap.write("SPACING " + str(spacing) + "\n")
        e_gridmap.write("NELEMENTS " + str(npts[0]) + " " + str(npts[1]) + " " + str(npts[2]) + "\n")
        e_gridmap.write("CENTER " + str(gridcenter[0]) + " " + str(gridcenter[1]) + " " + str(gridcenter[2]) + "\n")


        for i in gridmap_coords:
            potential, count = 0, 0
            for j in ligand_coords:
                r = np.linalg.norm(np.array(i)-np.array(j))
                loop_potential = atomic_partial_charge[count] / (-0.1465*r)
                potential += loop_potential
                count += 1
            # potential = (potential) * 5
            
            # if potential <= 0:
            #     potential = 0
            e_gridmap.write(f"{potential:.3f}\n")

def des_map(at, ligand, folder, spacing, npts, gridcenter, gridcoord, ligandcoord, filename):
    import numpy as np
    # extract data of an atom type
    # extract the coordiante and convert from string to float
    at_pdbqt = []
    for i in range(0, len(ligand)):
        if ligand[i][-1] == at:
            at_pdbqt.append(ligandcoord[i])

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
            for k in range(0, len(at_pdbqt)):
                if (gridcoord[i][0] - 1.54) <= at_pdbqt[k][0] <= (gridcoord[i][0] + 1.54) \
                        and (gridcoord[i][1] - 1.54) <= at_pdbqt[k][1] <= (gridcoord[i][1] + 1.54) \
                        and (gridcoord[i][2] - 1.54) <= at_pdbqt[k][2] <= (gridcoord[i][2] + 1.54):
                    # grid points collect energy based on distance to atoms
                    x_diff = gridcoord[i][0] - at_pdbqt[k][0]
                    y_diff = gridcoord[i][1] - at_pdbqt[k][1]
                    z_diff = gridcoord[i][2] - at_pdbqt[k][2]
                    dist = np.sqrt(x_diff ** 2 + y_diff ** 2 + z_diff ** 2)
                    energy = round(energy - 10 / dist, 8)

            # add positive value to the grid point with 0 energy
            gridmap.write(str(energy) + "\n")
