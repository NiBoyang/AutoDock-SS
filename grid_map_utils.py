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
        
class GridMapInformation:
    """
    A class to represent grid map information for a target ligand.
    """

    def __init__(self, target_gridmap_coords, target_atomtype, target_charge,
                 target_coords, target_ligand, target_center, target_size,
                 grid_center, spacing, npts):
        """
        Initialize the GridMapInformation object with the given parameters.
        """
        self.target_gridmap_coords = target_gridmap_coords
        self.target_charge = target_charge
        self.target_atomtype = target_atomtype
        self.target_coords = target_coords
        self.target_ligand = target_ligand
        self.target_center = target_center
        self.target_size = target_size
        self.grid_center = grid_center
        self.spacing = spacing
        self.npts = npts

    @classmethod
    def from_gridmap(cls, gridfile, pdbqt):
        """
        Create a GridMapInformation object from the given grid file and pdbqt file.
        """
        grid_center, spacing, npts = cls._parse_gridfile(gridfile)
        target_gridmap_coords = cls._generate_target_gridmap_coords(grid_center, spacing, npts)
        (target_ligand, target_charge, target_coords,
         target_atomtype, target_center, target_size) = cls._parse_pdbqt(pdbqt)

        return cls(target_gridmap_coords, target_atomtype, target_charge,
                   target_coords, target_ligand, target_center, target_size,
                   grid_center, spacing, npts)

    @staticmethod
    def _parse_gridfile(gridfile):
        with open(gridfile) as tempfile:
            for line in tempfile:
                if 'center' in line:
                    grid_center = [float(x) for x in ' '.join(line.split()).split(' ')[1:4]]
                elif 'spacing' in line:
                    spacing = float(' '.join(line.split()).split(' ')[1])
                elif 'npts' in line:
                    npts = [int(x) for x in ' '.join(line.split()).split(' ')[1:4]]
        return grid_center, spacing, npts

    @staticmethod
    def _generate_target_gridmap_coords(grid_center, spacing, npts):
        xmin = grid_center[0] - (npts[0] / 2) * spacing
        xcoord = [round(xmin + i * spacing, 4) for i in range(0, npts[0] + 1)]
        ymin = grid_center[1] - (npts[1] / 2) * spacing
        ycoord = [round(ymin + i * spacing, 4) for i in range(0, npts[1] + 1)]
        zmin = grid_center[2] - (npts[2] / 2) * spacing
        zcoord = [round(zmin + i * spacing, 4) for i in range(0, npts[2] + 1)]

        target_gridmap_coords = [[xcoord[k], ycoord[j], zcoord[i]]
                                 for i in range(0, npts[2] + 1) for j in range(0, npts[1] + 1) for k in range(0, npts[0] + 1)]
        return target_gridmap_coords

    @staticmethod
    def _parse_pdbqt(pdbqt):
        with open(pdbqt) as tempfile:
            target_ligand = [' '.join(line.split()).split(' ') for line in tempfile]
            target_charge = [float(x[-2]) for x in target_ligand]
            target_coords = np.array([[float(y) for y in x[-7:-4:1]] for x in target_ligand])
            target_atomtype = list(set([z[-1] for z in target_ligand]))

            x_coords, y_coords, z_coords = target_coords.T

            target_center = [round(sum(x_coords) / len(x_coords), 5),
                             round(sum(y_coords) / len(y_coords), 5),
                             round(sum(z_coords) / len(z_coords), 5)]

            max_ext = np.array([max(x_coords), max(y_coords), max(z_coords)])
            min_ext = np.array([min(x_coords), min(y_coords), min(z_coords)])
            target_size = (max_ext - min_ext).round(5)

        return target_ligand, target_charge, target_coords, target_atomtype, target_center, target_size

def affinity_map(ligand, folder, spacing, npts, gridcenter, gridcoord, ligandcoord, filename, at):
    at_pdbqt = np.array([ligandcoord[i] for i, atom in enumerate(ligand) if atom[-1] == at])

    # Generate map file of the specified atom type
    with open(f"{folder}{filename}.{at}.map", "w") as gridmap:

        # Write the document information into each map file
        gridmap.write("GRID_PARAMETER_FILE " + folder + "\ligand.gpf\n")
        gridmap.write(f"GRID_DATA_FILE {filename}.maps.fld\n")
        gridmap.write(f"MACROMOLECULE {filename}.pdbqt\n")
        gridmap.write("SPACING " + str(spacing) + "\n")
        gridmap.write("NELEMENTS " + str(npts[0]) + " " + str(npts[1]) + " " + str(npts[2]) + "\n")
        gridmap.write("CENTER " + str(gridcenter[0]) + " " + str(gridcenter[1]) + " " + str(gridcenter[2]) + "\n")

        # Vectorized distance calculation and energy computation
        energy_values = []
        threshold = 1 * 1.54
        for grid in gridcoord:
            grid = np.array(grid)
            distances = np.linalg.norm(at_pdbqt - grid, axis=1)
            dist_lis = distances[distances <= threshold]

            if len(dist_lis) > 0:
                sorted_dist = np.sort(dist_lis)[::-1]
                energy = 0
                for d in sorted_dist:
                    energy = 0.225 * (energy - 10) / np.sqrt(d)
                energy_values.append(energy)
            else:
                energy_values.append(0)

        # Write energy values to the file
        gridmap.write("\n".join(f"{energy:.3f}" for energy in energy_values))

def general_map(ligand, folder, spacing, npts, gridcenter, gridcoord, ligandcoord, filename, at):
    with open(f"{folder}{filename}.{at}.map", "w") as gridmap:
        gridmap.write("GRID_PARAMETER_FILE " + folder + "\ligand.gpf\n")
        gridmap.write(f"GRID_DATA_FILE {filename}.maps.fld\n")
        gridmap.write(f"MACROMOLECULE {filename}.pdbqt\n")
        gridmap.write("SPACING " + str(spacing) + "\n")
        gridmap.write("NELEMENTS " + str(npts[0]) + " " + str(npts[1]) + " " + str(npts[2]) + "\n")
        gridmap.write("CENTER " + str(gridcenter[0]) + " " + str(gridcenter[1]) + " " + str(gridcenter[2]) + "\n")

        for i in range(len(gridcoord)):
            gridmap.write("0.00\n")

def elec_map(folder, spacing, npts, gridcenter, gridmap_coords, ligand_coords, atomic_partial_charge, filename):
    with open(f"{folder}{filename}.e.map", "w") as e_gridmap:
        # Write the header information
        e_gridmap.write("GRID_PARAMETER_FILE " + folder + "\ligand.gpf\n")
        e_gridmap.write(f"GRID_DATA_FILE {filename}.maps.fld\n")
        e_gridmap.write(f"MACROMOLECULE {filename}.pdbqt\n")
        e_gridmap.write("SPACING " + str(spacing) + "\n")
        e_gridmap.write("NELEMENTS " + str(npts[0]) + " " + str(npts[1]) + " " + str(npts[2]) + "\n")
        e_gridmap.write("CENTER " + str(gridcenter[0]) + " " + str(gridcenter[1]) + " " + str(gridcenter[2]) + "\n")

        # Convert inputs to NumPy arrays
        gridmap_coords_array = np.array(gridmap_coords)
        ligand_coords_array = np.array(ligand_coords)
        atomic_partial_charge_array = np.array(atomic_partial_charge)

        # Initialize an array to store potentials
        potentials = []

        # Vectorized calculation of potentials
        for grid_point in gridmap_coords_array:
            r = np.linalg.norm(grid_point - ligand_coords_array, axis=1)
            with np.errstate(divide='ignore'):  # Ignore division by zero warnings
                loop_potential = atomic_partial_charge_array / (-0.1465 * r)
            potential = np.nansum(loop_potential)  # Sum, ignoring NaNs from division by zero
            potentials.append(f"{potential:.3f}")

        # Write all potentials at once
        e_gridmap.write("\n".join(potentials))
