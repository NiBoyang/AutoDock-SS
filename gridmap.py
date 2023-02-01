import numpy as np

class gridmap_information():

    def __init__(self, target_gridmap_coords, target_atomtype, target_charge,
                 target_coords, target_ligand, target_center, target_size,
                 gridcenter, spacing, npts):

        self.target_gridmap_coords = target_gridmap_coords
        self.target_charge = target_charge
        self.target_atomtype = target_atomtype
        self.target_coords = target_coords
        self.target_ligand = target_ligand
        self.target_center = target_center
        self.target_size = target_size
        self.gridcenter = gridcenter
        self.spacing = spacing
        self.npts = npts

    @classmethod
    def from_gridmap(cls, gridfile, pdbqt):

        with open(gridfile) as tempfile:

            for line in tempfile:
                if 'center' in line:
                    gridcenter = [float(x) for x in ' '.join(line.split()).split(' ')[1:4]]
                elif 'spacing' in line:
                    spacing = float(' '.join(line.split()).split(' ')[1])
                elif 'npts' in line:
                    npts = [int(x) for x in ' '.join(line.split()).split(' ')[1:4]]
        
            xmin = gridcenter[0] - (npts[0] / 2) * spacing
            xcoord = [round(xmin + i * spacing, 4) for i in range(0, npts[0] + 1)]
            ymin = gridcenter[1] - (npts[1] / 2) * spacing
            ycoord = [round(ymin + i * spacing, 4) for i in range(0, npts[1] + 1)]
            zmin = gridcenter[2] - (npts[2] / 2) * spacing
            zcoord = [round(zmin + i * spacing, 4) for i in range(0, npts[2] + 1)]

            target_gridmap_coords = [[xcoord[k], ycoord[j], zcoord[i]]
                    for i in range(0, npts[2] + 1) for j in range(0, npts[1] + 1) for k in range(0, npts[0] + 1)]
        
        with open(pdbqt) as tempfile:

            target_ligand = [' '.join(line.split()).split(' ') for line in tempfile]
            # extract partial charge from pdbqt
            target_charge = [float(x[-2]) for x in target_ligand]
            # extract target ligand's coordinates from pdbqt
            target_coords = np.array([[float(y) for y in x[-7:-4:1]] for x in target_ligand])
            # extract atom types of target ligand from pdbqt
            target_atomtype = list(set([z[-1] for z in target_ligand]))

            x_coords, y_coords, z_coords = target_coords.T

            target_center = [round(sum(x_coords) / len(x_coords), 5),  # simple centre of gravity
                             round(sum(y_coords) / len(y_coords), 5),
                             round(sum(z_coords) / len(z_coords), 5)]

            max_ext = np.array([max(x_coords), max(y_coords), max(z_coords)])  # max extent
            min_ext = np.array([min(x_coords), min(y_coords), min(z_coords)])  # min extent
            target_size = (max_ext - min_ext).round(5)

        return cls(target_gridmap_coords, target_atomtype, target_charge,
                   target_coords, target_ligand, target_center, target_size,
                   gridcenter, spacing, npts)


if __name__ == "__main__":
    gridfile = r"/home/s2263780/lbvs/dude_test/mapfile/ZINC617679grid.txt"
    pdbqtfile = "/home/s2263780/lbvs/dude_test/query_ligmap/pur2_ligand_receptor.pdbqt"
    test = gridmap_information.from_gridmap(gridfile, pdbqtfile)
    print(test.target_size/0.375)