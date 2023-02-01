import argparse

def map2pdb(mapfile, gpffile):
    with open(gpffile) as tempfile:
        for line in tempfile.readlines():
            if "spacing" in line:
                spacing = float(' '.join(line.split()).split(' ')[1])
            elif "gridcenter" in line:
                gridcenter = [float(x) for x in ' '.join(line.split()).split(' ')[1:4]]
            elif "npts" in line:
                ngridpoints = [int(x) for x in ' '.join(line.split()).split(' ')[1:4]]

    with open(mapfile) as tempfile:
        energies = [float(x.strip("\n")) for x in tempfile.readlines()[6:]]

    xmin = gridcenter[0] - ((ngridpoints[0])/2)*spacing
    ymin = gridcenter[1] - ((ngridpoints[1])/2)*spacing
    zmin = gridcenter[2] - ((ngridpoints[2])/2)*spacing
    atomno, x, y, z = 0, 0, 0, 0

    with open(f"{mapfile}.pdb", "w") as pdbfile:
        for energy in energies:
            while z < ngridpoints[2]+1:
                zcoord = zmin + z*spacing
                z += 1
                while y < ngridpoints[1]+1:
                    ycoord = ymin + y*spacing
                    y += 1
                    while x < ngridpoints[0]+1:
                        xcoord = xmin + x*spacing
                        if energies[atomno] <= -0.35:
                            pdbfile.write(f"ATOM{atomno:>7}  C   XXX     1    {xcoord:8.3f}{ycoord:8.3f}{zcoord:8.3f}  1.00{energies[atomno]:6.2f}\n")
                        atomno += 1
                        x += 1
                    x = 0
                y=0

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='map2pdb')
    parser.add_argument('mapfile')
    parser.add_argument('gpffile')
    args = parser.parse_args()

    mapfile = args.mapfile
    gpffile = args.gpffile

    map2pdb(mapfile, gpffile)