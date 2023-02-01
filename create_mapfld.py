from pathlib import Path

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

if __name__ == "__main__":

    pdbqt_path = "/root/autodl-tmp/lbvs/all_dude/aa2ar/docking_test/lig/lig.pdbqt"
    path = "/root/autodl-tmp/lbvs/all_dude/aa2ar/docking_test/lig/"

    create_mapfld(folder=path, filename=Path(pdbqt_path).stem, spacing=0.375, npts=[60,60,62], center=[72.899, 22.98, 31.04])
