# AutoDock-SS: Ligand-Based Virtual Screening Algorithm

## Introduction
AutoDock-SS is a novel algorithm designed for ligand-based virtual screening. It leverages advanced computational techniques to efficiently screen large compound libraries against a reference ligand, facilitating the discovery of potential drug candidates.

## Features
- Efficient screening of large compound libraries
- Utilization of `.pdbqt` format for reference ligands
- Support for `.sdf.gz` compressed SDF file format for compound libraries
- Integration with AutoDock-GPU for enhanced performance

## Installation
1. For AutoDock-GPU compilation, please refer to https://github.com/ccsb-scripps/AutoDock-GPU
2. Clone the AutoDock-SS repository:
```
git clone https://github.com/NiBoyang/AutoDock-SS.git
```
3. Create a conda environment using the provided `env.yaml` file:
```
conda env create -f env.yaml
```
4. Activate the environment:
```
conda activate adss
```

## Usage
1. Prepare your reference ligand in `.pdbqt` format and your compound library in a compressed `.sdf.gz` file.
2. Modify the common variables in `adss_main.py` to match your file paths and settings:
```
autodock_gpu_exec_path = "" # directory of Autodock-GPU exec file
lib_path = f'' # directory of your VS library
lig_path = f'' # directory of your reference ligand
path_of_scripts = "" # directory of Autodock-SS scripts
```
3. Run AutoDock-SS:
```
python adss_main.py
```

## Components
- [`adss_main.py`](https://github.com/NiBoyang/AutoDock-SS/blob/master/adss_main.py): Main script to run the virtual screening process. 
- [`env.yaml`](https://github.com/NiBoyang/AutoDock-SS/blob/master/env.yaml): Conda environment file with necessary dependencies.
- [`ligandprep.py`](https://github.com/NiBoyang/AutoDock-SS/blob/master/ligandprep.py): Script for ligand preparation. Called in the main file.
- [`roc_auc.py`](https://github.com/NiBoyang/AutoDock-SS/blob/master/roc_auc.py): Script for evaluating the performance of the screening. Called in the main file.
- Utility scripts: `geometry_utils.py`, `grid_map_utils.py`, `file_utils.py` for various computational tasks.
- [`DA.ipynb`](https://github.com/NiBoyang/AutoDock-SS/blob/master/util/DA.ipynb): A Jupyter notebook for data analysis.

## Example Results
The `Result_Examples` directory contains example outputs from the AutoDock-SS algorithm, showcasing its capabilities.

## Contributing
Contributions to AutoDock-SS are welcome. Please read the contribution guidelines before submitting pull requests.

## License
AutoDock-SS is released under [MIT License](https://opensource.org/licenses/MIT).

## Acknowledgements
Special thanks to all contributors and the open-source community for their support and contributions to this project.

