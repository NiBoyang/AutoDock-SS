import os

def extract_binding_energy(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if "Estimated Free Energy of Binding" in line:
                return float(line.split('=')[1].split()[0])
    return None

def main():
    base_path = '' #This is the lib file path after running AutoDock-SS
    base_value = None
    scores = {}

    for folder in os.listdir(base_path):
        folder_path = os.path.join(base_path, folder)
        if os.path.isdir(folder_path):
            for file in os.listdir(folder_path):
                if file.endswith('_best_position.pdbqt'):
                    file_path = os.path.join(folder_path, file)
                    binding_energy = extract_binding_energy(file_path)

                    if binding_energy is not None:
                        if folder.endswith('_lig'):
                            base_value += binding_energy
                        else:
                            scores[folder] = binding_energy

    if base_value is not None:
        sorted_scores = sorted([(molecule, energy / base_value) for molecule, energy in scores.items()], key=lambda x: x[1], reverse=True)

        # When running AutoDock-SS, please put your reference ligand(s) in the whole library, 
        # and name it with the suffix '_lig'. For example, if your original reference ligand is 'abc', then name it as 'abc_lig'.
        with open(f"{base_path.replace('lib','')}/scores.txt", 'w') as file: 
            for molecule, score in sorted_scores:
                print(f"Molecule: {molecule}, Score: {score}")
                file.write(f"{score}\n")
    else:
        print("Base value not found.")

if __name__ == "__main__":
    main()
