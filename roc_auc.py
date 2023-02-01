from importlib.resources import path
import os

def roc_auc_ef(ranked_file, actives_file, candidate_name):

    from sklearn import metrics
    import matplotlib.pyplot as plt
    import numpy as np

    with open(actives_file, "r") as tempfile:
        actives = [line.strip("\n") for line in tempfile.readlines()]

    with open(ranked_file, "r") as tempfile:
        total = [line.strip("\n") for line in tempfile.readlines()]

    num_of_actives = len(actives)
    num_of_decoys = len(total) - num_of_actives

    act, dcy = 0, 0
    fpr, tpr = [], []
    for candidate in total:
        if candidate in actives:
            act += 1
        else:
            dcy += 1
        fpr.append(dcy/num_of_decoys)
        tpr.append(act/num_of_actives)

    with plt.style.context(['science', "high-vis", "grid"]):
        pparam = dict(xlabel='False Positive Rate', ylabel='True Positive Rate')
        fig, ax = plt.subplots()
        ax.set_title(f"ROC of {candidate_name}")
        ax.plot(fpr, tpr, label="NCF-LBVS")
        ax.plot(fpr, fpr, label="Random")
        ax.legend(title=f"AUC = {metrics.auc(fpr, tpr):.5f}")
        ax.autoscale(tight=True)
        ax.set(**pparam)
        ax.set_ylim(bottom=-0.03, top=1.03)
        ax.set_xlim(left=-0.03, right=1.03)
        fig.savefig(f"ROC_{candidate_name}.svg")

    data_set_count = len(total)
    active_ratio = np.divide(num_of_actives, data_set_count)
    Subsets = [[int(round(data_set_count * 0.01)), 1 ],
               [int(round(data_set_count * 0.05)), 5 ],
               [int(round(data_set_count * 0.10)), 10],
               [int(round(data_set_count * 0.15)), 15],
               [int(round(data_set_count * 0.20)), 20],
               [int(round(data_set_count * 0.25)), 25],
               [int(round(data_set_count * 0.50)), 50],
               [int(round(data_set_count * 0.99)),100]]

    EF_active_found  = 0.0
    active_found     = 0.0
    active_not_found = 0.0
    Active_Found     = []

    with open(f"EF_{candidate_name}.txt", "w") as tempfile:
        tempfile.write(f"""\
Number of total ligands: {data_set_count}\n
Number of known actives: {num_of_actives}\n
Number of known decoys: {num_of_decoys}\n
AUC : {metrics.auc(fpr, tpr):.5f}\n
""")
        for data_set_index, compound in enumerate(total):
            Name = compound.split()

            if data_set_index <= num_of_actives and Name[0] in actives:
                active_found += 1.0
            elif data_set_index > num_of_actives and Name[0] in actives:
                active_not_found += 1.0

            ## For enrichment factor calculation
            if Name[0] in actives:
                EF_active_found += 1.0
                Active_Found.append(Name[0])

        # Calculate the Enrichment factor at certain subset ratio
            for Subset in Subsets:
                if data_set_index == Subset[0]:
                    EF = float(EF_active_found / data_set_index) / active_ratio
                    tempfile.write(f"EF at {Subset[1]} % : {EF:.2f}\n")
                    # print('  EF at {0} %: {1:.2f}'.format(Subset[1], EF))
                    
if __name__ == "__main__":

    dude_path = "/root/autodl-tmp/lbvs/all_dude/"
    candidates = [paths for paths in os.listdir(dude_path)]
    output_paths = [f"{dude_path}{i}/docking_test/" for i in candidates]
    print(output_paths)

    for ind, candidate_path in enumerate(output_paths):
        os.chdir(candidate_path)
        roc_auc_ef("ranked.txt", "actives.txt", candidates[ind])
