# Note that windows users might need to install the following to make things run faster
#conda install -c conda-forge m2w64-toolchain

# Load functions from pyBATMAN
from pybatman.functions import train, peptide2index

# Packages for downstream analysis of BATMAN outputs
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

'''
Train on test data and infer full asymmetric AA matrix and TCR-specific weights
'''
inferred_weights, inferred_AA_matrix = train('test_input.csv','full')

'''
Generate peptide2index distances for a selected TCR (TCR1)
'''
#Read data
peptide_data = pd.read_csv('test_input.csv')

# Extract data for a particular TCR
tcr_name = 'TCR1'

index_peptide = peptide_data[peptide_data.tcr==tcr_name]['index'].tolist()

mutant_peptide_list = peptide_data[peptide_data.tcr==
                                   tcr_name]['peptide'].tolist()

# Inferred positional weight profile for the selected TCR
weight_profile = inferred_weights.loc[[tcr_name]].to_numpy()

# generate peptide-to-index distances with pyBATMAN
peptide_distance = peptide2index(index_peptide,
                                 mutant_peptide_list,
                                 inferred_AA_matrix,
                                 weight_profile)

'''
plot histogram of peptides-to-index distances for 3 classes
'''
activation_level = peptide_data[peptide_data.tcr
                                ==tcr_name]['activation'].to_numpy()
bins = np.linspace(0, 6, 20)
plt.hist([peptide_distance[activation_level == 2],
          peptide_distance[activation_level == 1],
          peptide_distance[activation_level == 0]], 
         bins, label=['strong activation', 'weak activation', 'no activation'])
plt.xlabel('Peptide to index distance')
plt.ylabel('Count')
plt.legend(loc='upper right')
plt.rc('font', size=50)
plt.show()








