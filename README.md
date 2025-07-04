[![DOI:10.1101/2024.01.22.576714](http://img.shields.io/badge/DOI-10.1101/2024.01.22.576714-B31B1B.svg?style=flat-square)](https://doi.org/10.1101/2024.01.22.576714)&emsp;[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/meyer-lab-cshl/BATMAN/blob/main/run_batman/pyBATMAN_Tutorial.ipynb) &emsp;[![PyPI version](https://badge.fury.io/py/pybatman.svg?icon=si%3Apython&style=flat-square)](https://badge.fury.io/py/pybatman)&emsp;[![Pepy Total Downloads](https://img.shields.io/pepy/dt/pybatman?style=flat-square&color=A52A2A)](https://pypistats.org/packages/pybatman)


# BATMAN: <ins>B</ins>ayesian Inference of <ins>A</ins>ctivation of <ins>T</ins>CR by <ins>M</ins>utant <ins>An</ins>tigens

A single T Cell Receptor (TCR) can recognize a diverse variety of peptides, an essential property known as TCR cross-reactivity.
Predicting which peptides a TCR cross-reacts to is critical for numerous applications, including predicting viral escape,
cancer neoantigen immunogenicity, autoimmunity, and off-target toxicity of T-cell-based therapies. But predicting TCR activation
is challenging due to the lack of both unbiased benchmarking datasets and computational methods that are sensitive to small
mutations to an epitope.

To address this, we curated the BATCAVE database, encompassing complete single amino acid mutational assays of over 22,000
TCR-peptide pairs, centered around 25 immunogenic epitopes, across both major histocompatibility complex classes, against 151 human and mouse TCRs.

We developed BATMAN, an interpretable Bayesian model, trained on BATCAVE, for predicting the peptides that activate
a TCR, and an active learning extension, which can efficiently map targets of a novel TCR by selecting a few peptides to
assay. We show that BATMAN outperforms existing methods, reveals structural and biochemical predictors of TCR-peptide
interactions, and can predict polyclonal T cell responses and TCR targets with high sequence dissimilarity.

<div align='center'>
<img src="BATMAN_schematic.jpg"   alt="BATMAN schematic diagram showing that it integrates mutational
  scan datasets across many TCRs to build a hierarchical Bayesian inference model. BATMAN infers hyperparameters from the training database
  and uses them to generate prior distributions for cross-TCR AA distance and TCR-specific
positional weights, which are multiplied and used as a predictor of TCR activation by a given mutant."/>
</div>

BATMAN predicts TCR activation by mutant peptides based on their distances to the TCR's index peptide. The peptide-to-index
distance is a product of a learned positional weight profile vector, corresponding to effects of mutated residues at different
positions in the sequence, and a learned AA substitution distance from the index peptide amino acid to the mutant amino acid.

BATMAN can be trained in two modes: (1) within-TCR, where the train and test peptides are associated with the same TCR, and
BATMAN-inferred positional weight profiles are TCR-specific, and (2) leave-one-TCR-out, where peptides are tested for activation
of a TCR left out of the training data, and BATMAN-inferred positional weight profile is common across all TCRs. For more information,
refer to our preprint [![DOI:10.1101/2024.01.22.576714](http://img.shields.io/badge/DOI-10.1101/2024.01.22.576714-B31B1B.svg)](https://doi.org/10.1101/2024.01.22.576714) (now accepted in Cell Systems). 

# Interactive tutorial of pyBATMAN
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/meyer-lab-cshl/BATMAN/blob/main/run_batman/pyBATMAN_Tutorial.ipynb)

For an interactive tutorial on pyBATMAN usage, refer to our [jupyter notebook](https://github.com/meyer-lab-cshl/BATMAN/blob/main/run_batman/pyBATMAN_Tutorial.ipynb).

# Installing and running pyBATMAN
[![PyPI version](https://badge.fury.io/py/pybatman.svg)](https://pypi.org/project/pybatman/)

It is advisable to install and run the Python implementation of BATMAN ('pyBATMAN') in a Conda environment with Python v=3.11. For instruction on creating and activating Conda environments, please refer to the [Conda user guide](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#). For example, execute the following in your Anaconda prompt to create and activate a new environment.

```
conda create -n BATMAN-env python=3.11
conda activate BATMAN-env
```

In the newly created Conda environment named 'BATMAN-env' with Python v3.11 installed, run the following command to install BATMAN.

```
pip install pybatman
```
Once installed successfully, you should be able to run a [test script](https://github.com/meyer-lab-cshl/BATMAN/blob/main/run_batman/test_script.py) on [sample input data](https://github.com/meyer-lab-cshl/BATMAN/blob/main/run_batman/test_input.csv), both available to download from this repository.

BATMAN combines TCR-pMHC data across multiple TCRs to infer TCR-specific positional weight profiles and an amino acid distance matrix. To train on an input dataset, pyBATMAN requires an amino acid matrix prior, and a training mode to be specified. 

Conveniently, pyBATMAN comes with more than 50 pre-saved conventional amino acid substitution distance matrices (e.g., hamming, blosum62, blosum100, pam30, pam100, dayhoff, gonnet etc), and we can use one of them. Alternatively, we can also pass a custom amino acid substitution matrix to pyBATMAN as a 20-by-20 Pandas dataframe with amino acid letters as row and column names. 

pyBATMAN runs in 3 training modes: (1) "full", which we will use below and is the most general one working without any constraint, (2) "symm", where we constrain the inferred amino acid matrix to be symmetric, and (3) "weight_only" where we only infer positional weights and directly use the amino acid distance matrix supplied without inferring amy amino acid matrix. pyBATMAN performs approximate Bayesian inference with automatic differentiation variational inference (ADVI) method to learn the weight profiles and amino acid matrices. You can optionally specify a random seed for reproducibility, and the number of steps in this method (default=20,000) and check that the loss function converges. An application of the "train" function of pyBATMAN looks like the following.

```
from pybatman.functions import train
inferred_weight_profiles, inferred_AA_matrix = train('test_input.csv','full',
                                             'blosum62',
                                             steps=50000, seed=10)
```

The input TCR-pMHC dataset must be a csv file containing 4 columns with the following names and structure: "tcr" and "peptide" columns referring to the name of the TCR and sequence of the peptide respectively in a TCR-pMHC pair. The "index" column refers to the index peptide of the TCR. Finally, the "activation" column refers to ordered categorical activation level of the TCR, when it interacts with the indicated peptide: 0 means no activation, 1 means weak activation, and 2 means strong activation. Note that whereas we standardized BATMAN to work with 2 or 3 TCR activation levels, you can have an arbitrary number of activation levels in your data, but they must start from 0 (non or weakest activation) and be consecutive integers denoting increasing activation without any missing level. Refer to the [sample input data](https://github.com/meyer-lab-cshl/BATMAN/blob/main/run_batman/test_input.csv) as an example. 

Once pyBATMAN outputs the positional weight profiles and the amino acid substitution matrix, we can multiply them to calculate a distance between an index peptide and its mutant, and predict if the mutant peptide activates the TCR or not based on how far it is from the index peptide. The distance between peptides can be conveniently calculated with pyBATMAN using the following function.

```
from pybatman.functions import peptide2index
peptide_distance = peptide2index(index_peptide,
                                 mutant_peptide_list,
                                 inferred_AA_matrix,
                                 inferred_weight_profiles)
```

For an interactive tutorial on different functions available with pyBATMAN, please refer to our [jupyter notebook](https://github.com/meyer-lab-cshl/BATMAN/blob/main/run_batman/pyBATMAN_Tutorial.ipynb). The Jupyter notebook trains and validates pyBATMAN on the test data and visualizes the results.

# Downloading BATMAN dataset
The fully curated database of publicly-available TCR-pMHC interactions can be downloaded from the [database folder](https://github.com/meyer-lab-cshl/BATMAN-paper/tree/main/results_batman/tcr_epitope_datasets/mutational_scan_datasets/database) in the publication repository.

# Citation
Our BATMAN [preprint](https://www.biorxiv.org/content/10.1101/2024.01.22.576714v3) is now accepted at Cell Systems. If you use or refer to BATMAN or BATCAVE in your work, please cite us as

```
@article{banerjee2025comprehensive,
  title={Comprehensive epitope mutational scan database enables accurate T cell receptor cross-reactivity prediction},
  author={Banerjee, Amitava and Pattinson, David J and Wincek, Cornelia L and Bunk, Paul and Axhemi, Armend and Chapin, Sarah R and Navlakha, Saket and Meyer, Hannah V},
  journal={bioRxiv},
  pages={2024--01},
  year={2025}
}
```
