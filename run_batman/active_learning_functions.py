import os
# Get the path of the current module
current_module_path = os.path.dirname(__file__)
# Define the path to your data directory
data_directory = os.path.join(current_module_path, 'data')

# Access your data files using relative paths
#dataset1_path = os.path.join(data_directory, 'dataset1.csv')
#dataset2_path = os.path.join(data_directory, 'dataset2.json')
# You can now read or manipulate your data

##############################################################
'''
Function to generate features of a list of mutant peptides, using their
position-dependent distances from the index peptide,
based on a given AA distance matrix.

Index peptide can be a string (common for all mutants) or a list (potentially
different for all mutants) with length equal to that of mutant list. All input 
sequences must have equal length.

AA matrix can be the name of a stored matrix (e.g., "BLOSUM100"), or be custom 
defined. Custom AA matrix must be a 20-by-20 Pandas Dataframe object 
with AAs as row and column names.

'''

def generate_mutant_features(index_peptide, # Single index peptide or list 
                 mutant_peptide_list, #List of mutant peptide sequences
                 aa_matrix): #Named or user-defined AA matrix   
    
    import pandas as pd
    import numpy as np
    import sys
    import os
    
    #AA names
    amino_acid_list=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R',
                     'S','T','V','W','Y']
    
    # Check if mutant peptide list is a list
    if (isinstance(mutant_peptide_list,list) == False):
        sys.exit("Mutant peptides must be supplied as a list")
        
    # Check if there is only one index or as many as mutants
    if ((isinstance(index_peptide, str)==False) and 
        (len(mutant_peptide_list)!=len(index_peptide))):
        sys.exit("Supply either a unique index peptide or one for each mutant")
     
    # Check if index peptides form a list
    if ((isinstance(index_peptide, str)==False) and 
        ((isinstance(index_peptide,list) == False))):
        sys.exit("More than one index peptide must be supplied as a list")
        
    
    # Check if all mutants and index peptides have the same length
    if isinstance(index_peptide, str): #One index peptide provided
        if len(np.unique(np.char.str_len([index_peptide]+mutant_peptide_list)))!=1:
            sys.exit("All index and mutant sequences must have equal length")
    else: #Index peptide list provided
        if len(np.unique(np.char.str_len(index_peptide+mutant_peptide_list)))!=1:
            sys.exit("All index and mutant sequences must have equal length")
    
    
    
    # Check if input sequences contain non-AA and wildcard characters
    if isinstance(index_peptide, str): #One index peptide provided
        if set(list("".join(([index_peptide] +
                             mutant_peptide_list)))).issubset(
                                 set(amino_acid_list)) == False:
           sys.exit("Discard sequences with non-AA characters, including X and *")
    else: #Index peptide list provided
        if set(list("".join((index_peptide +
                             mutant_peptide_list)))).issubset(
                                 set(amino_acid_list)) == False:
           sys.exit("Discard sequences with non-AA characters, including X and *")        
    
    
    
    # If AA matrix provided, check size is 20-by-20 and it has row and column names
    if isinstance(aa_matrix, str)==False: #Name of AA matrix not provided
        if ((isinstance(aa_matrix, pd.DataFrame) == False) or
            (aa_matrix.shape != (20,20)) or
            (set(aa_matrix.index) != set(amino_acid_list)) or
            (set(aa_matrix.columns) != set(amino_acid_list))):
            sys.exit("".join(["Custom AA matrix must be a 20-by-20 Pandas ",
                        "Dataframe object with AAs as row and column names."]))
        else: #If all goes well, load custom AA matrix
            aa_distance_matrix = aa_matrix
    
    # Load named AA distance matrix if the name is provided
    if isinstance(aa_matrix, str): #Name of AA matrix provided
        # Check if input AA distance matrix name exists, if it is provided
        filename = "".join([aa_matrix,'.csv'])
        if os.path.isfile(os.path.join(data_directory,'AA_matrices',filename))==False:
            sys.exit("".join(['AA matrix ', aa_matrix, ' does not exist. ',
                              'Check spelling and use all lower cases. ',
                       'You can also define your custom matrix as a 20-by-20 ',
                        'pandas DataFrame with AAs as row and column names.']))
        else: #load stored AA matrix
            aa_distance_matrix = pd.read_csv(os.path.join(data_directory,
                                                          'AA_matrices',
                                                          filename),
                                             index_col=0);    
    
    
    # save row and column AA orders for row to column AA substitutions
    from_AA = np.array(aa_distance_matrix.index).astype(str)
    to_AA = np.array(aa_distance_matrix.columns).astype(str)
    
    # Convert to np matrix for easy indexing
    aa_distance_matrix = aa_distance_matrix.to_numpy()
    
    #infer peptide length
    if isinstance(index_peptide, str): #One index peptide provided
        peptide_length = len(index_peptide)
    else: #List of index peptides provided
        peptide_length = len(index_peptide[0]) 
    
    # To locate mutation positions together in all mutant sequences, join them all
    if isinstance(index_peptide, str): #One index peptide provided
        index_peptide_repeated = np.array(list("".join(list(np.repeat(index_peptide,
                                            len(mutant_peptide_list))))))
    else: #List of index peptide provided
        index_peptide_repeated = np.array(list("".join(index_peptide)))
        
    mutant_peptides_joined = np.array(list("".join(mutant_peptide_list)))
    
    # Compare mutant seqs to index peptide seq to find mutation locations
    is_mutation_location = np.compare_chararrays(list(index_peptide_repeated),
                            list(mutant_peptides_joined),
                            "!=",'True')
    
    # Locations of WT AAs in from_AA array
    from_AA_index = np.nonzero(index_peptide_repeated[is_mutation_location,None] ==
                               from_AA)[1]
    # Locations of mutated AAs in to_AA array
    to_AA_index = np.nonzero(mutant_peptides_joined[is_mutation_location,None] ==
                               to_AA)[1]
    
    # Collect distance matrix elements corresponding to mismatches
    aa_distance = np.zeros(len(index_peptide_repeated)) #initialize
    aa_distance[is_mutation_location] = aa_distance_matrix[from_AA_index,
                                                           to_AA_index]
    
    # Reshape AA distance array to dim #mutants-by-peptide_length
    aa_distance = np.reshape(aa_distance, 
                             (len(mutant_peptide_list),peptide_length))
    
    
    return aa_distance

##########################################################################

'''
Function to generate peptide to index distance by multiplying positional weight
profile with features of a list of mutant peptides, created using their
position-dependent distances from the index peptide,
based on a given AA distance matrix.

Index peptide can be a string (common for all mutants) or a list (potentially
different for all mutants) with length equal to that of mutant list. All input 
sequences must have equal length.

AA matrix can be the name of a stored matrix (e.g., "BLOSUM100"), or be custom 
defined. Custom AA matrix must be a 20-by-20 Pandas Dataframe object 
with AAs as row and column names.

Weight profile should be array of size peptide_length-by-1 (common for all) or 
peptide_length-by-#mutants (different weights for different mutants)

'''
def peptide2index(index_peptide, # Single index peptide or list 
                 mutant_peptide_list, #List of mutant peptide sequences
                 aa_matrix, #Named or user-defined AA matrix
                 weight_profile): #array of weights
    
    import numpy as np
    import sys
    
    
    # Create mutant features
    mutant_features = generate_mutant_features(index_peptide, 
                                               mutant_peptide_list, 
                                               aa_matrix)
    
    # Check if weight profile is an array with desired size
    peptide_length = len(mutant_peptide_list[0])    
    if ((weight_profile.dtype not in ('int','float')) or
        (np.shape(weight_profile) not in ((1,peptide_length),
                                (len(mutant_peptide_list),peptide_length)))):
        sys.exit("".join(["Weight profile must be a numerical array of shape ",
                  "(1,peptide length) or (number of mutants,peptide length)"]))
    
    # Calculate peptide-to-index distances by multiplication
    if (np.shape(weight_profile)==(1,peptide_length)):
        distance = mutant_features.dot(weight_profile.transpose())[:,0]
    else:
        distance = np.multiply(mutant_features,weight_profile).sum(axis=1)
    
    return distance
#############################################################################
'''
Function that outputs sampled positional weights after pooled inference with peptide
features and TCR activation categories
'''


def pooled_inference_weights_only(TCR_index, # TCR index for each peptide
                                  peptide_feature, # features for each peptide
                                  peptide_activation_category,#activation category
                                  seed, steps): #Random seed and #steps for sampler
    import arviz as az
    import pymc as pm
    import numpy as np
    
    # Infer number of TCRs, peptide length, and number of TCR activation levels
    n_tcr = len(np.unique(TCR_index))
    peptide_length = np.shape(peptide_feature)[1]
    n_level = 1 + max(peptide_activation_category)
    
    # Build Bayesian classifier
    with pm.Model() as peptide_classifier_model:        
        
        # TCR-specific parameters
        
        # positional weights, pooled over TCRs and positions
        weights = pm.Beta("weights",
                          alpha = pm.Gamma("alpha_w", mu=1.0, sigma=0.5),
                          beta = pm.Gamma("beta_w", mu=5.0, sigma=1.0),
                          shape=(n_tcr,peptide_length))        
        
        # Hyperprior parameters for TCR-specific intercept        
        mu_bar = pm.Normal("mu_bar", mu=0, sigma=2)
        sigma_bar = pm.HalfNormal("sigma_bar", sigma=2)
        normal = pm.Normal("normal", mu=0, sigma=1, shape=n_tcr) 
        
        intercepts=pm.Deterministic("intercepts",mu_bar+sigma_bar*normal)  
             
        
        # Full Predictor
        eta = intercepts[TCR_index] - pm.math.sum(
            peptide_feature*weights[TCR_index,:],axis=1)
        
        # Binomial Regression
        # Generate cutpoints
        cutpoints=pm.Normal("cutpoints", 
                            mu=0.0, sigma=2.0, shape=[n_level-1],
         transform=pm.distributions.transforms.univariate_ordered,
         initval=np.linspace(0.1,0.5,n_level-1))
        
        peptide_activation_category_obs = pm.OrderedLogistic(
                                        "peptide_activation_category_obs",
                                        eta=eta,cutpoints=cutpoints,
                                        observed=peptide_activation_category,
                                        compute_p=False)
        
    # Sampling with approximate posterior
    with peptide_classifier_model:
         posterior_draws=pm.fit(n=steps,method="advi",
                                random_seed=seed,progressbar=True)
         inferred_params = az.summary(posterior_draws.sample(50000))
         
         
    # Extract position-dependent weights of TCRs         
    inferred_weights=np.reshape(inferred_params.iloc[
        (n_tcr+3):(n_tcr+3+n_tcr*peptide_length),0].to_numpy(),
        newshape=(n_tcr,peptide_length),order='C') #TCR-by-position
    
    # Extract sd of position-dependent weights of TCRs         
    inferred_weights_sd=np.reshape(inferred_params.iloc[
        (n_tcr+3):(n_tcr+3+n_tcr*peptide_length),1].to_numpy(),
        newshape=(n_tcr,peptide_length),order='C') #TCR-by-position
    
    return inferred_weights,inferred_weights_sd

#######################################################################

'''
Function that outputs sampled positional weights after pooled inference with peptide
features, specified priors and TCR activation categories
'''


def inference_weights_only(peptide_feature, # features for each peptide
                                  peptide_activation_category,#activation category
                                  prior_mean_sd,#beta prior mean and sd
                                  seed, steps): #Random seed and #steps for sampler
    import arviz as az
    import pymc as pm
    import numpy as np
    import sys
    
    # Infer number of TCRs, peptide length, and number of TCR activation levels
    peptide_length = np.shape(peptide_feature)[1]
    n_level = 1 + max(peptide_activation_category)
    
    # Raise shape error if prior size does not match peptide length
    if ((prior_mean_sd.dtype not in ('int','float')) or
        (np.shape(prior_mean_sd)!=(2,peptide_length))):
        sys.exit("".join(["prior array must be a numerical array of shape ",
                  "(2,peptide length)"]))
    
    # Build Bayesian classifier
    with pm.Model() as peptide_classifier_model:        
        
        # TCR-specific parameters
        
        # positional weights, pooled over TCRs and positions
        weights = pm.Normal("weights",
                          mu = prior_mean_sd[0,:],
                          sigma = pm.Exponential("sd",lam=1,
                                                 shape=(1,peptide_length))) 
        
       
        # Hyperprior parameters for TCR-specific intercept        
        mu_bar = pm.Normal("mu_bar", mu=0, sigma=2)
        sigma_bar = pm.HalfNormal("sigma_bar", sigma=2)
        normal = pm.Normal("normal", mu=0, sigma=1) 
        
        intercepts=pm.Deterministic("intercepts",mu_bar+sigma_bar*normal)  
             
        
        # Full Predictor
        eta = intercepts - pm.math.sum(
            peptide_feature*weights,axis=1)
        
        # Binomial Regression
        # Generate cutpoints
        cutpoints=pm.Normal("cutpoints", 
                            mu=0.0, sigma=2.0, shape=[n_level-1],
         transform=pm.distributions.transforms.univariate_ordered,
         initval=np.linspace(0.1,0.5,n_level-1))
        
        peptide_activation_category_obs = pm.OrderedLogistic(
                                        "peptide_activation_category_obs",
                                        eta=eta,cutpoints=cutpoints,
                                        observed=peptide_activation_category,
                                        compute_p=False)
        
    # Sampling with approximate posterior
    with peptide_classifier_model:
         posterior_draws=pm.fit(n=steps,method="advi",
                                random_seed=seed,progressbar=True)
         inferred_params = az.summary(posterior_draws.sample(50000))
    #print(inferred_params)
    # Extract position-dependent weights of TCRs         
    inferred_weights=np.reshape(inferred_params.iloc[
        0:(0+peptide_length),0].to_numpy(),
        newshape=(1,peptide_length),order='C')
    
    # Extract sd of position-dependent weights of TCRs         
    inferred_weights_sd=np.reshape(inferred_params.iloc[
        0:(0+peptide_length),1].to_numpy(),
        newshape=(1,peptide_length),order='C')
        
    return inferred_weights,inferred_weights_sd

#######################################################################
###################################################################


'''
Function that outputs sampled positional weights and AA matrix after pooled
inference with peptide features and TCR activation categories
'''


def pooled_inference(TCR_index, # TCR index for each peptide
                                  peptide_feature, # features for each peptide
                                  peptide_activation_category,#activation category
                                  aa_change_index,#indexing which AA to which
                                  mode, # symmetric (default) or full
                                  seed, steps): #Random seed and #steps for sampler
    import arviz as az
    import pymc as pm
    import numpy as np

    # Infer number of TCRs, peptide length, and number of TCR activation levels
    n_tcr = len(np.unique(TCR_index))
    peptide_length = np.shape(peptide_feature)[1]
    n_level = 1 + max(peptide_activation_category)
    
    # Renumber AA change index matrix with unique AA change indices
    unique_indices,_, aa_change_index_unique = np.unique(aa_change_index, 
                                   return_index=True,
                                   return_inverse=True);
    aa_change_index_unique = aa_change_index_unique.reshape(
                                                     np.shape(peptide_feature))
    
    
    # Build Bayesian classifier
    with pm.Model() as peptide_classifier_model: 
        
        # TCR-indepedent common amino acid distance matrix multiplier flattened        
        aa_distance_multiplier = pm.math.concatenate([[0],
                                          pm.Normal("aa_distance_multiplier",
                                mu=pm.Normal("mu", mu=0, sigma=0.5),
                                sigma=pm.Exponential("sigma",lam=1),
                                shape=len(unique_indices)-1)], 
                                          axis=0)
        #0 at beginning put for No AA substituion
        
        # TCR-specific parameters
        
        # positional weights, pooled over TCRs and positions
        weights = pm.Beta("weights",
                          alpha = pm.Gamma("alpha_w", mu=1.0, sigma=0.5),
                          beta = pm.Gamma("beta_w", mu=5.0, sigma=1.0),
                          shape=(n_tcr,peptide_length))        
        
        # Hyperprior parameters for TCR-specific intercept        
        mu_bar = pm.Normal("mu_bar", mu=0, sigma=2)
        sigma_bar = pm.HalfNormal("sigma_bar", sigma=2)
        normal = pm.Normal("normal", mu=0, sigma=1, shape=n_tcr) 
        
        intercepts=pm.Deterministic("intercepts",mu_bar+sigma_bar*normal)  
             
        
        # Full Predictor
        eta = intercepts[TCR_index] - pm.math.sum(
            weights[TCR_index,:]*peptide_feature*
            (1 + aa_distance_multiplier[aa_change_index_unique]),
            axis=1) #positional weight * D *(1+multiplier) summed over positions
                    
        
        # Binomial Regression
        # Generate cutpoints
        cutpoints=pm.Normal("cutpoints", 
                            mu=0.0, sigma=2.0, shape=[n_level-1],
         transform=pm.distributions.transforms.univariate_ordered,
         initval=np.linspace(0.1,0.5,n_level-1))
        
        peptide_activation_category_obs = pm.OrderedLogistic(
                                        "peptide_activation_category_obs",
                                        eta=eta,cutpoints=cutpoints,
                                        observed=peptide_activation_category,
                                        compute_p=False)
        
    # Sampling with approximate posterior
    with peptide_classifier_model:
         posterior_draws=pm.fit(n=steps,method="advi",
                                random_seed=seed,progressbar=True)
         inferred_params = az.summary(posterior_draws.sample(50000))     
         
    # Extract position-dependent weights of TCRs         
    inferred_weights=np.reshape(inferred_params.iloc[
        (n_tcr+len(unique_indices)+4):
        (n_tcr+len(unique_indices)+4+n_tcr*peptide_length),0].to_numpy(),
        newshape=(n_tcr,peptide_length),order='C') #TCR-by-position
        
    # Extract errors in position-dependent weights of TCRs
    inferred_weights_sd=np.reshape(inferred_params.iloc[
        (n_tcr+len(unique_indices)+4):
        (n_tcr+len(unique_indices)+4+n_tcr*peptide_length),1].to_numpy(),
        newshape=(n_tcr,peptide_length),order='C') #TCR-by-position
        
    # Extract inferred AA distance matrix multipliers 
    aa_multiplier = inferred_params.iloc[1:len(unique_indices),0].to_numpy()
    aa_multiplier = np.insert(aa_multiplier,0,0)
    #Insert 0 for no substitution
    
    # Reconstruct full AA factor matrix (to be multiplied with regularizer)
    # Initialize with inferred mean of the pooling distribution for AA multiplier
    inferred_aa_matrix = np.zeros((20,20))+(1+inferred_params.iloc[0,0])
    
    if mode=='symm': #Construct symmetric AA matrix multiplier
        for aa1 in np.arange(0,20,1):
            for aa2 in np.arange(0,20,1):
               if 20*min(aa1,aa2)+max(aa1,aa2) in unique_indices:
                   inferred_aa_matrix[aa1,aa2] = 1+aa_multiplier[
                       np.where(unique_indices==20*min(aa1,aa2)+max(aa1,aa2))[0][0]]                   
                   # Takes care of AA reindexing done at the beginning
      
    if mode=='full': #Construct full AA matrix multiplier
        for aa1 in np.arange(0,20,1):
            for aa2 in np.arange(0,20,1):
               if 20*aa1+aa2 in unique_indices:
                   inferred_aa_matrix[aa1,aa2] = 1+aa_multiplier[
                       np.where(unique_indices==20*aa1+aa2)[0][0]]                   
                   # Takes care of AA reindexing done at the beginning
 
    return inferred_weights,inferred_weights_sd,inferred_aa_matrix
##############################################################################

'''
Function to train BATMAN, using TCR activation data provided in a csv file.
Refer to the example file to see format of the input TCR activation datafile.

TCR activation levels must start from 0 (weakest) and are integers 
(no missing level allowed)

AA matrix can be the name of a stored matrix (e.g., "BLOSUM100"), or be custom 
defined. Custom AA matrix must be a 20-by-20 Pandas Dataframe object 
with AAs as row and column names.

Runs in one of 3 modes: weights_only, full, or symm, based on if only weights
are inferred (AA matrix being the one specified) or if both weights and AA matrices
are inferred (symmetric or full AA matrix)

Also takes #steps for sampling and a seed argument for reproduction 

'''

def train(filename,#Path to file containing TCR data (see example for format)
          mode,# weights_only, symm, or full
          # Default values of optional parameters
          aa_matrix='BLOSUM100',#Named or user-defined AA matrix used for regularization
          seed=100,#seed for sampling
          steps=20000):# number of steps for sampling


    import pandas as pd
    import numpy as np
    import os.path, sys
    
    
    # Mode must be specified as one of 3 options
    if mode not in ('weights_only','symm','full'):
        sys.exit("Mode must be one of: 'weights_only','symm', 'full'")
    
    # Check if file exists and is csv
    if os.path.isfile(filename)==False:
        sys.exit("File does not exist. Check filename and/or directory again.")
    
    if filename.endswith('.csv')==False:
        sys.exit("File must be of csv format. See example input for details.")
    
    # Check that file has correct headers for relevant columns
    if {'activation', 'index', 'peptide', 'tcr'}.issubset(
            set(pd.read_csv(filename).columns))==False:
        sys.exit(''.join(["Input file must have these 4 headers:", 
                 "'activation', 'index', 'peptide', 'tcr'.",
                 "Check spelling and have all headers in lower case."]))

    # Check that 'activation' column has 0 and integer data
    if ((pd.read_csv(filename)['activation'].to_numpy().dtype.kind not in ('i','u'))
    or (min(pd.read_csv(filename)['activation'].to_numpy()))!=0):
        sys.exit("activation levels must be 0 and positive integer(s) (smaller=weaker)")
    
    # Check that all activation levels between 0 to K (max) are present in data
    if sum(np.unique(pd.read_csv(filename)['activation'])!=
           np.arange(0,1+max((pd.read_csv(filename)['activation'].to_numpy())),1))!=0:
        sys.exit("One or more missing activation levels in data")
    
    
    # Read and featurize peptide data
    peptide_data = pd.read_csv(filename)
    
    # Assign indices to unique TCRs in data
    TCR_names, TCR_index = np.unique(peptide_data['tcr'], return_inverse=True)
    
    # Make list of index peptide
    index_list = peptide_data['index'].tolist()
    # Make list of mutant peptide
    peptide_list = peptide_data['peptide'].tolist()
    
    # Make peptide features
    peptide_feature = generate_mutant_features(index_list,
                                               peptide_list,
                                               aa_matrix)
    
    # Peptide activation categories
    peptide_activation_category = peptide_data['activation'].to_numpy()
    
    if mode=='weights_only':
        # Run sampling for Inference of only weights
        inferred_weights,inferred_weights_sd = pooled_inference_weights_only(TCR_index, 
                                                         peptide_feature, 
                                                         peptide_activation_category, 
                                                         seed=seed, steps=steps)
        # Add TCR names
        inferred_weights = pd.DataFrame(inferred_weights,index=TCR_names)
        inferred_weights_sd = pd.DataFrame(inferred_weights_sd,index=TCR_names)
       
        return inferred_weights,inferred_weights_sd  

    
    
    else: 
        
        #AA name list
        aa_list=np.array(['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R',
                         'S','T','V','W','Y'])
        # List which AA substitutions are in data (based on full or symm AA matrix)
        # Indices of index and mutant AAs in AA list (flattened)
        index_aa = np.where(np.array(list(''.join(index_list)))[:, None] == 
                            aa_list[None, :])[1]
        peptide_aa = np.where(np.array(list(''.join(peptide_list)))[:, None] == 
                              aa_list[None, :])[1]
        
        if mode=='symm':
            #Symmetric case, indicate a single number to index AA substitution, 
            #and reshape flattened array
            aa_subs = (20*np.minimum(index_aa,peptide_aa)+ np.maximum(index_aa,peptide_aa))
            aa_subs[index_aa==peptide_aa]=0 #assign 0 to positions with unchanged AA 
            aa_change_index = aa_subs.reshape(np.shape(peptide_feature))
            
        if mode=='full':
            #full AA matrix, indicate a single number to index AA substitution, 
            #and reshape flattened array
            aa_subs = (20*index_aa + peptide_aa)
            aa_subs[index_aa==peptide_aa]=0 #assign 0 to positions with unchanged AA 
            aa_change_index = aa_subs.reshape(np.shape(peptide_feature))   
        
        # Run Inference
        inferred_weights,inferred_weights_sd,inferred_aa_matrix = pooled_inference(TCR_index,
                                                               peptide_feature, 
                                                               peptide_activation_category, 
                                                               aa_change_index, 
                                                               mode=mode, 
                                                               seed=seed, steps=steps)
        # Add TCR names
        inferred_weights = pd.DataFrame(inferred_weights,index=TCR_names)
        inferred_weights_sd = pd.DataFrame(inferred_weights_sd,index=TCR_names)

        
        # multiply to regularizer AA matrices
        # Load named AA distance matrix data if the name is provided
        if isinstance(aa_matrix, str): #Name of AA matrix provided
            filename = "".join([aa_matrix,'.csv'])
            aa_matrix_prior = pd.read_csv(os.path.join(data_directory,
                                                          'AA_matrices',
                                                          filename),
                                          index_col=0)
        else: #full AA matrix provided
            aa_matrix_prior = aa_matrix
          
        inferred_aa_matrix = aa_matrix_prior*inferred_aa_matrix 
        
    return inferred_weights,inferred_weights_sd,inferred_aa_matrix

##############################################################################
'''
Function for one active learning step.
Takes a train and test dataframe, priors, #steps, and outputs auc and weights
'''
def active_learning_cycle(active_train_data,
          active_test_data,
          aa_matrix,
          prior_mean_sd,
          steps=20000,
          seed=101):
    
    import numpy as np
    import pandas as pd
    import pymc as pm
    import arviz as az
    from sklearn.metrics import roc_curve, auc

    peptide_feature = generate_mutant_features(active_train_data['index'].tolist(),
                                               active_train_data['peptide'].tolist(),
                                               aa_matrix)
    peptide_activation_category = active_train_data.loc[:,'activation']
    
    # Train BATMAN to infer TCR-independent weights and AA distance matrix
    inferred_weights_al, _ = inference_weights_only(peptide_feature,
                                             peptide_activation_category,
                                             prior_mean_sd,
                                             steps = steps,
                                             seed = seed)
    
    # Normalize learned weights
    inferred_weights_al = (inferred_weights_al-np.min(inferred_weights_al))/(
        np.max(inferred_weights_al)-np.min(inferred_weights_al))        
    
    # Get mutant distance
    mutant_distance_al = peptide2index(active_test_data['index'].tolist(), 
                              active_test_data['peptide'].tolist(),
                              aa_matrix, 
                              inferred_weights_al)
    
    # Calculate pairwise AUC for AL-inferred distances
    fpr = dict()
    tpr = dict()
    auc_al = np.zeros((1,3))
    for i in np.arange(0,3):
        fpr[i], tpr[i], _ = roc_curve(active_test_data.activation[active_test_data.activation !=i],
                                      -mutant_distance_al[active_test_data.activation !=i],
                                      pos_label=max(set(np.arange(0,3)).difference({i})))
        auc_al[0,i] = auc(fpr[i], tpr[i])
        
    auc_mean = np.nanmean(auc_al) 
    
    return inferred_weights_al,auc_mean

