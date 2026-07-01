#!/Users/evolutioneco/anaconda3/bin/python
import importlib
import subprocess
import random
from hmmlearn import hmm
import pyreadr,random,time,h5py,gzip,sys,os
import allel; 
print('scikit-allel', allel.__version__)
import numpy as np
import pandas as pd
from hmmlearn import hmm
import seaborn as sns
#import bcolz
# %reload_ext memory_profiler
from itertools import compress,groupby
from functools import reduce
from sklearn.cluster import KMeans
from collections import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

# def load_pkgs():
#     import subprocess
#     import pyreadr,random,time,h5py,gzip,sys,os
#     import allel; 
#     print('scikit-allel', allel.__version__)
#     import numpy as np
#     import pandas as pd
#     from hmmlearn import hmm
#     random.seed(42)
#     import seaborn as sns
#     #import bcolz
#     # %reload_ext memory_profiler
#     from itertools import compress,groupby
#     from functools import reduce
#     from sklearn.cluster import KMeans
#     # from collections import *
#     import matplotlib as mpl
#     import matplotlib.pyplot as plt
#     import matplotlib.gridspec as gridspec
#     from matplotlib.lines import Line2D
#     from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#     from mpl_toolkits.axes_grid1 import make_axes_locatable
#     import scipy
#     from scipy.cluster.hierarchy import dendrogram, linkage
#     from scipy.spatial.distance import squareform
#     return None

######## function for panel E
def plot_pca_coords(coords, model, pc1, pc2, ax, sample_population):
    sns.despine(ax=ax, offset=5)
    x = coords[:, pc1]
    y = coords[:, pc2]
    for pop in populations:
        flt = (sample_population == pop)
        ax.plot(x[flt], y[flt], marker='o', linestyle=' ', color=pop_colours[pop], 
                label=pop, markersize=6, mec='k', mew=.5)
    ax.set_xlabel('PC%s (%.1f%%)' % (pc1+1, model.explained_variance_ratio_[pc1]*100))
    ax.set_ylabel('PC%s (%.1f%%)' % (pc2+1, model.explained_variance_ratio_[pc2]*100))
    

def fig_pca(coords, model, title, sample_population=None):
    if sample_population is None:
        sample_population = df_samples.population.values
    # plot coords for PCs 1 vs 2, 3 vs 4
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 2, 1)
    plot_pca_coords(coords, model, 0, 1, ax, sample_population)
    ax = fig.add_subplot(1, 2, 2)
    plot_pca_coords(coords, model, 2, 3, ax, sample_population)
    ax.legend(bbox_to_anchor=(1, 1), loc='upper left')
    fig.suptitle(title, y=1.02)
    fig.tight_layout()
    
def window_GT(callset, window, mode="SNP", start=0, step=None):
    '''
    split the genotype into windows based on the num of SNPs or physical chromosome position.
    default start from 0
    '''
    genotype = allel.GenotypeChunkedArray(callset['calldata/GT'])
    if mode == "SNP":
        while start < len(genotype) - window:
            start += window 
            gt_subset = genotype[start : start + window]
            yield gt_subset 
        else:
            gt_subset = genotype[start : ] #final should be 31497110
            yield gt_subset
            
    elif mode == "BP":
        start = 0
        end = start + 100
        med = (start + end) / 2
        pos_dict = defaultdict(list)
        for i in pos:
            #filter(lambda x: x not in subset_of_A, A)
            if end > i > start:
                pos_dict[med].append(i)
                #print("process 1 " + str(i))
            else:
                while i >= end:
                    start += 100
                    end  += 100
                    med = (start + end) / 2
            pos_dict[med].append(i) if i not in pos_dict[med] else pos_dict[med]
            
        # temp pos of pos array
        temp_pos = 0 
        for k,v in pos_dict.items():
            temp_pos_start = temp_pos
            temp_pos_end = temp_pos_start + len(v)
            temp_pos = len(v)
            gt_subset_BP = genotype[temp_pos_start : temp_pos_end]
            yield gt_subset_BP
            
    else:
        print("please specify the mode of sliding window: mode = ['SNP', 'BP']")
        

def vcf_file_handle(vcffile, data_path, POP, CHROM, START, END):
    """
    Reads the VCF file and returns the callset.
    """
    chrom = CHROM if "female_genome" in data_path else "LR8806{}.1".format(44 + int(CHROM.split("LG")[1]))
    callset = allel.read_vcf(vcffile, region="{}:{}-{}".format(chrom, START, END))
    return callset
    

def cal_pca_het_cluster(callset, START=0, END=None, n_components=10, scaler='patterson'):
    """
    Perform genotyping from a VCF file.

    Parameters:
    - callset: The allel callset object containing the VCF data
    - start: The start position for SNP filtering (default: 0)
    - end: The end position for SNP filtering (default: None, meaning no specific end)
    - n_components: The number of principal components to calculate (default: 10)
    - scaler: The scaling method for PCA (default: 'patterson')

    Returns:
    - coords: PCA coordinates
    - model: PCA model
    """
    # Extract genotype data
    g = allel.GenotypeArray(callset['calldata/GT'])
    
    # Count allele frequencies
    ac = g.count_alleles()[:]
    
    # SNP filtering
    flt = (ac.max_allele() == 1) & (ac[:, :2].min(axis=1) > 1)
    gf = g.compress(flt, axis=0)
    gn = gf.to_n_alt() # fill=-1 with missing data?
    # cal pca
    coords1, model1 = allel.pca(gn, n_components=n_components, scaler=scaler)
    #cal_het
    sample_list = []
    for i in range(len(callset['samples'])):
        het_count = g[:,i].count_het()
        prop_het_sample = het_count*100 / (END - START)
        sample_list.append([callset['samples'][i], coords1[i][0], coords1[i][1], prop_het_sample])
    return gn, gf, coords1, model1, sample_list

def cal_pca(gn):
        # Perform PCA
    coords1, model1 = allel.pca(gn, n_components=n_components, scaler=scaler)
    return coords, model

def cal_het(callset, coords1, START, END):
    sample_list = []
    g = allel.GenotypeArray(callset['calldata/GT'])
    for i in range(len(callset['samples'])):
        het_count = g[:,i].count_het()
        prop_het_sample = het_count*100 / (END - START)
        sample_list.append([callset['samples'][i], coords1[i][0], coords1[i][1], prop_het_sample])
    df = pd.DataFrame(sample_list, columns=["samples", "PC1", "PC2", "hetero"])
    return df

def clutering(df, n_clusters=3):
    #### clustering 
    X = df[['PC1']]
    kmeans = KMeans(n_clusters=n_clusters, n_init=10)
    kmeans.fit(X)
    df_cl = X.copy()
    # add one more column
    df_cl['cluster'] = kmeans.predict(X)
    df["cluster"] = kmeans.predict(X)
    df['sex'] = ["F" if "F" in i else "M" for i in df["samples"]]
    df['pop'] = [i.split('-')[0].upper() for i in df["samples"]]
    
    # sort by PC1 
    sort_df = df.sort_values(by=["PC1"])
    
    # Assuming you have the Counter object
    cnt = Counter(sort_df['cluster'])
    
    # Get the unique cluster values and their counts
    unique_clusters = list(cnt.keys())
    counts = list(cnt.values())
    
    # Map the original cluster values to new ones based on their frequency
    new_cluster_values = dict(zip(unique_clusters, range(len(unique_clusters))))
    
    # Create a new column 'new_cluster' based on the mapping
    sort_df['cluster_2'] = sort_df['cluster'].map(new_cluster_values)
    
    sort_df_sex = sort_df.sort_values(by=["cluster_2",  "pop", "PC1"])
    sort_df = sort_df_sex
    return sort_df


def run_vcftools(vcffile, POP, CHROM, START, END, sort_df):
    prefix = "{}/{}_{}_{}-{}Mb".format(os.path.dirname(vcffile), POP, CHROM, str(int(START)/1e6),str(int(END)/1e6))
    cluster0 = prefix + ".cluster0"
    sort_df[sort_df["cluster_2"]==0]["samples"].to_csv(
        cluster0, 
        index=False
    )
    
    cluster1 = prefix + ".cluster1"
    sort_df[sort_df["cluster_2"]==1]["samples"].to_csv(
        cluster1, 
        index=False
    )
    cluster2 = prefix + ".cluster2"
    sort_df[sort_df["cluster_2"]==2]["samples"].to_csv(
        cluster2, 
        index=False
    )
    
    ##### run linux cmd #########################
    cmd = "vcftools --gzvcf {} \
                    --weir-fst-pop {} \
                    --weir-fst-pop {} \
                    --fst-window-size 10000 --fst-window-step 5000 \
                    --out {} ".format(vcffile, cluster0, cluster2, prefix)
    
    # Execute the command and capture the output
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    # Check if the command ran successfully and print output or error message
    if result.returncode == 0:
        print("VCFtools ran successfully. Output saved to:", prefix)
    else:
        print("VCFtools encountered an error:", result.stderr)
    
    return None


def get_mds_20(rds_df, num_mds, flip=None):
    """
    Extracts and computes MDS differences for the specified DataFrame within rds_df.
    
    Parameters:
    - rds_df: Dictionary of DataFrames
    - num_mds: Number of MDS columns to retain
    - flip: Optional list of columns to flip
    
    Returns:
    - DataFrame with selected MDS columns and computed differences
    """
    # Extract first DataFrame in rds_df
    df = next(iter(rds_df.values())).copy()
    
    # Calculate differences for mds01 and mds02
    df["mds01_diff"] = df["mds01"].diff()
    df["mds02_diff"] = df["mds02"].diff()
    
    # Select relevant columns dynamically up to specified MDS count
    mds_columns = [f"mds{str(i).zfill(2)}" for i in range(1, num_mds + 1)]
    selected_columns = ["chrom", "start", "end", "n", "mid"] + mds_columns + ["mds02_diff", "mds01_diff"]
    df_mds = df[selected_columns].iloc[1:]  # Drop the first row for NaNs from diff
    
    return df_mds

def hmm_mds(df_mds, n_components=3, covariance_type="diag", n_iter=50, random_state=42): 
    """ Applies a Gaussian Hidden Markov Model (HMM) to MDS data to detect hidden states. 
    Parameters: 
    - df_mds: DataFrame containing MDS measurements as features (mds01-mds40). 
    - n_components: Number of hidden states for the HMM (default=3). 
    - covariance_type: Covariance type for the HMM ('diag', 'spherical', 'full', 'tied'). 
    - n_iter: Maximum number of iterations for model fitting (default=50). 
    - random_state: Seed for random number generator (default=42). Returns: 
    - states: Unique hidden states predicted by the model. - Z: Predicted hidden states sequence. - model: Fitted HMM model. 
    """ 
    # Extract MDS features from the DataFrame 
    X = df_mds.iloc[:, 5:].values 
    # Initialize the HMM model 
    model = hmm.GaussianHMM(n_components=n_components, covariance_type=covariance_type, n_iter=n_iter, random_state=random_state) 
    # Fit the model to the data 
    model.fit(X) 
    # Predict the hidden states for the observed MDS data 
    Z = model.predict(X) 
    states = pd.unique(Z) 
    return states, Z, model
    
def flip_regions(df_mds, flip_status=True, flip_start=29305433, flip_end=32378932):
    """
    this script will filp 
    1. MDS data
    2. VCF file
    3. LD data
    """
    if flip_status:
        # Identify the subset of rows to flip
        condition = (df_mds['start'] > flip_start) & (df_mds['start'] < flip_end)
        subset = df_mds[condition]
        # # Columns to leave unchanged
        unchanged_columns = ['chrom','start', 'mid', 'end', 'n']
        # Columns to flip
        columns_to_flip = [col for col in df_mds.columns if col not in unchanged_columns] 
        # Flip the subset while leaving 'start', 'mid', and 'end' unchanged
        subset_flipped = subset.copy()
        subset_flipped[columns_to_flip] = subset[columns_to_flip].iloc[::-1].values   
        # Create a copy of the original dataframe
        df_flipped = df_mds.copy()  
        # Replace the original rows with the flipped subset in the new dataframe
        df_flipped.loc[condition, columns_to_flip] = subset_flipped[columns_to_flip]   
        df2= df_flipped      
    else:
        df2 = df_mds
    return df2

# def plot_mds(df_mds, states, Z, data_path, POP, CHROM):
#     plt.figure(figsize = (16,8))
#     plt.subplot(2,1,1)
#     for i in states:
#         want = (Z == i)
#         x = df_mds["mid"].iloc[want]/1e6
#         y = df_mds["mds01"].iloc[want]
#         plt.plot(x, y, '.')
#     plt.legend(states, fontsize=16)
#     plt.grid(True)
#     plt.xlabel("{}_{}".format(POP, CHROM), fontsize=16)
#     plt.ylabel("MDS01", fontsize=16)
#     plt.subplot(2,1,2)
#     for i in states:
#         want = (Z == i)
#         x = df_mds["mid"].iloc[want]/1e6
#         y = df_mds["mds02"].iloc[want]
#         plt.plot(x, y, '.')
#     plt.legend(states, fontsize=16)
#     plt.grid(True)
#     plt.xlabel("{}_{}".format(POP, CHROM), fontsize=16)
#     plt.ylabel("MDS02", fontsize=16)
#     plt.savefig("{}{}_{}.pdf".format(data_path, POP, CHROM), format='pdf')
#     return None

def plot_mds(df_mds, states, Z, data_path, POP, CHROM):
    plt.figure(figsize=(16, 8))

    # ---------- Plot MDS01 ----------
    plt.subplot(2, 1, 1)
    for i in states:  # loop over unique states
        want = (Z == i)  # 1D boolean mask
        x = df_mds["mid"].to_numpy()[want] / 1e6
        y = df_mds["mds01"].to_numpy()[want]
        plt.plot(x, y, '.', label=f"State {i+1}")
    plt.legend(fontsize=16)
    plt.grid(True)
    plt.xlabel(f"{POP}_{CHROM}", fontsize=16)
    plt.ylabel("MDS01", fontsize=16)

    # ---------- Plot MDS02 ----------
    plt.subplot(2, 1, 2)
    for i in states:
        want = (Z == i)
        x = df_mds["mid"].to_numpy()[want] / 1e6
        y = df_mds["mds02"].to_numpy()[want]
        plt.plot(x, y, '.', label=f"State {i+1}")
    plt.legend(fontsize=16)
    plt.grid(True)
    plt.xlabel(f"{POP}_{CHROM}", fontsize=16)
    plt.ylabel("MDS02", fontsize=16)

    # ---------- Save ----------
    os.makedirs(data_path, exist_ok=True)
    plt.savefig(os.path.join(data_path, f"{POP}_{CHROM}.pdf"), format='pdf')
    plt.show()



def print_consec_rows(df_mds, Z, want_stat=0, min_consec_len=1):
    '''
    # Find the indices of consecutive 2s/0s/1s
    '''
    want_stat = want_stat
    consecutive_twos = np.where((Z[:-1] == want_stat) & (Z[1:] == want_stat))[0]
    
    rows_corresponding_to_twos = df_mds.iloc[consecutive_twos]
    
    lst = rows_corresponding_to_twos['n']
    
    result = []
    for k, g in groupby(enumerate(lst), lambda x: x[0] - x[1]):
        consecutive = list(map(lambda x: x[1], g))
        result.append(consecutive)
    
    # output consecutive rows
    seq_list = []
    for sequence in result:
        seq_list.append(sequence)
    for num,seq in enumerate(seq_list):
        if len(seq) > min_consec_len:
            print(num,seq)
            data_tmp = df_mds[df_mds['n'].isin(seq)]
            START = int(data_tmp.iloc[0]['start'])
            END = int(data_tmp.iloc[-1]['end'])
            print("start: {}, end: {}".format(START, END))
        else:
            continue
    return seq_list
 
# This block ensures that the script can be imported without executing any code
if __name__ == "__main__":
    print("This script is intended to be imported as a module. No direct execution.")