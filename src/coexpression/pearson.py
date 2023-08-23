#Boiler plate
import sys
if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)
         
import os
import numpy as np
import concurrent.futures as cf
import multiprocessing as mp
import scipy

def precalc(expmat_path, delimiter="\t"):
    """Precalculates denominators and nominators for all genes to be used for PCC calculation."""
    genes = []
    nominators = []
    denominators = []
    with open(expmat_path, 'r') as fin:
        for idx, line in enumerate(fin):
            if idx != 0: # assuming first line is header
                parts = line.rstrip().split(delimiter)
                all_values = np.array([float(i) for i in parts[1:]]) # 1 being the geneid so we exclude
                genes.append(parts[1]) # keep geneid info
                nomi = all_values - np.array([(np.sum(all_values)/len(all_values))])
                denomi = np.sqrt(np.sum(nomi**2))
                nominators.append(nomi)
                denominators.append(denomi)
            if idx % 10000:
                print("Precalculations for",idx, "genes completed.")
    
    #convert nominators and denominators into numpy arrays
    nominators = np.array(nominators)
    denominators = np.array(denominators)

    return genes, nominators, denominators

def calc_one_v_all(GOI , genes, nominators, denominators, rank_cutoff = 1000):
    """Give it one gene, it will calculate the PCC vs every gene. Returns corvalues as well as their ranks
     use rank_cutoff=0 for to return values and ranks of all genes.
     rank_cutoff =1000 to only return values and ranks of top 1000 genes."""
    GOI_idx= genes.index(GOI)
    cor_values = np.sum( nominators[GOI_idx] * nominators, axis=1)/(denominators[GOI_idx] * denominators)
    cor_ranks =  scipy.stats.rankdata(cor_values, method="min", nan_policy= "omit") # Smaller cor have smaller (higher) ranking here. 
    cor_ranks =  np.nanmax(cor_ranks) - cor_ranks +1 # flip the ranking . Larger cor have smaller (higher) ranking here. Tied corrs will be assigned maximum possible rank.
    cor_genes = np.array(genes)
    if rank_cutoff==0:
        return cor_values, cor_ranks , cor_genes
    else:
        return cor_values, cor_ranks[cor_ranks<=rank_cutoff] , cor_genes[cor_ranks<=rank_cutoff]


