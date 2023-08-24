#Boiler plate
import sys
if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)
         
import os
import numpy as np
import concurrent.futures as cf
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
                genes.append(parts[0]) # keep geneid info
                nomi = all_values - np.array([(np.sum(all_values)/len(all_values))])
                denomi = np.sqrt(np.sum(nomi**2))
                nominators.append(nomi)
                denominators.append(denomi)
            if idx % 10000==0:
                print("Precalculations for",idx, "genes completed.")
    
    #convert nominators and denominators into numpy arrays
    nominators = np.array(nominators)
    denominators = np.array(denominators)

    return genes, nominators, denominators

def calc_one_v_all(idx, GOI , genes, nominators, denominators, rank_cutoff = 1000):
    """Give it one gene, it will calculate the PCC vs every gene. Returns corvalues as well as their ranks
     use rank_cutoff=0 for to return values and ranks of all genes.
     rank_cutoff =1000 to only return values and ranks of top 1000 genes."""
    GOI_idx= genes.index(GOI)
    cor_values = np.sum( nominators[GOI_idx] * nominators, axis=1)/(denominators[GOI_idx] * denominators)
    cor_ranks =  scipy.stats.rankdata(cor_values, method="min", nan_policy= "omit") # Smaller cor have smaller (higher) ranking here. 
    cor_ranks =  np.nanmax(cor_ranks) - cor_ranks +1 # flip the ranking . Larger cor have smaller (higher) ranking here. Tied corrs will be assigned maximum possible rank.
    cor_genes = np.array(genes)
    if rank_cutoff==0:
        return idx, cor_values, cor_ranks , cor_genes
    else:
        return idx, cor_values[cor_ranks<=rank_cutoff], cor_ranks[cor_ranks<=rank_cutoff] , cor_genes[cor_ranks<=rank_cutoff]


def calc_all_v_all(networkdir, genes, nominators, denominators, rank_cutoff = 1000):
     """Repeating one_v_all for all genes """
     for idx, GOI in enumerate(genes, start = 1):
        idx, cor_values, cor_ranks , cor_genes = calc_one_v_all(GOI , genes, nominators, denominators, rank_cutoff = rank_cutoff)
        with open(os.path.join(networkdir, GOI), "w") as f:
             f.write("Target\tPCC\tRank\n")
             for cv, cr, cg in zip(cor_values,cor_ranks,  cor_genes):
                f.write(f"{cg}\t{cv}\t{cr}\n")
        if idx % 1000==0:
             print("PCCs for",idx , "genes completed.")

def calc_one_v_all_mp(idx, GOI , genes, nominators, denominators, networkdir ,rank_cutoff = 1000):
    """Give it one gene, it will calculate the PCC vs every gene. Returns corvalues as well as their ranks
     use rank_cutoff=0 for to return values and ranks of all genes.
     rank_cutoff =1000 to only return values and ranks of top 1000 genes.
     *FOR MULTIPLE PROCESSING*"""
    GOI_idx= genes.index(GOI)
    cor_values = np.sum( nominators[GOI_idx] * nominators, axis=1)/(denominators[GOI_idx] * denominators)
    cor_ranks =  scipy.stats.rankdata(cor_values, method="min", nan_policy= "omit") # Smaller cor have smaller (higher) ranking here. 
    cor_ranks =  np.nanmax(cor_ranks) - cor_ranks +1 # flip the ranking . Larger cor have smaller (higher) ranking here. Tied corrs will be assigned maximum possible rank.
    cor_genes = np.array(genes)
    
    if rank_cutoff==0:
        pass
    else:
        cor_values, cor_ranks, cor_genes = cor_values[cor_ranks<=rank_cutoff], cor_ranks[cor_ranks<=rank_cutoff] , cor_genes[cor_ranks<=rank_cutoff]
    with open(os.path.join(networkdir, GOI), "w") as f:
        f.write("Target\tPCC\tRank\n")
        for cv, cr, cg in zip(cor_values,cor_ranks,  cor_genes):
            f.write(f"{cg}\t{cv}\t{cr}\n")
    if idx % 1000==0:
        print("PCCs for",idx , "genes completed.")

def calc_all_v_all_mp(networkdir, genes, nominators, denominators, rank_cutoff = 1000, workers = 4):
     """Repeating one_v_all for all genes.
     supports multiple processes"""
     with cf.ProcessPoolExecutor(max_workers=workers) as executor:
         results = [executor.submit(calc_one_v_all_mp, idx, GOI , genes, nominators, denominators, networkdir, rank_cutoff = rank_cutoff) for idx, GOI in enumerate(genes, start = 1)]

     
