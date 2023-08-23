#Boiler plate
import sys
if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)
import os
import scipy
import numpy as np
import math

def load_PCC_Rank(networkdir):
    network = {}
    for idx, source in enumerate(os.listdir(networkdir), start =1):
        network[source]={}
        with open(os.path.join(networkdir, source), "r") as f:
            for line_no, line in enumerate(f):
                 if line_no != 0 and line != "": # skip first and last line
                      target, PCC, Rank = line.split("\t")
                      PCC = float(PCC)
                      Rank = int(Rank.split("\n")[0])
                      network[source][target] = [PCC, Rank] # 
        if idx%1000 == 0:
            print(idx,"genes loaded")
    return network

def add_HRR_MR(network, cutoff= 1000):
    for source in network.keys():
        for target in network[source].keys():
            T_S_vals = network.get(target, {}).get(source, [])
            if len(T_S_vals)==2: # not calculated. both source and target perspective present.
                HRR , MR =  np.nanmax([T_S_vals[1], network[source][target][1]]) , scipy.stats.gmean([T_S_vals[1], network[source][target][1]])
                network[source][target].extend([HRR, MR])
                network[target][source].expend([HRR, MR])
            elif len(T_S_vals) >2: # already calulated. skip
                pass
            elif len(T_S_vals) <2: # not present from target's perspective. That is okay
                HRR , MR =  np.nanmax([cutoff, network[source][target][1]]) , scipy.stats.gmean([cutoff, network[source][target][1]])
                network[source][target].extend([HRR, MR])
                network[target][source] = [ math.nan,math.nan,HRR, MR]
    return network

            
network ={"A":{"B":[0.7,1], "C":[0.6,2]},
          "B":{"A":[0.7,2], "C":[0.8, 1]},
          "C": {"B":[0.8,1]}}
network = add_HRR_MR(network, cutoff= 5)