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
                      Rank = float(Rank.split("\n")[0])
                      network[source][target] = [PCC, Rank] # 
        if idx%10000 == 0:
            print(idx,"genes loaded")
    return network

def load_PCC_Rank_HRR_MR(networkdir):
    network = {}
    for idx, source in enumerate(os.listdir(networkdir), start =1):
        network[source]={}
        with open(os.path.join(networkdir, source), "r") as f:
            for line_no, line in enumerate(f):
                 if line_no != 0 and line != "": # skip first and last line
                      target, PCC, Rank, HRR, MR = line.split("\t")
                      PCC = float(PCC)
                      Rank = float(Rank)
                      HRR = float(HRR)
                      MR = float(MR.split("\n")[0])
                      network[source][target] = [PCC, Rank , HRR, MR]  
        if idx%10000 == 0:
            print(idx,"genes loaded")
    return network

def add_HRR_MR(network, cutoff= 1000):
    for idx, source in enumerate(list(network.keys())):
        for target in network[source].keys():
            T_S_vals = network.get(target, {}).get(source, [])
            if len(T_S_vals)==2: # not calculated. both source and target perspective present.
                HRR , MR =  np.nanmax([T_S_vals[1], network[source][target][1]]) , scipy.stats.gmean([T_S_vals[1], network[source][target][1]])
                network[source][target].extend([HRR, MR])
                network[target][source].extend([HRR, MR])
            elif len(T_S_vals) <2: # not present from target's perspective. That is okay
                HRR , MR =  np.nanmax([cutoff, network[source][target][1]]) , scipy.stats.gmean([cutoff, network[source][target][1]])
                network[source][target].extend([HRR, MR])
                network[target][source] = [ network[source][target][0],math.nan,HRR, MR]
        if idx % 1000 == 0:
            print("Added MR and HRR for", idx ,"genes")
    return network

def write_HRR_MR(networkdir, PCC_network, headers = ["Target","PCC","Rank", "HRR", "MR"], cutoff= 1000):
    #gonna speed this up somehow
    n=0
    for source , neighbourhood in PCC_network.items():
        if neighbourhood != {}: 
            neighbourhood_unpacked = np.array(list((neighbourhood.values())))
            source_PCCs = neighbourhood_unpacked[:,0]
            source_ranks = neighbourhood_unpacked[:,-1]
            target_ranks = []
            targets = list((neighbourhood.keys()))
            for target in targets:
                try:
                    target_ranks.append(PCC_network[target][source])
                except:
                    target_ranks.append([cutoff,cutoff])
            source_target_ranks = np.array([source_ranks, np.array(target_ranks)[:,-1]])
            HRRs = np.nanmax(source_target_ranks, axis = 0)
            MRs = scipy.stats.mstats.gmean(source_target_ranks, axis = 0, nan_policy="omit")
            with open(os.path.join(networkdir, source) , "w") as f:
                f.write("\t".join(headers)+"\n")
                for target,PCC ,rank, HRR, MR in zip(targets , source_PCCs , source_ranks , HRRs , MRs):
                    f.write(f"{target}\t{PCC}\t{rank}\t{HRR}\t{MR}\n")

        else:
            with open(os.path.join(networkdir, source) , "w") as f:
                f.write("\t".join(headers)+"\n")
        
        n+=1
        if n % 100 ==0:
            print("Added MR and HRR for", n ,"genes")




def dump(network, network_dir, headers = ["Target","PCC","Rank", "HRR", "MR"]):
     for source, neighbourhood in network.items():
          with open(os.join(network_dir, source) , "w") as f:
               f.write("\t".join(headers)+"\n")
               for target, values in neighbourhood.items():
                    f.write( target + "\t"+ "\t".join(values)+"\n")
     

def write_GOIs_table(GOIs, network, table_path, info_dict = {"PCC" : 0, "HRR": 2, "MR": 3 }):
    with open(table_path, "w") as f:     
        f.write("source\ttarget\t" + "\t".join(list(info_dict.keys())) + "\n")
        for idx, source in enumerate(GOIs, start =1):

            for target in GOIs[idx:]:
                info_to_write=[source, target]
                for info_idx in list(info_dict.values()):
                    try:
                        info_to_write.append(str(network[source][target][info_idx]))
                    except:
                        try:
                            info_to_write.append(str(network[target][source][info_idx]))
                        except:
                            pass
                if len(info_to_write) ==5:
                    f.write("\t".join(info_to_write) + "\n")



if __name__ == "__main__":
    #I use this to test out my functions.
    network ={"A":{"B":[0.7,1], "C":[0.6,2]},
          "B":{"A":[0.7,2], "C":[0.8, 1]},
          "C": {"B":[0.8,1]}}
    
    print("Before adding MR_HRR:")
    print(network)

    network = add_HRR_MR(network, cutoff= 5)
    print("After adding MR_HRR:")
    print(network)
    GOIs = ["A", "B"]
    table_path = "/home/ken/table.txt"
    write_GOIs_table(GOIs, network, table_path, info_dict = {"PCC" : 0, "HRR": 2, "MR": 3 })


    