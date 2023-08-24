#Boiler plate
import sys
if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)
         
import os
import argparse
#mutwilab module imports
from data_processing import read_write , network

if __name__ == "__main__":
    parser= argparse.ArgumentParser(description="Get_HRR_MR.py.\n\
                                   Load PCC network and calculate HRR and MR.")
    

    parser.add_argument("-o", "--output_dir", type=str, metavar="", required = True,
    help = "Directory to read PCC network. Will also output data there." )
    
    parser.add_argument("-rank", "--rank_to_retain", type= int, metavar="", default = 10000,
    help = "Must be same as the paramenter you used in  PCC_net_from_expmat.py")

    parser.add_argument("-GOIs", "--path_to_GOIs", type= str, metavar="", default = "",
    help = "path to text file containing GOIs to generate network table for cytoscape. If not specified, the network table will not be generated.")

    parser.add_argument("-skip", "--skip_HRR_MR", action="store_true",
    help = "Skip calculating MR and HRR, and load network with precalculated MR and HRR values. Will just generate table instead.")

    args=parser.parse_args()
    output_dir ,rank_to_retain , path_to_GOIs =  args.output_dir , args.rank_to_retain, args.path_to_GOIs
    skip_HRR_MR = args.skip_HRR_MR
    
    PCC_network_dir = os.path.join(output_dir, "PCC_network")
    MR_HRR_network_dir = os.path.join(output_dir, "PCC_network_w_MR_HRR")
    read_write.establish_dir(MR_HRR_network_dir, isdir= True)

    if skip_HRR_MR:
        print("Skipping adding of MR, HRR. Loading network that has MR and HRR already calculated into memory...")
        MR_HRR_network = network.load_PCC_Rank_HRR_MR(MR_HRR_network_dir)
    else:    
        print("Loading PCC_network into memory...")
        PCC_network = network.load_PCC_Rank(PCC_network_dir)
        print("Calculating MR and HRR values and writing to output directory...")
        network.write_HRR_MR(MR_HRR_network_dir, PCC_network, headers = ["Target","PCC","Rank", "HRR", "MR"], cutoff= rank_to_retain)
        #MR_HRR_network = network.add_HRR_MR(PCC_network, cutoff= rank_to_retain)
        print("Loading network with MR and HRR into memory... ")
        MR_HRR_network = network.load_PCC_Rank_HRR_MR(MR_HRR_network_dir)

    
    if path_to_GOIs != "":
        with open(path_to_GOIs, "r") as f:
            GOIs = f.read()
        GOIs = GOIs.split("\n")
        GOIs = [GOI for GOI in GOIs if GOI !=""]
        
        print(f"\nNumber of GOIs in {path_to_GOIs}:\n{len(GOIs)}\n")
        print(f"Example of GOIs specified by user:\n{GOIs[:3]}\n")

        table_path = os.path.join(output_dir, path_to_GOIs.split("/")[-1].split(".")[0]+"_network_table.tsv")
        print("Building network table from network...")
        network.write_GOIs_table(GOIs, MR_HRR_network, table_path, info_dict = {"PCC" : 0, "HRR": 2, "MR": 3 })
        print("Get_HRR_MR.py complete.\n\nif you want to generate another network table for this species, you can skip the calculation of MR and HRR ")




