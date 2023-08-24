#Boiler plate
import sys
if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)
         
import os
import argparse
#mutwilab module imports
from coexpression import pearson  
from data_processing import read_write

if __name__ == "__main__":
    parser= argparse.ArgumentParser(description="PCC_net_from_expmat.py.\n\
                                    Generate PCC network from expression matrix.")
    
    parser.add_argument("-in", "--input_expmat_path", type= str, metavar="", required= True,
    help = "Input expression matrix. where columns: samples and rows: tpm values of genes).")

    parser.add_argument("-o", "--output_dir", type=str, metavar="", required = True,
    help = "Directory to output network into." )

    parser.add_argument("-w", "--workers", type=int, metavar="", default = 20,
    help = "Number of concurrent calculations. Note: larger number of concurrent processes will consume more RAM." )
    
    parser.add_argument("-de", "--delimiter", type= str, metavar="", default = "t", choices=["t", "c"],
    help = "Delimiter for expression matrix. -de=\"t\" for tab seperated (.tsv). -de=\"c\" for comma seperated (.csv). TSV by default.")

    parser.add_argument("-rank", "--rank_to_retain", type= int, metavar="", default = 10000,
    help = "If rank= 1000, only the first 1000 PCCs of a given gene will be written into the network. Directly limits the HRR cutoff in Get_HRR_MR.py")

    args=parser.parse_args()
    
    output_dir, delimiter ,rank_to_retain , expmat_path =  args.output_dir , args.delimiter, args.rank_to_retain, args.input_expmat_path
    workers = args.workers

    if delimiter == "t":
        delim = "\t"
    else:
        delim = ","

    #establish output and sub directories
    sub_outdir = os.path.join(output_dir , "PCC_network")
    read_write.establish_dir(sub_outdir , isdir=True)

    #print info and start
    print(f"Building PCC network with these specifications:\n\
          input_expmat_path={expmat_path}\n\
          rank_to_retain={rank_to_retain}\n\n\
            And saving it to {sub_outdir}\n")
    
    print("Running precalculations...")
    genes, nominators, denominators = pearson.precalc(expmat_path, delimiter=delim)

    print("Calculating PCCs for every gene...")
    pearson.calc_all_v_all_mp(sub_outdir, genes, nominators, denominators, rank_cutoff = 1000, workers = workers)
    
