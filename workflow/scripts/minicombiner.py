#!/usr/bin/python

import sys
import argparse
from pathlib import Path
import pandas as pd

def quit(message=None):
    """Exit the program with an optional message."""
    if message:
        print(message)
    sys.exit(1)

def usage():
    """ Usage of the script"""
    print('------------------------------------------------------------------------------------------------------------------------------------------------------')
    print('Harry Chown h.chown@imperial.ac.uk')
    print('')
    print('Combining Miniprot GFF to enable synteny analysis with DAGChainer')
    print('')
    print('Usage:')
    print('mandatory:')
    print('-r/--reference - reference Miniprot output GFF')
    print('-d/--denovo - denovo Miniprot output GFF')
    print('-o/--outfilename - output file name')
    print('optional:')
    print('-v/--version for the version')
    print('------------------------------------------------------------------------------------------------------------------------------------------------------')


'''
Read GFF and convert to dictionary
'''
def extractLocation(in_file):
    gff_lines = open(in_file).readlines()
    gene_location = {}
    for line in gff_lines:
        if "#" in line:
            continue
        separated_line = line.split('\t')
        if "mRNA" in separated_line[2]:
            contig=separated_line[0]
            start=separated_line[3]
            stop=separated_line[4]
            gene=separated_line[8].split(" ")[0].split("=")[-1]
            location=[contig,gene,start,stop]
            if gene_location.get(gene, "new") == "new":
                gene_location[gene] = [location]
            else:
                gene_location[gene] = gene_location.get(gene) + [location]
    return(gene_location)

'''
Identify matching gene pairs and save output
'''
def generatePairs(dict1, dict2, out_filename):
    pair_list = []
    for k, rv in dict1.items():
        if k in dict2.keys():
            dv = dict2.get(k)
            for r_loc in rv:
                for d_loc in dv:
                    #Add reference prefix
                    r_loc[1] = "ref_" + r_loc[1]
                    pair = r_loc + d_loc + ["0"]
                    pair_list.append(pair)
                    
    pair_df = pd.DataFrame(pair_list)
    pair_df.to_csv(out_filename, sep='\t', header=False, index=False)
    



def main():
    #Parse all arguments
    parser = argparse.ArgumentParser(description='Runs FUR')
    parser.add_argument('-r', '--reference', type=str, default=None, required=True, help='Input reference genome GFF in the formatted output of Miniprot')
    parser.add_argument('-d', '--denovo', type=str, default=None, required=True, help='Input denovo genome GFF in the formatted output of Miniprot')
    parser.add_argument('-o', '--outfilename', type=str, default=None, required=True, help='Filename of output')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.0.1')

    args = parser.parse_args()
    

    # Set the input/output
    reference_gff = Path(args.reference)
    denovo_gff = Path(args.denovo)
    outfilename = Path(args.outfilename)
    
    # Calculate gene pairs
    reference_loc = extractLocation(reference_gff)
    denovo_loc = extractLocation(denovo_gff)

    generatePairs(reference_loc, denovo_loc, outfilename)




if __name__ == '__main__':
    main()