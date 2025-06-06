#!/usr/bin/python

import sys
import argparse
from pathlib import Path
import pandas as pd
import math
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
    print('Summarising Tandem Repeat Finder Results')
    print('')
    print('Usage:')
    print('mandatory:')
    print('-i/--input - Tandem repeat finder output in ".dat" format')
    print('optional:')
    print('-o/--outdir - Output directory (default: current working directory')
    print('-abs/--absolute - Absolute length cutoff (>x) (default: 23)')
    print('-unb/--unit - Repeat unit boundary between short/long (default: 40)')
    print('-cb/--copy - Unit copy number boundary between low/high (default: 5)')

    print('-p/--prefix - Output prefix')
    print('-v/--version for the version')
    print('------------------------------------------------------------------------------------------------------------------------------------------------------')


# Filename

def reformatTRF(infile, nt_cutoff=23, unit_boundary=40, copy_boundary=5):
    with open(infile, "r") as fh:
        input_data = fh.readlines()
    # Init gff line list
    trf_list = list()
    # Preprocess incomming data to remove blank lines
    lines = [line for line in (ln.rstrip() for ln in input_data) if line]
    for line in lines:
        # Split line in space delimited elements
        ele = line.strip().split(" ")
        # Check if begining of new sequence and update current name
        if line.startswith("Sequence:"):
            contig = str(ele[1])
        # Else check if data line and unpack fields
        elif ele[0][0].isdigit():
            # Dat file fields
            # rep_start, rep_end, period_size, copy_no, pattern_size,
            # percent_match, percent_indel, align_score, a_percent,
            # c_percent, g_percent, t_percent, entropy, consensus, repeat
            [
                start,
                stop,
                period,
                copies,
                consensus_size,
                perc_match,
                perc_indels,
                align_score,
                perc_A,
                perc_C,
                perc_G,
                perc_T,
                entropy,
                cons_seq,
                repeat_seq,
            ] = ele
    
            # Changing format
            # Calculate the absolute sequence length
            absolute_length=math.floor(int(period) * float(copies))
            # Define category
            repeat_unit_category = "NA"
            repeat_unit_copy_number_category = "NA"
            category = "NA"
            if int(period) == 1:
                category = "Homopolymer"
            elif (int(period) > 2) & (absolute_length > nt_cutoff):
                if (int(period) < unit_boundary) & (float(copies) < copy_boundary):
                    category = "Short-Low"
                    repeat_unit_category = "Short"
                    repeat_unit_copy_number_category = "Low"
                elif (int(period) >= unit_boundary) & (float(copies) < copy_boundary):
                    category = "Long-Low"
                    repeat_unit_category = "Long"
                    repeat_unit_copy_number_category = "Low"
                elif (int(period) < unit_boundary) & (float(copies) >= copy_boundary):
                    category = "Short-High"
                    repeat_unit_category = "Short"
                    repeat_unit_copy_number_category = "High"
                elif (int(period) >= unit_boundary) & (float(copies) >= copy_boundary):
                    category = "Long-High"
                    repeat_unit_category = "Long"
                    repeat_unit_copy_number_category = "High"
    
                
            # Convert to gff line
            df_row = {
                "contig":contig,
                "start":int(start),
                "stop":int(stop),
                "repeat_unit_length":int(period),
                "copy_number":float(copies),
                "absolute_length":int(absolute_length),
                "perc_match":float(perc_match),
                "perc_indels":float(perc_indels),
                "align_score":float(align_score),
                "entropy":float(entropy),
                "cons_seq":cons_seq,
                "repeat_seq":repeat_seq,
                "category":category,
                "repeat_unit_category":repeat_unit_category,
                "repeat_unit_copy_number_category":repeat_unit_copy_number_category
                
            }
            trf_list.append(df_row)
    out_df=pd.DataFrame(trf_list)
    return(out_df)



def main():
    #Parse all arguments
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--input', type=str, default=None, required=True, help='Tandem repeat finder output in ".dat" format')
    parser.add_argument('-o', '--outdir', type=str, default='.', required=False, help='Directory for output files (default is current working directory')
    parser.add_argument('-abs', '--absolute', type=int, default=23, required=False, help='Absolute length cutoff (>x) (default: 23)')
    parser.add_argument('-unb', '--unit', type=int, default=40, required=False, help='Repeat unit boundary between short/long (default: 40)')
    parser.add_argument('-cb', '--copy', type=int, default=5, required=False, help='Unit copy number boundary between low/high (default: 5)')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.0.1')
    parser.add_argument('-p', '--prefix', type=str, default=None, required=False, help='Prefix for results')

    args = parser.parse_args()
        
    
    # Set the input/output
    infile = Path(args.input)
    outdir = Path(args.outdir)
        
    # Reformat TRF output
    trf_df=reformatTRF(infile, nt_cutoff=args.absolute, unit_boundary=args.unit, copy_boundary=args.copy)
    # Calculate categorical counts
    trf_filtered=trf_df[trf_df['category'] != 'NA']
    trf_summary=trf_filtered['category'].value_counts()
    # Output summary and reformatted file
    summary_df=pd.DataFrame(trf_summary).transpose()
    if args.prefix is not None:
        summary_df.to_csv(str(outdir) + "/" + args.prefix + ".TRF-summary.csv", index=False, sep=',')
        trf_df.to_csv(str(outdir) + "/" + args.prefix + ".TRF-reformat.csv", index=False, sep=',')
    else:
        summary_df.to_csv(str(outdir) + "/TRF-summary.csv", index=False, sep=',')
        trf_df.to_csv(str(outdir) + "/TRF-reformat.csv", index=False, sep=',')


if __name__ == "__main__":
    main()
