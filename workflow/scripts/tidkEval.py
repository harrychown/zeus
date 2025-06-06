#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import argparse
from pathlib import Path
import pandas as pd
import pandas as pd
import numpy as np

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
    print('Quanitative analysis of telomeres using TIDK output')
    print('')
    print('Usage:')
    print('mandatory:')
    print('-d/--denovo - denovo TIDK output TSV')
    print('optional:')
    print('-o/--outdir - output directory (default = current directory')
    print('-r/--reference - reference TIDK output TSV')
    print('-c/--cutoff - Minimum number of repeats to confirm telomere end status (default = 7)')
    print('-w/--window - Minimum window size (default = 9999)')
    print('-u/--useRef - Use reference file to select for minimum repeat number (recommended if you have not calculated this prior)')

    print('-v/--version for the version')
    print('------------------------------------------------------------------------------------------------------------------------------------------------------')


"""
Identify contig ends
"""
def evalTelomeres(tidk_table, repeat_cutoff=0, window_cutoff=500):
    # Sum values for either forward / reverse count
    tidk_table['count'] = tidk_table[['forward_repeat_number', 'reverse_repeat_number']].sum(axis=1)
    unique_contig_id = set(tidk_table["id"])
    # Calculate repeats at ends of each contig
    contig_ends = []
    for i in unique_contig_id:
        contig_results = tidk_table[tidk_table['id'] == i]
        contig_start = contig_results[contig_results['window'] == min(contig_results['window'])]
        contig_finish = contig_results[contig_results['window'] == max(contig_results['window'])]
        contig_ends.append(contig_start)
        contig_ends.append(contig_finish)
    tidk_ends = pd.concat(contig_ends)
    
    # Report number of contigs with no telomeric repeats
    telomeres_present = tidk_ends[(tidk_ends['count'] > repeat_cutoff) & (tidk_ends['window'] > window_cutoff)]
    num_zero_telomeres = len(tidk_ends) - len(telomeres_present)
    # Calculate total repeat count
    total_telomeric_repeat = sum(telomeres_present['count'])
    total_repeat = sum(tidk_ends['count'])
    # Calculate minimum repeat size
    min_repeat = min(telomeres_present['count'])
    # Calculate median repeat count
    med_repeat = np.median(telomeres_present['count'])
    # Output results
    telomere_results = {"num_telomere_ends":[len(telomeres_present)], 
                        "num_non-telomere_ends":[num_zero_telomeres], 
                        "total_repeats":[total_repeat], 
                        "total_telomeric_repeats":[total_telomeric_repeat], 
                        "minimum_repeat_size":[min_repeat], 
                        "median_repeat_size":[med_repeat]}
    return(telomere_results)


"""
Run local
"""
"""
reference_file = "/Users/user/OneDrive/Documents/assembly/scripts/sandbox/tidk_check/ref_telomeric_repeat_windows.tsv"
reference_tidk_df = pd.read_csv(reference_file, sep='\t')
reference = evalTelomeres(reference_tidk_df, repeat_cutoff=0, window_cutoff=9999)
"""

def main():
    #Parse all arguments
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-d', '--denovo', type=str, default=None, required=True, help='Input denovo TIDK TSV')
    parser.add_argument('-o', '--outdir', type=str, default='.', required=False, help='Directory for output files (default is current working directory')
    parser.add_argument('-r', '--reference', type=str, default=None, required=False, help='Input reference TIDK TSV')
    parser.add_argument('-c', '--cutoff', type=int, default=7, required=False, help='Minimum number of repeats to confirm telomere end status (default = 7)')
    parser.add_argument('-w', '--window', type=int, default=9999, required=False, help='Minimum window size (default = 9999)')
    parser.add_argument('-u', '--useRef', action='store_true', required=False, help='Use reference file to select for minimum repeat number (recommended if you have not calculated this prior)')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.0.1')
    parser.add_argument('-p', '--prefix', type=str, default=None, required=False, help='Prefix for denovo results')

    args = parser.parse_args()
    

    # Set the input/output
    denovo_file = Path(args.denovo)
    outdir = Path(args.outdir)
    # Open input
    denovo_tidk_df = pd.read_csv(denovo_file, sep='\t')
    # Evaluate reference telomeres
    if args.reference is not None:
        reference_file = Path(args.reference)
        reference_tidk_df = pd.read_csv(reference_file, sep='\t')
        if args.useRef:
            reference = evalTelomeres(reference_tidk_df, repeat_cutoff=0, window_cutoff=args.window)
        else:
            reference = evalTelomeres(reference_tidk_df, repeat_cutoff=args.cutoff, window_cutoff=args.window)
        # Store reference output
        ref_df = pd.DataFrame(reference)
        ref_df.to_csv(str(outdir) + "/reference.telomeres.csv", index=False, sep=',')
        print(ref_df)
    
    # Evaluate denovo telomeres
    if args.useRef:
        denovo = evalTelomeres(denovo_tidk_df, repeat_cutoff=(reference["minimum_repeat_size"][0] - 1), window_cutoff=args.window)
    else:
        denovo = evalTelomeres(denovo_tidk_df, repeat_cutoff=args.cutoff, window_cutoff=args.window)
    denovo_df = pd.DataFrame(denovo)
    if args.prefix is not None:
        denovo_df.to_csv(str(outdir) + "/" + args.prefix + ".denovo.telomeres.csv", index=False, sep=',')
    else:
        denovo_df.to_csv(str(outdir) + "/denovo.telomeres.csv", index=False, sep=',')
    #print(denovo_df)
if __name__ == '__main__':
    main()
