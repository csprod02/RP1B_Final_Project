import argparse

import pandas as pd
import numpy as np


# read parsnp vcf output file and remove any snps that are presesnt in 
def singletons_from_parsnp_vcf(file_path, job_name):
    
    with open(file_path) as in_file, open(f'singletons_{job_name}.csv', 'w') as out_file: # file path = parsnp.vcf

        # read all extra header lines until the actual header line
        # then stop reading the lines and move on (break)
        for line in in_file:
            if line.startswith('#CHROM'):
                out_file.write(line)
                break

        snps = 0

        # each line represents a snp
        for line in in_file:

            snps += 1

            # this filters to include ONLY snps with one variant (i.e. removes multialleleic snps)
            if len(line.split()[4].strip()) > 1:
                continue
            else:
                samples = line.split()[10:] # these are all the sample columns for each line
                snp_chrom = line.split()[0] # this is the contig where the snp is found
                snp_pos = line.split()[1] # this is the position on the contig where the snp is found

                # for each sample, if it differs from the reference (0 for all snps), it will have a value of 1
                # so if only one sample has the snp then the sum of values for that snp will be 1
                # if all samples have a 1 then this is still a singleton because the ref seq will have a 0
                if ((sum(map(int, samples)) == 1) and (line.split()[6] == 'PASS')):
                    # keep these
                    out_file.write(line)

    #            # CONVERTING the 999 lines so that 0 and 1 are swapped so they're in the same format as the other lines
                elif ((sum(map(int, samples)) == 999) and (line.split()[6] == 'PASS')):
                    modified_line = line[:9].join('1' if c == '0' else '0' for c in line[9:])
                    lines.append(modified_line)

                else:
                    # remove
                    continue
#        print(snps)
#        print()
                    
    return


# Creates a dictionary containing all snps from the singletons.csv file
# {'sample1': {'contig1': [array of pos],
#               'contig2': [array of pos],
#               'contig3': [array of pos]},
#  'sample2': {'contig1': [array of pos]},
#  'sample3': {'contig3': [array of pos]}
def create_snp_dict(file):

    samples = file.readline().split()[10:]

    snp_dict = {}

    for sample in samples:
        snp_dict[sample] = {}


    for line in file:
        l = line.split()
        contig = l[0]
        values = l[10:] # samples start from line[10] onwards
        pos = l[1]

        for i, current_sample in enumerate(samples):
            if values[i] == '1':
                if contig not in snp_dict[current_sample]:
                    snp_dict[current_sample][contig] = []

                snp_dict[current_sample][contig].append(pos)
                    
    return snp_dict


# Uses a sliding window to indentify areas of recombination
# Looks through the entire vcf
# Scans every sample (and every contig for every sample)
def identify_recombination(snp_dict, num_snp_threshold, window):
    
    recomb_snps = []
    recomb_sets = set()

    for sample in snp_dict:
        window_size = int(window*1000)
        num_snps_threshold = int(num_snp_threshold)

        for contig in snp_dict[sample]:
            positions = snp_dict[sample][contig]
            index = 0

            # This makes sure that the index will never exceed the number of snps for the specific contig (for that specific sample)
            # Equivalent to using a end_of_contig = True/False situation- it automatically is False and moves on as soon as the index exceeds the end of the contig
            while index < len(positions):
                # Start a new window at current index
                window_start = int(positions[index])
                window_end = window_start + window_size
                start_index = index # This keeps track of the start of this window so that it can be used when compiling array of all the snps in a cluster

                # Now move the index along, one snp at a time to count all the snps within the window
                # While the current index is within the window AND the position at that index is within the window, increase the index by 1
                while index < len(positions) and int(positions[index]) <= window_end:
                    index += 1

                # update index so it starts at the first snp outside of the current window
                snp_count = index - start_index

                # store the cluster if num snps is larger than 
                if snp_count >= num_snps_threshold:
                    recomb_snps.append({
                        'numSNPs': snp_count,
                        'sample': sample,
                        'chrom': contig,
                        'pos': positions[start_index:index]})
                    
                    recomb_sets.add((contig, ))

                # Move window forward by ONE SNP (not to the end of the cluster)
                else:
                    index = start_index + 1
                
    return recomb_snps


# Creates a set with all the unique contig, pos combinations that have been flagged as recombinant snps
# This is used later when filtering the singletons.csv to exclude these reocmbinant snps
def create_snp_sets(recomb_snps):
    
    recomb_sets = set()

    for recomb_group in recomb_snps:
        contig = recomb_group['chrom']
        pos_arr = recomb_group['pos']

        for pos in pos_arr:
            recomb_sets.add((contig, pos))
            
    return recomb_sets
    

# Filters the singletons.csv file and writes the non-recombinant snps to a new file    
def write_filtered_snps_to_file(recomb_sets, in_file, out_file):
    
    header_line = in_file.readline()

    out_file.write(header_line)

    for line in in_file:
        contig = line.split()[0]
        pos = line.split()[1]

        if (contig, pos) not in recomb_sets:
            out_file.write(line)
            
    return
                


# reads the ref seq file and returns the sequence as well as a dict containing the sequence covered by each contig for ref sequences that are made of multiple contigs
def read_refseq_file(file_path):
    ref_dict = {}
    
    with open(file_path, 'r') as file:
        
        sequence = ''
        
        for line in file:
            if not line.startswith('>'):

                sequence += line.strip()
                
        # resetting the cursor
        # creating dict with contig name as keys and sequence as values
        # this will be used later on, but creating this now so that I don't have to open/close the file again later
        file.seek(0)
        
        for line in file:
            if line.startswith('>'):
                
                if line not in ref_dict:
                    prev_contig = line[1:].split()[0].strip()
                    ref_dict[prev_contig] = ''
                    
            else:
                ref_dict[prev_contig] += line.strip()
                

                    
    return sequence, ref_dict


# Counts the number of each nucleotide in the ref sequence and stores as a dict
def count_refseq_nucleotides(ref_seq):
    
    counts = {}
    
    for n in ref_seq:
        
        if n not in counts:
            counts[n] = 1
        else:
            counts[n] += 1
            
            
    return counts



# calculates the frequencies of each of the bases in the ref seq
def ref_base_frequency(ref_base_counts):

    freqs = {}

    total_ref_bases = sum(ref_base_counts.values())
    
    for base in ref_base_counts:
        freqs[base] = ref_base_counts[base]/total_ref_bases

    return freqs
    

    

# finds all the triplets in the ref sequence and counts how many of each there are.
# Stores this as a dict
def get_ref_triplets(ref_seq):

    ref_3mers = {}


    for i in range(len(ref_seq)-2):
        triplet = ref_seq[i:i+3]

        if triplet not in ref_3mers:
            ref_3mers[triplet] = 1
        else:
            ref_3mers[triplet] += 1


    return ref_3mers
        

    

# returns the reverse complement of a DNA string
def revC(dna_seq):
    # Flips the entire String i.e. hello -> olleh
    reverse_seq = dna_seq[::-1]

    # Creates an empty string to hold the reverse complement sequence/String
    reverse_complement_seq = ''

    # Createa a dictionary containing the complementary base
    complement = {'A': 'T',
                 'T': 'A',
                 'G': 'C',
                 'C': 'G'}

    # For each base in the reversed sequence, add its complement to the reverse_complement String
    for i in reverse_seq:
        reverse_complement_seq += complement[i]

    return reverse_complement_seq




def get_context_for_df(row, ref_dict):
    
    chrom = row['#CHROM']
    pos = int(row['POS'])
    
    contig_seq = ref_dict[chrom]

    return contig_seq[pos-2:pos+1]


def get_norm_for_row(row, ref_n_counts):

    ref = row['collapsed_group'][0]
    alt = row['collapsed_group'][2]
    count = row['count']


    if ref in ['C', 'G']:
        ref_total = int(ref_n_counts['C']) + int(ref_n_counts['G'])

    elif ref in ['T', 'A']:
        ref_total = int(ref_n_counts['T']) + int(ref_n_counts['A'])

    norm = (count/ref_total) *100

    return norm


def get_context_norm_for_row(row, ref_triN_counts):

    count = row['count']
    
    # each collapsed group is like this: W[X>Y]Z
    group = row['context_collapsed_group']
    
    five_prime_base = group[0]
    centre_base = group[2]
    three_prime_base = group[-1]

    ref_triN = five_prime_base + centre_base + three_prime_base

    complement_triN = revC(ref_triN)

    total_ref_triN = ref_triN_counts[ref_triN] + ref_triN_counts[complement_triN]

    norm = (count/total_ref_triN) *100

    return norm
    
    
def create_sample_df(ref_dict, ref_n, ref_triN, job):

    # reads the singletons_no_recombination csv into a pandas dataframe
    df = pd.read_csv(f'singletons_no_recombination_{job}.csv', sep='\t')

    # cleaning up the dataframe.  Dropping columns that are not used
    ref_col = df.columns[9]
    df.drop(columns=ref_col, inplace=True)
    df.drop(['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], axis=1, inplace=True)


    # meta_cols are the columns that I want to keep in the dataframe
    # the sample cols are all the columns to the right of the cols in meta_cols
    meta_cols = df.columns[0:4]
    sample_cols = df.columns[4:]

    # melts the dataframe so that for each snp has its own row in the dataframe for every sample
    # The has_snp col is the binary 0/1 from the vcf
    df_melted = pd.melt(df,
                        id_vars=meta_cols,
                        value_vars=sample_cols,
                        var_name='sample',
                        value_name='count'
                       )

    # now drop the rows where the has_snp is 0/ only keeping rows with snps
    # This means that the dataframe is now every snp + sample pairing as its own row
    df_melted = df_melted[df_melted['count'] == 1]

    snp_matrix = df_melted.pivot_table(
        index='REF',        # rows
        columns='ALT',      # columns
        values='count',   # values to aggregate
        aggfunc='sum',      # sum the SNP counts
        fill_value=0        # replace missing combos with 0
        )

    print(snp_matrix)
    print()
    
    snp_matrix.to_csv(f'snp_matrix_{job}.csv')

    # adding a context column by using the get_context_for_df function on each row of the dataframe
    df_melted['context'] = df_melted.apply(get_context_for_df, ref_dict=ref_dict, axis=1)


    conditions = [
        # C>X
        (df_melted['REF'] == 'C') & (df_melted['ALT'] == 'A'),
        (df_melted['REF'] == 'G') & (df_melted['ALT'] == 'T'),

        (df_melted['REF'] == 'C') & (df_melted['ALT'] == 'G'),
        (df_melted['REF'] == 'G') & (df_melted['ALT'] == 'C'),

        (df_melted['REF'] == 'C') & (df_melted['ALT'] == 'T'),
        (df_melted['REF'] == 'G') & (df_melted['ALT'] == 'A'),

        # T>X
        (df_melted['REF'] == 'T') & (df_melted['ALT'] == 'A'),
        (df_melted['REF'] == 'A') & (df_melted['ALT'] == 'T'),

        
        (df_melted['REF'] == 'T') & (df_melted['ALT'] == 'C'),
        (df_melted['REF'] == 'A') & (df_melted['ALT'] == 'G'),

        (df_melted['REF'] == 'T') & (df_melted['ALT'] == 'G'),
        (df_melted['REF'] == 'A') & (df_melted['ALT'] == 'C')
        ]

    values = ['C>A', 'C>A',
              'C>G', 'C>G',
              'C>T', 'C>T',

              'T>A', 'T>A',
              'T>C', 'T>C',
              'T>G', 'T>G'              
              ]


    df_melted['collapsed_group'] = np.select(conditions, values, default='ERROR')
    


    conditions = [
        df_melted['context'].str[1].isin(['C', 'T']),
        (df_melted['context'].str[1] == 'G'),
        (df_melted['context'].str[1] == 'A')
        ]

    values = [(df_melted['context'].str[0] + '[' + df_melted['context'].str[1] + '>' + df_melted['ALT'] + ']' + df_melted['context'].str[2]),
              
              (df_melted['context'].apply(revC).str[0] + '[' + 'C' + '>' + df_melted['ALT'].apply(revC) + ']' + df_melted['context'].apply(revC).str[2]),
              
              (df_melted['context'].apply(revC).str[0] + '[' + 'T' + '>' + df_melted['ALT'].apply(revC) + ']' + df_melted['context'].apply(revC).str[2])]


    df_melted['context_collapsed_group'] = np.select(conditions, values, default='ERROR')


    print('DATA FRAME')
    print(df_melted.head(15))
    print()

    # creating a new dataframe using groupby to group the collapsed group
    # this also sums the number of snps (has_snp col)
    norm_df = df_melted.groupby('collapsed_group', as_index=False)[['count']].sum()

    print('NORM DF')
    print(norm_df.head(15))
    print()

    # adds a column with the normalised values after applying the get_norm_for_row function to the dataframe
    norm_df['norm'] = norm_df.apply(get_norm_for_row, ref_n_counts=ref_n, axis=1)

    # adds a column with the fraction of the normalised total mutations that each mutation type accounts for
    norm_df['frac'] = (norm_df['norm']/(norm_df['norm'].sum())) *100

    # adding transition/transversion column to dfs
    transversion = ['C>A', 'C>G', 'T>A', 'T>G']
    transition = ['C>T', 'T>C']

    conditions = [
        norm_df['collapsed_group'].isin(transversion),
        norm_df['collapsed_group'].isin(transition)
        ]

    values = ['transversion', 'transition']

    norm_df['mutation_type'] = np.select(conditions, values, default='ERROR')

    print('NORM DF')
    print(norm_df.head(15))
    print()

    # applying the same process to the context dataframe
    context_norm_df = df_melted.groupby('context_collapsed_group', as_index=False)[['count']].sum()

    print('CONTEXT NORM DF')
    print(context_norm_df.head(15))
    print('---------------')


    context_norm_df['norm'] = context_norm_df.apply(get_context_norm_for_row, ref_triN_counts=ref_triN, axis=1)

    context_norm_df['frac'] = (context_norm_df['norm']/context_norm_df['norm'].sum()) *100

    # adding transition v transversion column
    conditions = [
        context_norm_df['context_collapsed_group'].str[2:5].isin(transversion),
        context_norm_df['context_collapsed_group'].str[2:5].isin(transition)
        ]

    values = ['transversion', 'transition']

    context_norm_df['mutation_type'] = np.select(conditions, values, default='ERROR')

    # adding some extra colums for easier plotting later on
    context_norm_df['flanking'] = context_norm_df['context_collapsed_group'].str[0] + '.' + context_norm_df['context_collapsed_group'].str[-1]
    context_norm_df['fivePrime'] = context_norm_df['flanking'].str[0]
    context_norm_df['threePrime'] = context_norm_df['flanking'].str[-1]
    context_norm_df['sb_mutation'] = context_norm_df['context_collapsed_group'].str[2:5]
    context_norm_df['flanking_equal'] = context_norm_df.apply(lambda x: x['flanking'][0] == x['flanking'][-1], axis=1)
    context_norm_df['flanking_equal_mut'] = context_norm_df.apply(lambda x: (x['flanking_equal'] == True) and (x['flanking'][2] == x['flanking'][0]), axis=1)
    

    print('CONTEXT NORM DF')
    print(context_norm_df.head(15))
    

    # saving all 3 dfs to csv
    df_melted.to_csv(f'melted_df_{job}.csv', index=False)
    norm_df.to_csv(f'snp_count_df_{job}.csv', index=False)
    context_norm_df.to_csv(f'context_count_df_{job}.csv', index=False)
    

    # additional columns added so that I can normalise the values in pivot table easier
    df_melted['opportunity'] = df_melted['context'].map(ref_triN)
    df_melted['norm'] = (df_melted['count']/df_melted['opportunity']) *100
    

    # creating pivot table for pca etc
    spectrum = df_melted.pivot_table(
        index='sample', 
        columns='context_collapsed_group', 
        values='count',
        aggfunc='sum',
        fill_value=0
        )

    spectrum = spectrum.reset_index()

    spectrum_freq = df_melted.pivot_table(
        index='sample', 
        columns='context_collapsed_group', 
        values='norm',
        aggfunc='sum',
        fill_value=0
        )

    
    spectrum_freq.to_csv(f'mutation_freqs_per_sample_{job}.csv')
    
    print(spectrum)

    spectrum.to_csv(f'mutation_counts_per_sample_{job}.csv', index=True)

    print(spectrum_freq)

    spectrum_freq
    

    return df_melted

    
    
# parses arguments from terminal
def parse_args():
    parser = argparse.ArgumentParser(description='Count SNPs')
    
    parser.add_argument('-i', '--input_vcf', required=True, help='Path to input multisample vcf')
    parser.add_argument('-s', '--snp_threshold', required=False, default=4, help='Threshold value for number of SNPs per window for recombination')
    parser.add_argument('-w', '--window', required=False, default=2, help='Window size measured in kb to use for identifying recombinant snps')
    parser.add_argument('-r', '--ref_seq', required=True, help='Reference Sequence file path')
    parser.add_argument('-j', '--job', required=True, help='Job name to use in output file names')
        
    return parser.parse_args()




def main():
    args = parse_args()
    
    ref_path = args.ref_seq

    
    # Filters parsnp vcf and creates .csv with only singletons
    singletons_from_parsnp_vcf(args.input_vcf, args.job)  
    
    with open(f'singletons_{args.job}.csv', 'r') as singles, open(f'singletons_no_recombination_{args.job}.csv', 'w') as singles_no_recomb:
        
        # creates dict with all snps
        all_snps = create_snp_dict(singles)
        singles.seek(0) # reset the cursor (just in case)
        
        # identifies recombination using dict
        recomb_dict = identify_recombination(all_snps, args.snp_threshold, args.window)
        
        # creates a set from the recomb dict
        recomb_set = create_snp_sets(recomb_dict)
        
        # Filters the singletons to remove recombinant snps, writes to a new file
        write_filtered_snps_to_file(recomb_set, singles, singles_no_recomb)
        
    # contructs the ref seq from the fasta.ref file or FNA file
    ref_seq, ref_dict = read_refseq_file(ref_path)
    
    # returns a dict with the total number of each base found in the ref seq
    ref_bases = count_refseq_nucleotides(ref_seq)
    
    # returns a dict with all the triplets/3mers in the ref seq and the total number found in the ref seq for each
    ref_3mers = get_ref_triplets(ref_seq)
    

    # creates the dataframes that will be used to create visualisations
    create_sample_df(ref_dict, ref_bases, ref_3mers, args.job)

    
      
    
if __name__ == '__main__':
    main()
    
    
                
