import msprime
import sys
import subprocess
import shutil
import os
import math
import argparse
import random

########################
# Important parameters #
########################

parser = argparse.ArgumentParser(description='Simulate chromosomes that are compatible with gargammel.')

# Age of ancient individuals
parser.add_argument('--gens', nargs=1, help='Age of ancient individuals to be simulated', required=True)
# Number of bases to simulate
parser.add_argument('--bases', nargs=1, help='Number of base pairs to be simulated', required=True)
# Time of mass migration / population split
parser.add_argument('--split', nargs='?', const='25', default='0')

# Number of contamination individuals
contamination_individuals = 1
# Number of present day individuals, used for ref.fa
present_individuals = 1
# Number of ancient individuals
ancient_individuals = 1

args = parser.parse_args()
gens = int(args.gens[0])
num_bases = int(args.bases[0])
split = int(args.split[0])

#####################################
# Simulate chromosomes with msprime #
#####################################

# Simulate mass migration with the given split argument
pop_configs = [
	msprime.PopulationConfiguration(initial_size=1e4),
	msprime.PopulationConfiguration(initial_size=1e4)
]

dem_events = [
	msprime.MassMigration(time=split, source=0, destination=1, proportion=1.0)
]

# Create sample list
samples = []

for _ in range(contamination_individuals + present_individuals):
        samples.append(msprime.Sample(0, 0))
        samples.append(msprime.Sample(0, 0))

# If there was a population split, we'll use 2 different populations
if(split > 0):
	ancient_pop = 1
else:
	ancient_pop = 0

for _ in range(ancient_individuals):
	samples.append(msprime.Sample(ancient_pop, gens))
	samples.append(msprime.Sample(ancient_pop, gens))

# Run simulation and extract results
if(split > 1):
	tree_seq = msprime.simulate(
		samples=samples, recombination_rate=2e-8,
        	mutation_rate=2e-8, length=num_bases,
		population_configurations=pop_configs,
		demographic_events=dem_events)
else:
	tree_seq = msprime.simulate(
		samples=samples, recombination_rate=2e-8,
        	mutation_rate=2e-8, length=num_bases, Ne=1e4)

############################################
# Transform data for seq-gen compatibility #
############################################

tree_filepath = 'tree_data'
tree_seq.dump(tree_filepath)

# Get Newick format tree and partitions
newick_filepath = 'newick_tree'
newick_file = open(newick_filepath, 'w')
subprocess.run(['msp', 'newick', '--precision', '14', tree_filepath], stdout=newick_file)
newick_file.close()

# Get each tree's interval, this needs to be appended to the beginning
# of each Newick tree described in the file. Intervals are used by seq-gen
# to merge the multiple trees that result from recombination
intervals = []
for tree in tree_seq.trees():
    length = tree.get_length()
    intervals.append(int(length))

# Fix rounding error
diff = num_bases - sum(intervals)

if diff != 0:
    intervals[len(intervals) - 1] += diff

# Get number of partitions and add intervals
partitions = 0
added_intervals = []
with open(newick_filepath, 'r') as newick_file:
    
    for line, interval in zip(newick_file, intervals):
        added_intervals.append('[' + str(interval) + '] ' + line)
        partitions += 1

# Overwrite Newick file with added intervals
with open(newick_filepath, 'w') as newick_file:
    newick_file.writelines(added_intervals)

##############################
# Run seq-gen on Newick tree #
##############################

# Magic numbers obtained through trial and error
if gens <= 100:
    branch_scale = 80e-9
elif gens <= 200:
    branch_scale = 73e-9
else:
    branch_scale = 60e-9

seqgen_filepath = 'sequence_data'

with open(seqgen_filepath, 'w') as seqgen_file:
    subprocess.run(['seq-gen', '-mHKY', '-l' + str(num_bases),
        '-s' + str(branch_scale), '-p', str(partitions),
        newick_filepath], stdout=seqgen_file)

# Sort sequences, msprime does not output chromosomes in order.
# We will also remove the header, since this sequence file will be split
# into several files representing our individuals

# Small auxiliary function for sorting with key
def get_key(s):
    return int(s.split()[0])

chr_sequences = []
with open(seqgen_filepath, 'r') as seqgen_file:
    chr_sequences = seqgen_file.readlines()

# Remove header and sort
chr_sequences.pop(0)
chr_sequences.sort(key=get_key)

# Write only sequence to file
with open(seqgen_filepath, 'w') as seqgen_file:
    for line in chr_sequences:
        seqgen_file.write(line.split()[1] + '\n')

##############################################
# Split seq-gen output into individual files #
##############################################

# Directory constants
con_dir = 'cont/'
anc_dir = 'endo/'

# Create clean directories
if os.path.exists(con_dir):
    shutil.rmtree(con_dir)

if os.path.exists(anc_dir):
    shutil.rmtree(anc_dir)

os.makedirs(con_dir)
os.makedirs(anc_dir)

# Split individuals
with open(seqgen_filepath, 'r') as seqgen_file:
    chr_sequences = seqgen_file.readlines()

chr_index = 0
    
#####################
# Split contaminant #
#####################

for _ in range(contamination_individuals):

    # Individual string
    ind_string = 'cont'
        
    # Write two chromosomes
    for i in range(2):
        filename = ind_string + '.{}.fa'.format(i + 1)

        with open(con_dir + filename, 'w') as f:
            # Write header
            f.write('>cont_1\n')
            # Write sequence
            f.write(chr_sequences[chr_index])

        chr_index += 1

#################################
# Split present day individuals #
#################################

for _ in range(present_individuals):

    # Individual string
    ind_string = 'ref'
        
    # Write two chromosomes
    for i in range(2):
        filename = ind_string + '.{}.fa'.format(i + 1)

        with open(filename, 'w') as f:
            # Write header
            f.write('>reference_1\n')
            # Write sequence
            f.write(chr_sequences[chr_index])

        chr_index += 1

#############################
# Split ancient individuals # 
#############################

for _ in range(ancient_individuals):

    # Individual string
    ind_string = 'endo'
        
    # Write two chromosomes
    for i in range(2):
        filename = ind_string + '.{}.fa'.format(i + 1)

        with open(anc_dir + filename, 'w') as f:
            # Write header
            f.write('>ancient_1\n')
            # Write sequence
            f.write(chr_sequences[chr_index])

        chr_index += 1

##################
# Write segsites #
##################

# Individual string
ind_string = 'endo'

# Build both filepaths
chr_path_1 = anc_dir + ind_string + '.1.fa'
chr_path_2 = anc_dir + ind_string + '.2.fa'

# Segregating sites filepath
seg_path = anc_dir + 'segsites'

# Lines that will be written to segsites file
out = []

# Open both chromosome files
with open(chr_path_1, 'r') as chr1, open(chr_path_2, 'r') as chr2:
    # Extract sequences
    line_1 = chr1.readlines()[1]
    line_2 = chr2.readlines()[1]

    # Compare base by base
    site = 1
    for x, y in zip(line_1, line_2):

        if x != y:
            out.append('>ref_1\t' + str(site) + '\t' + x.upper() +
                    '\t' + y.upper() + '\n')

        site += 1

# Write to file
with open(seg_path, 'w') as f:
    f.writelines(out)

###################
# Merge reference #
###################

# Build both reference chromosome filepaths
chr_path_1 = 'ref.1.fa'
chr_path_2 = 'ref.2.fa'

# Segregating sites filepath
seg_path = 'ref.fa'

# Lines that will be written to segsites file
out = []
out.append('>ref_1\n')

# Open both chromosome files
with open(chr_path_1, 'r') as chr1, open(chr_path_2, 'r') as chr2:
    # Extract sequences
    line_1 = chr1.readlines()[1]
    line_2 = chr2.readlines()[1]

    # Compare base by base
    for x, y in zip(line_1, line_2):

        if x != y:
                if random.random() < 0.5:
                    out.append(x.upper())
                else:
                    out.append(y.upper())
        else:
            out.append(x.upper())

# Write to file
with open(seg_path, 'w') as f:
    f.writelines(out)

###########
# Cleanup #
###########

# Remove Newick, HDF5, and Seq-Gen files
if os.path.exists(newick_filepath):
    os.remove(newick_filepath)
if os.path.exists(tree_filepath):
    os.remove(tree_filepath)
if os.path.exists(seqgen_filepath):
    os.remove(seqgen_filepath)

# Remove unmerged reference chromosomes
os.remove('ref.1.fa')
os.remove('ref.2.fa')
