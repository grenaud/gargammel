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
demo = msprime.Demography()
demo.add_population(initial_size=1e4, name="ancient")
demo.add_population(initial_size=1e4, name="modern")

demo.add_mass_migration(time=split, source="modern", dest="ancient", proportion=1.0)

# Create sample list
samples = []

if split > 0:
    samples.append(msprime.SampleSet(contamination_individuals + present_individuals, "modern", time=0))
else:
    samples.append(msprime.SampleSet(contamination_individuals + present_individuals, "ancient", time=0))
    
samples.append(msprime.SampleSet(ancient_individuals, "ancient", time=gens))

# Run simulation and extract results

tree_seq = msprime.sim_ancestry(samples=samples, demography=demo,
    	                        sequence_length=num_bases, recombination_rate=2e-8)

############################################
# Transform data for seq-gen compatibility #
############################################

# Get Newick format tree and partitions
newick_filepath = 'newick_tree'

# Get each tree's interval, this needs to be appended to the beginning
# of each Newick tree described in the file. Intervals are used by seq-gen
# to merge the multiple trees that result from recombination

# Get number of partitions and add intervals
partitions = tree_seq.num_trees

with open(newick_filepath, 'w') as newick_file:
    for tree in tree_seq.trees():
        length = str(int(tree.get_length()))
        newick_file.write('[' + length + '] ' + tree.newick() + "\n")

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

# Simulate sequences given ancestry trees

with open(seqgen_filepath, 'w') as seqgen_file:
    subprocess.run(['seq-gen', '-mHKY', '-l' + str(num_bases),
        '-s' + str(branch_scale), '-p', str(partitions),
        newick_filepath], stdout=seqgen_file)

with open(seqgen_filepath, 'r') as seqgen_file:
    chr_sequences = seqgen_file.readlines()

# Remove header and create {id => sequence} dictionary

chr_sequences.pop(0)
chr_dict = dict((int(chr_seq.split()[0]), chr_seq.split()[1]) for chr_seq in chr_sequences)

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

chr_index = 1
    
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
            f.write(chr_dict[chr_index])

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
            f.write(chr_dict[chr_index])

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
            f.write(chr_dict[chr_index])

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
if os.path.exists(seqgen_filepath):
    os.remove(seqgen_filepath)

# Remove unmerged reference chromosomes
os.remove('ref.1.fa')
os.remove('ref.2.fa')
