#!/usr/bin/env python
import fileinput
import sys
import re
import os
import gzip
import logging
from random import randint
logformat = "{}\t{}\t{}\t{}".format
listformat = "{}\t{}\n".format
species_filter = "{}\t{} below abundance cut off: {}\n".format

currentpwd = os.getcwd()
logging.basicConfig(filename=currentpwd+"/Microbial_ID.log",
                    level=logging.INFO,
                    filemode="w",
                    format="%(message)s")

logging.info(logformat("#Microbial_Strain", "ID", "ID_asm", "FTP"))

sum_abundance = 0
species_lst = []
metanames = {}
metanames_path = os.path.join(os.path.dirname(__file__),"species.ftp.metanames.txt.gz")


with gzip.open(metanames_path, "r") as fin:
    for line in fin:
        ID, ftp, species = re.split("\s+", line.decode().rstrip("\n"))
        if species in metanames.keys():
            metanames[species].append((ID, ftp))
        else:
            metanames[species] = [(ID, ftp)]

try:
    abundance_cut_off = float(sys.argv[1])
except IndexError:
    abundance_cut_off = 0.1  # in percentages (1/1000) min relative abundance
count_species_filter = 0

for line in fileinput.input():
    species, abundance = re.split("\s+", line.decode().rstrip())
    if float(abundance) < abundance_cut_off:
        count_species_filter += 1
        # sys.stderr.write(species_filter(species, abundance,
        #                                 abundance_cut_off))
        continue
    try:
        hits = metanames[species]
        out = randint(0, len(hits)-1)
        ID, ftp = hits[out]
        ID_asm = ftp.split("/")[-1]
        ftp_fullpath = ftp+"/"+ID_asm+"_genomic.fna.gz"
        #print("ftp_fullpath");
        #print(ftp_fullpath);
        #sys.exit(1);
        logging.info(logformat(species, ID, ID_asm,ftp_fullpath))
        sys.stdout.write(ftp_fullpath+"\n")
        sum_abundance += float(abundance)
        species_lst.append((ID_asm+"_genomic.fna", float(abundance)))
    except KeyError:
        sys.stderr.write("{} is not present in the database\n".format(species))

sys.stderr.write(str(len(species_lst)) + " genomes will be downloaded from ncbi.\n" + str(count_species_filter) + " microbe abundances were below" + " the minimum relative abundance cut off on " + str(abundance_cut_off) + "%\n")

listfile = currentpwd+"/fasta/list"

with open(listfile, "w") as fh_out:
    for species_name, abundance in species_lst:
        fh_out.write(listformat(species_name, abundance/sum_abundance))
