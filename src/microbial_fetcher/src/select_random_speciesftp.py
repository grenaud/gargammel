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


currentpwd = str(os.getcwd());
#sys.stderr.write("#"+currentpwd+"#");
#sys.exit(1);
logging.basicConfig(filename=currentpwd+"/Microbial_ID.log",
                    level=logging.INFO,
                    filemode="w",
                    format="%(message)s")

logging.info(logformat("#Microbial_Strain", "ID", "ID_asm", "FTP"))

sum_abundance = 0
species_lst = []
metanames = {}
metanames_path = os.path.join(os.path.dirname(__file__),
                              "species.ftp.metanames.txt.gz")


with gzip.open(metanames_path, "r") as fin:
    for line in fin:
        ID, ftp, species = re.split("\s+", line.rstrip("\n"))

        if species in metanames.keys():
            metanames[species].append((ID, ftp))
        else:
            metanames[species] = [(ID, ftp)]

for line in fileinput.input():
    species, abundance = re.split("\s+", line.rstrip())

    try:

        hits = metanames[species]
        out = randint(0, len(hits)-1)
        ID, ftp = hits[out]
        ID_asm = ftp.split("/")[-1]
        ftp_fullpath = ftp+"/"+ID_asm+"_genomic.fna.gz"

        # divide by 100 to get fraction instead of percentage
        logging.info(logformat(species, ID, ID_asm,
                               ftp_fullpath))
        sys.stdout.write(ftp_fullpath+"\n")
        sum_abundance += float(abundance)
        species_lst.append((ID_asm+"_genomic.fna", float(abundance)))
    except KeyError:
        sys.stderr.write("{} is not present"
                         " in the database\n".format(species))

listfile = ""+currentpwd+"/fasta/list";
#sys.stderr.write(listfile);

with open( listfile, "w") as fh_out:
    for species_name, abundance in species_lst:
        fh_out.write(listformat(species_name, abundance/sum_abundance))
