#!/usr/bin/python
# Get_EVS_Data.py - This will take in the EVS_TXT file and a file containing positions of interest (various formats), and spit out EVS data for ONLY positions in the position file.
# Usage: python Get_EVS_Data.py /PATH/TO/Var_File /PATH/TO/EVS/TXTs /PATH/TO/Output_Dir
       # python ${GET_EVS_DATA_PY} ${VAR_FILE} ${EVS_DIR} ${OUTPUT_DIR}
# Output: File with EVS data for variants contained in a variant file

# By Regina Nguyen and Mollee Jain (TSRI - July 2013)
# Modified by Tristan Carland (TSRI - Jan 2014)
# Modified by Kristopher Standish (UCSD/JCVI - March/April 2014)

from sys import argv # standard packages to allow imput params
from sys import exit
import gzip
import datetime

# Checks for the correct number of input params, displays usage info if needed
try:
    script, var_file, evs_dir, out_dir, = argv
except ValueError:
    print "\nScript used to create subsets of variants from vcf files that correspond to each of the snps lists."
    exit("cmd_usage$> script  var_file  evs_dir  out_dir\n")
# Input params:
# script - this script
# var_file - whatever file has the list of variants...could be one of multiple formats
# evs_dir - directory in which the EVS text files are in
# out_dir - directory to which the output file will be written

# Make this "GetVar" command that will be handy...helps with indexing
getVar = lambda searchList, ind: [searchList[i] for i in ind]

CHROMS = range(2,23)
CHROMS.append("X")
CHROMS.append("Y")
print CHROMS

###########################################################
## Open up files and get things queued up

# Open Variant File
if "gz" in var_file:
  var = gzip.open(var_file)
if "gz" not in var_file:
  var = open(var_file)
print "Variant File Open - %s" % datetime.datetime.now().time().isoformat()
var_line = var.next()

# Skip all the comment lines up top on the Variant file
while var_line.startswith("#"):
  var_line = var.next()
print "Comments Skipped in Variant File - %s" % datetime.datetime.now().time().isoformat()

# Open up Chr1 EVS File
EVS_Path_Shell = "%s/ESP6500SI-V2-SSA137.updatedProteinHgvs.chrXYZ.snps_indels.txt" % evs_dir
print EVS_Path_Shell
EVS_Path = EVS_Path_Shell.replace("XYZ","1")
print EVS_Path
evs_file = open(EVS_Path)
evs_line = evs_file.next()
while evs_line.startswith("##"):
  evs_line = evs_file.next()
print "Comments Skipped in EVS File - %s" % datetime.datetime.now().time().isoformat()

# Open Files to Write to
Out_E = "%s/EVS_Data.txt" % out_dir
OUT_E = open(Out_E, 'w')

Out_P = "%s/Pos_Data.txt" % out_dir
OUT_P = open(Out_P, 'w')

# Write header line of EVS file to output
OUT_E.write(evs_line)

###########################################################
## Make Dictionary for Variant File
print "Making Library for Var File - %s" % datetime.datetime.now().time().isoformat()
LOCS = {}
RSID = {}
if "vcf" in var_file:
  for line in var:
    temp = line.split("\t", 3)
    key = "%s:%s" % (temp[0],temp[1])
    LOCS[key] = 1
    key = "%s" % temp[2]
    RSID[key] = 1
if "bim" in var_file:
  for line in var:
    temp = line.split("\t", 4)
    key = "%s:%s" % (temp[0],temp[3])
    LOCS[key] = 1
    key = "%s" % temp[1]
    RSID[key] = 1
print "DONE Making Library for Var File - %s" % datetime.datetime.now().time().isoformat()


###########################################################
## Go through EVS file and pull out the variants in dictionary
for chrom in CHROMS:
  print "Going through EVS File for Chrom %s - %s" % (chrom, datetime.datetime.now().time().isoformat() )
  for evs_line in evs_file:
    temp = evs_line.split(" ", 2)
    evs_loc = temp[0]
    evs_rsid = temp[1]
    if ( evs_loc in LOCS or evs_rsid in RSID ):
      OUT_E.write(evs_line)
      OUT_P.write("%s\t%s\n" % (evs_loc, evs_rsid) )
  #    print evs_line
  EVS_Path = EVS_Path_Shell.replace("XYZ","%s" % chrom)
  evs_file = open(EVS_Path)
  evs_line = evs_file.next()
  while evs_line.startswith("#"):
    evs_line = evs_file.next()
  print "Loaded/Skipped Comments in EVS File for Chrom %s - %s" % (chrom, datetime.datetime.now().time().isoformat() )
#  print evs_line

OUT_E.close()
OUT_P.close()
evs_file.close()
var.close()
exit()