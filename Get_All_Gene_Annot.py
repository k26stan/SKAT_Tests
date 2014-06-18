#!/usr/bin/python
# Get_All_Gene_Annots.py - This will take in a large annotation file from Cypher and a file containing position (ranges) of interest (chr:pos1-pos2), and spit out annotations for ONLY positions in the position file.
  # Also spit out position files making use of the annotations (Coding, Damaging, etc...)
# Usage: python Get_All_Gene_Annot.py /PATH/TO/Gene_Coords.txt /PATH/TO/ANNOT.file /PATH/TO/BASE_DIR
       # python ${GET_ALL_GENE_ANNOT_PY} ${GENE_COORDS} ${ANNOTS} ${ASSOC}

# By Regina Nguyen and Mollee Jain (TSRI - July 2013)
# Modified by Tristan Carland (TSRI - Jan 2014)
# Modified by Kristopher Standish (UCSD/JCVI - March 2014)

from sys import argv # standard packages to allow imput params
from sys import exit
import gzip
import datetime

# Checks for the correct number of input params, displays usage info if needed
try:
    script, gene_coords, annot_file, base_dir, = argv
except ValueError:
    print "\nScript used to create subsets of variants from vcf files that correspond to each of the snps lists."
    exit("cmd_usage$> script  gene_coords  annot_file  base_dir\n")
# Input params:
# script - this script
# gene_coords - file that contains gene names. Will be used to search/open additional files
# annot_file - file with annotations at a bunch of variant positions
# base_dir - base directory containing all the gene directories

# Make this "GetVar" command that will be handy...helps with indexing
getVar = lambda searchList, ind: [searchList[i] for i in ind]

###########################################################
## Open up file and get things queued up

# open gene_list file
GENE_COORDS = open(gene_coords)
print "Gene List open - %s" % datetime.datetime.now().time().isoformat()

# open Annotation file
if "gz" in annot_file:
  annot = gzip.open(annot_file)
if "gz" not in annot_file:
  annot = open(annot_file)
print "Annotation File Open - %s" % datetime.datetime.now().time().isoformat()

# Get columns of interest in Annotation File
an_line = annot.next()
HEADER_LINE = an_line
Head_Out = "%s/Header.txt" % base_dir
HEAD_OUT = open(Head_Out, 'w')
HEAD_OUT.write(HEADER_LINE)
HEAD_OUT.close()

print HEADER_LINE
HEADERS = an_line.rstrip().split("\t")
CHR_COL = HEADERS.index("Chromosome")
CRD_COL = HEADERS.index("End")
GEN_COL = HEADERS.index("Gene")
LOC_COL = HEADERS.index("Location")
COD_COL = HEADERS.index("Coding_Impact")

# Skip all the comment lines up top on the Annotation file
while an_line.startswith("#"):
  an_line = annot.next()
print "Comments Skipped - %s" % datetime.datetime.now().time().isoformat()

###########################################################
## Set up 2 Genes at a time for which to pull out variants
GENE_COORDS.next() # Skip header line
line_count = 0
for g_line in GENE_COORDS:
  line_count += 1
GENE_COORDS.seek(0) # Back to first line
print line_count
GENE_COORDS.next() # Skip header line
# Name, Chrom, Start, Stop
GOI1 = GENE_COORDS.next().rstrip().split("\t")
GOI2 = GENE_COORDS.next().rstrip().split("\t")
GOI3 = GENE_COORDS.next().rstrip().split("\t")
GOI4 = GENE_COORDS.next().rstrip().split("\t")
print "Initial Genes Loaded"
print GOI1
print GOI2
print GOI3
print GOI4

# Open Files to Write Annotations to
An_Out_1 = "%s/%s/%s.Annot.txt" % (base_dir, GOI1[0], GOI1[0])
AN_OUT_1 = open(An_Out_1, 'w')
AN_OUT_1.write(HEADER_LINE)
An_Out_2 = "%s/%s/%s.Annot.txt" % (base_dir, GOI2[0], GOI2[0])
AN_OUT_2 = open(An_Out_2, 'w')
AN_OUT_2.write(HEADER_LINE)
An_Out_3 = "%s/%s/%s.Annot.txt" % (base_dir, GOI3[0], GOI3[0])
AN_OUT_3 = open(An_Out_3, 'w')
AN_OUT_3.write(HEADER_LINE)
An_Out_4 = "%s/%s/%s.Annot.txt" % (base_dir, GOI4[0], GOI4[0])
AN_OUT_4 = open(An_Out_4, 'w')
AN_OUT_4.write(HEADER_LINE)
print "First Gene Annotation Files Opened"

###########################################################
## Go through Annotation File and search for Genes/Vars
CHROM = str(1)
# for an_line in annot:
#   temp = an_line.split("\t")[0:COD_COL]
#   if temp[CHR_COL][3:]!=CHROM:
#     CHROM = temp[CHR_COL][3:]
#     print temp[CHR_COL][3:]
#     print CHROM

# exit()
which_line = 4
DONE_1 = 0
DONE_2 = 0
DONE_3 = 0
DONE_4 = 0
for an_line in annot:
  temp = an_line.split("\t",COD_COL)[0:COD_COL]
  # Get out info for GOI1
  if ( temp[CHR_COL][3:]==GOI1[1] and DONE_1==0 ):
    if ( int(temp[CRD_COL]) > int(GOI1[2]) and int(temp[CRD_COL]) < int(GOI1[3]) ):
      AN_OUT_1.write(an_line)
    if int(temp[CRD_COL]) > int(GOI1[3]):
      print "Done with %s" % GOI1[0]
      AN_OUT_1.close()
      if which_line < line_count:
        which_line += 1
        next_line = GENE_COORDS.next()
        print len(next_line)
        GOI1 = next_line.rstrip().split("\t")
        An_Out_1 = "%s/%s/%s.Annot.txt" % (base_dir, GOI1[0], GOI1[0])
        AN_OUT_1 = open(An_Out_1, 'w')
        AN_OUT_1.write(HEADER_LINE)
        print "On to %s" % GOI1[0]
      else:
        DONE_1 = 1
  # Get out info for GOI2
  if (temp[CHR_COL][3:]==GOI2[1] and DONE_2==0 ):
    if ( int(temp[CRD_COL]) > int(GOI2[2]) and int(temp[CRD_COL]) < int(GOI2[3]) ):
      AN_OUT_2.write(an_line)
    if int(temp[CRD_COL]) > int(GOI2[3]):
      print "Done with %s" % GOI2[0]
      AN_OUT_2.close()
      if which_line < line_count:
        which_line += 1
        next_line = GENE_COORDS.next()
        print len(next_line)
        GOI2 = next_line.rstrip().split("\t")
        An_Out_2 = "%s/%s/%s.Annot.txt" % (base_dir, GOI2[0], GOI2[0])
        AN_OUT_2 = open(An_Out_2, 'w')
        AN_OUT_2.write(HEADER_LINE)
        print "On to %s" % GOI2[0]
      else:
        DONE_2 = 1
  # Get out info for GOI3
  if (temp[CHR_COL][3:]==GOI3[1] and DONE_3==0 ):
    if ( int(temp[CRD_COL]) > int(GOI3[2]) and int(temp[CRD_COL]) < int(GOI3[3]) ):
      AN_OUT_3.write(an_line)
    if int(temp[CRD_COL]) > int(GOI3[3]):
      print "Done with %s" % GOI3[0]
      AN_OUT_3.close()
      if which_line < line_count:
        which_line += 1
        next_line = GENE_COORDS.next()
        print len(next_line)
        GOI3 = next_line.rstrip().split("\t")
        An_Out_3 = "%s/%s/%s.Annot.txt" % (base_dir, GOI3[0], GOI3[0])
        AN_OUT_3 = open(An_Out_3, 'w')
        AN_OUT_3.write(HEADER_LINE)
        print "On to %s" % GOI3[0]
      else:
        DONE_3 = 1
  # Get out info for GOI4
  if (temp[CHR_COL][3:]==GOI4[1] and DONE_4==0 ):
    if ( int(temp[CRD_COL]) > int(GOI4[2]) and int(temp[CRD_COL]) < int(GOI4[3]) ):
      AN_OUT_4.write(an_line)
    if int(temp[CRD_COL]) > int(GOI4[3]):
      print "Done with %s" % GOI4[0]
      AN_OUT_4.close()
      if which_line < line_count:
        which_line += 1
        next_line = GENE_COORDS.next()
        print len(next_line)
        GOI4 = next_line.rstrip().split("\t")
        An_Out_4 = "%s/%s/%s.Annot.txt" % (base_dir, GOI4[0], GOI4[0])
        AN_OUT_4 = open(An_Out_4, 'w')
        AN_OUT_4.write(HEADER_LINE)
        print "On to %s" % GOI4[0]
      else:
        DONE_4 = 1

  if ( DONE_1==1 and DONE_2==1 and DONE_3==1 and DONE_4==1 ):
    break

######### NEEDS WORK ##################

AN_OUT_1.close()
AN_OUT_2.close()
AN_OUT_3.close()
AN_OUT_4.close()
exit()