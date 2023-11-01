import os, sys, time
from optparse import OptionParser
from scipy.stats import norm
from Bio import SeqIO
#from Bio.Alphabet import generic_dna

################################################################################################################################
# 
# Create NEXUS alignment with sampling dates and tip calibration priors for sequences without exact sampling dates that can be 
# directly read into BEAUti2.
#
# Assumes sequence ids are well-behaved and all composed of the same number of blocks with a separator between them. 
#
# The input option dateidx is used to specify the index of the sampling date or date range (numbering from 0). The option fasep 
# specifies the character used to separate the blocks in the fasta sequence ids. There are 3 possible types of priors that can be 
# placed on sequences without exact sampling dates: Hard uniform priors, soft normal priors or wide uniform priors. Wide uniform 
# priors are not related to the upper and lower limts of the sampling dates, but can be used to test whether or not clock-like 
# evolution is consistent with the sampling date limits.
#
# The input option cladeidx can be used to specify the index of the clade a taxon belongs to. If specified, taxon sets are added
# for all clades (by default BEAUti2 makes these monophyletic and logs the tMRCAs).
#
# The output file has the same name as the input file, but with extension ".nexus" and saved to outputpath
#
################################################################################################################################
#
# Examples:
#
# If the sequence id format is: "<name>|<accession>|<clade>|<date>" OR "<name>|<accession>|<clade>|<date-lower>-<date-upper>"
# Then:
#
#     - python fasta2nexus.py -i <inputfile> -f "|" -d 3 
#       
#       If there is no range the sampling date is fixed to <date>
#       Otherwise a uniform prior between <date_lower> and <date_upper> is placed on the sampling date
#       and the starting date is set to the middle of the range.
#
#     - python fasta2nexus.py -i <inputfile> -f "|" -d 3 -n
#
#       If there is no range the sampling date is fixed to <date>
#       Otherwise a normal prior is placed with mean equal to the middle of the range and with standard 
#       deviation set such that the 2.5 quantile is at <date_lower> and the 97.5 quantile at <date_upper>
#       (i.e. the date range spans 95% of the density)
#       
#     - python fasta2nexus.py -i <inputfile> -f "|" -L 1000 -U 2016 -d 3 -w
#
#       If there is no range the sampling date is fixed to <date>
#       Otherwise <date_lower> and <date_upper> are ignored and a uniform prior between 1000 and 2016 is
#       placed on the sequence date (seq = uniform(1000,2016)).
#
#     - python fasta2nexus.py -i <inputfile> -f "|" -d 3 -p
#
#       All sampling dates are fixed to either <date> or the middle of the date range.
#
################################################################################################################################
# Parameters

usage = "usage: %prog [option]"
parser = OptionParser(usage=usage)

parser.add_option("-i","--inputfile",
                  dest = "inputfile",
                  default = "",
                  metavar = "path",
                  help = "Path to input csv file [required]")

parser.add_option("-o","--outputpath",
                  dest = "outputpath",
                  default = "",
                  metavar = "path",
                  help = "Path to save output file in [required]")

parser.add_option("-f","--fastasep",
                  dest = "fasep",
                  default = "|",
                  metavar = "string",
                  help = "Separator between fields in fasta sequence id's [default: %default] [optional]")

parser.add_option("-d","--dateidx",
                  dest = "dateidx",
                  default = "-1",
                  metavar = "integer",
                  help = "Index of the field representing the sampling date or date range in the sequence id [default: %default] [optional]")

parser.add_option("-c","--cladeidx",
                  dest = "cladeidx",
                  default = "-1",
                  metavar = "integer",
                  help = "Index of the field representing the clades. If specified monophyly constraints are added [default: %default] [optional]")

parser.add_option("-L","--lower",
                  dest = "lower",
                  default = "0",
                  metavar = "float",
                  help = "Lower limit for uniform prior to be placed on sequences with unknown sampling dates [default: %default] [optional]")

parser.add_option("-U","--upper",
                  dest = "upper",
                  default = "2000",
                  metavar = "float",
                  help = "Upper limit for uniform prior to be placed on sequences with unknown sampling dates  [default: %default] [optional]")

parser.add_option("-n","--normal",
                  dest = "normal",
                  default = False,
                  action = "store_true",
                  help = "Use normally distributed priors around the middle of the date range [default: %default] [optional]")

parser.add_option("-w","--wide",
                  dest = "wide",
                  default = False,
                  action = "store_true",
                  help = "Use uniformly distributed priors between upper and lower [default: %default] [optional]")

parser.add_option("-p","--point",
                  dest = "point",
                  default = False,
                  action = "store_true",
                  help = "Use point priors between upper and lower (fixed sampling dates) [default: %default] [optional]")

(options,args) = parser.parse_args()

inputfile   = os.path.abspath(options.inputfile)
outputpath  = os.path.abspath(options.outputpath)+"/"
fasep       = options.fasep
dateidx     = int(options.dateidx)
cladeidx    = int(options.cladeidx)
lower       = float(options.lower)
upper       = float(options.upper)
normal      = options.normal
wide        = options.wide
point       = options.point

################################################################################################################################

start = time.time()

if (dateidx < 0):
      sys.stdout.write("Error: Please specify indices of the sampling dates or date ranges\n")
      sys.exit()
      
if (wide + normal + point > 1):
      sys.stdout.write("Error: Only one of wide, normal or point can be true\n")
      sys.exit()


if (wide):     
      sys.stdout.write("Using uniform priors between %.2f and %.2f\n" % (lower, upper))
elif (normal):
      sys.stdout.write("Using normal priors covering date range for each sequence\n")
elif (point):
      sys.stdout.write("Using point priors (fixed sampling dates) in the middle of the date range for each sequence\n")
else:
      sys.stdout.write("Using uniform priors across specified date range for each sequence\n")


if (not os.path.exists(outputpath)):
    os.mkdir(outputpath)
outputfile = outputpath+inputfile[inputfile.rfind('/')+1:inputfile.rfind(".")]+".nexus"

# Write sequences to NEXUS (switch to generic_alphabet)
#seqs = list(SeqIO.parse(inputfile, "fasta", alphabet=Alphabet.generic_dna))
#SeqIO.write(seqs, outputfile, "nexus")
seqs = list(SeqIO.parse(inputfile, "fasta"))
SeqIO.convert(inputfile, "fasta", outputfile, "nexus", "DNA")


# Get calibrations
calibrations = []
clades       = dict()
for seq in seqs: 
      idparts  = seq.id.split(fasep)

      # Parse tip dates or calibrations
      seqdates = list(map(float, idparts[dateidx].strip().split("-")))

      if (len(seqdates) == 1):
            seqdate = seqdates[0]
            sequpper = seqdate
            seqlower = seqdate
      elif (len(seqdates) == 2):            
            (seqlower, sequpper) = seqdates
            seqdate = 0.5*(seqlower + sequpper)
      else:
            sys.stdout.write("Error: Too many or too few dates in date range\n")

      # Output tip dates or calibrations
      if (seqlower == sequpper):
            sys.stdout.write("%20s = %d\n" % (idparts[0], seqdate))
            outstr = "%s = fixed(%d)" % (seq.id, seqdate)
      else:            
            sys.stdout.write("%20s = %d (%d-%d)\n" % (idparts[0], seqdate, seqlower, sequpper))

            if (normal):
                  m = 0.5*(sequpper + seqlower)
                  s = (sequpper - m)/norm.ppf(0.975)

                  outstr = "%s = normal(%.5f,%.5f)" % (seq.id, m, s)      

            elif (wide):
                  outstr = "%s = uniform(%.5f,%.5f)" % (seq.id, lower, upper)      

            elif (point):
                  outstr = "%s = fixed(%.5f)" % (seq.id, 0.5*(seqlower + sequpper))

            else:
                  outstr = "%s = uniform(%.5f,%.5f)" % (seq.id, seqlower, sequpper)      

      calibrations.append("\tcalibrate %s" % outstr)


      # Parse clades
      if (cladeidx >=0):
            clade = idparts[cladeidx].strip()
            if (clade in clades.keys()):
                  clades[clade].append(seq.id)
            else:
                  clades[clade] = [seq.id]




# Add NEXUS block with dates and calibrations
nex = open(outputfile, "a")
nex.write("\n\nbegin assumptions;\n\toptions scale = years;\n\n")
nex.write(",\n".join(calibrations))
nex.write(";\n")
if (cladeidx >= 0):

      for clade in clades.keys():
            nex.write("\n\ttaxset %s = %s;" % (clade, " ".join(clades[clade])))

nex.write("\nend;\n\n")


end = time.time()
sys.stdout.write("Total time taken: "+str(end-start)+" seconds\n")    