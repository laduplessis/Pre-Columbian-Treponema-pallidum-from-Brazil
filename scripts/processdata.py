import os, sys, re, time, csv
from optparse import OptionParser
from Bio import Seq, SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord


# Process data
################################################################################################################################
# Pre-process datafiles
#
# Input:  Database (csv) file and alignment (fasta) file
# Output: Database (csv) file, alignment (fasta - .fas) file, unaligned sequence (fasta - .fa) file, 
#         list of skipped sequences, list of not found sequences
#
# Options:
#   - Change dbsep or fasep to modify the separators in the output database and fasta file sequence ids
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

parser.add_option("-a","--alignfile",
                  dest = "alignfile",
                  default = "",
                  metavar = "path",
                  help = "Path to input alignment file [required]")

parser.add_option("-o","--outputpath",
                  dest = "outputpath",
                  default = "",
                  metavar = "path",
                  help = "Path to save output file in [required]")

parser.add_option("-p","--prefix",
                  dest = "prefix",
                  default = "output",
                  metavar = "string",
                  help = "Prefix for the output file(s) [default: %default] [required]")

parser.add_option("-d","--dbsep",
                  dest = "dbsep",
                  default = ",",
                  metavar = "string",
                  help = "Separator between columns in output database [default: %default] [optional]")

parser.add_option("-f","--fastasep",
                  dest = "fasep",
                  default = "|",
                  metavar = "string",
                  help = "Separator between fields in fasta sequence id's [default: %default] [optional]")

parser.add_option("-s","--seqidfields",
                  dest = "seqidfields",
                  default = "",
                  metavar = "comma-separated list of strings",
                  help = "Comma-separated list of fields in the input csv file to use in the output sequence ids, by default "
                         "uses the ids from the input alignment [optional]")

parser.add_option("-c","--clades",
                  dest = "clades",
                  default = "",
                  metavar = "comma-separated list of strings",
                  help = "Comma-separated list of clades (Subspecies column) to extract. Leave blank to extract all. [optional]")

parser.add_option("-e","--exclude",
                  dest = "exclude",
                  default = "",
                  metavar = "comma-separated list of strings/path",
                  help = "Comma-separated list of sequence names/accessions to exclude or csv file with sequence ids in first column [optional]")



(options,args) = parser.parse_args()

inputfile   = os.path.abspath(options.inputfile)
alignfile   = os.path.abspath(options.alignfile) if options.alignfile != "" else ""
outputpath  = os.path.abspath(options.outputpath)+"/"
prefix      = options.prefix
dbsep       = options.dbsep
fasep       = options.fasep
seqidfields = [] if options.seqidfields == "" else options.seqidfields.split(",")
clades      = [] if options.clades == "" else options.clades.split(",")
exclude     = [] if options.exclude == "" else options.exclude.split(",")


################################################################################################################################

def readAlignment(inputfile, fmt, sep="|"):
      """Read alignment from file

      """

      alignfile = open(inputfile,"r")
      alignment = AlignIO.read(alignfile,fmt)
      alignfile.close()

      seqdict = dict()
      for seqrecord in alignment:
          i = seqrecord.id.find(sep)
          if (i > 0):
              seqid = seqrecord.id[:]
          else:
              seqid = seqrecord.id
          
          seqdict[seqid] = seqrecord.seq

      return seqdict
#

def getMapping(header, fields):
    """Figure out the mapping between the database columns and columns we are processing

    """ 
    
    mapping = dict()
    for field in fields:        
        for i in range(0,len(header)):                    
            if (header[i].strip() == field):
                mapping[field] = i

    return mapping
#


################################################################################################################################


start = time.time()

if (not os.path.exists(outputpath)):
    os.mkdir(outputpath)


# Read alignment sequences
if (alignfile != ""):
    msaseqs   = readAlignment(alignfile, "fasta", fasep)
    msaseqids = list(msaseqs.keys())
else:
    sys.stdout.write("WARNING: No alignment file specified, exiting...\n")
    sys.exit()

# Get column mapping
db        = csv.reader(open(inputfile,'r',  encoding="utf-8"))
header    = db.__next__()


if (len(seqidfields) == 0):
    seqidfields = header

mapping = getMapping(header, seqidfields)
sys.stdout.write("\nNEW MAPPING:\n")
for field in seqidfields:
    sys.stdout.write("%15s\t: " % field)
    if (field in mapping.keys()):
        sys.stdout.write(header[mapping[field]]+" (%d)\n" % mapping[field])
    else: 
        sys.stdout.write("NA\n")


i = 0
sequences = []
skipped   = []
notfound  = []
outputdb  = open(outputpath+prefix+".csv", "w")
outputdb.write(dbsep.join(seqidfields)+"\n")
for line in db:
    if (len(line) > 0):

        i += 1

        # Map and process data from database
        row = dict()
        for field in mapping.keys():
           
            # Read field from mapping and remove spaces
            data = line[mapping[field]].strip()

            # Change all except sequence id and accession to upper case
            if (field != "name" and field != "accession" and field != "Taxon"):
                data = data.upper()              

            row[field] = data


        # Exclude specific sequences
        if (row["Taxon"] in exclude):
            sys.stdout.write("Skipping %s (in exclude list)...\n" % row["Taxon"])
            skipped.append(dbsep.join([row[i] for i in seqidfields]) + ",in exclude list\n")
            continue


        # Select clades (Subspecies) to include
        if (len(clades) > 0 and row["Subspecies"] not in clades):
            sys.stdout.write("Skipping %s (%s clade not included)...\n" % (row["Taxon"], row["Subspecies"]))
            skipped.append(dbsep.join([row[i] for i in seqidfields]) + ",%s clade not included\n" % row["Subspecies"])        
            continue

    
        # Get and check the sequence for the fasta file
        if (alignfile != ""):

            # Need to search in two ways since different fasta files have different sequence id conventions
            seqid = row["Taxon"]
            
            # Found the sequence - continue
            if (seqid in msaseqids):                
                msaseq = msaseqs[seqid]                
                
                # We have seen the sequence
                msaseqids.remove(seqid)

            else:
                sys.stdout.write("ERROR! %s not found in alignment file!\n" % row["Taxon"])
                skipped.append(dbsep.join([row[i] for i in seqidfields]) + ",not in alignment\n")
                notfound.append(seqid)
                continue        
        

        # New sequence id for fasta file  
        newseqid = [] 
        for field in seqidfields: 
            newseqid.append(row[field])

        #newseqid = fasep.join([row["Taxon"], row["accession"], row["date"]])        
        sequences.append(SeqRecord(msaseq, id=fasep.join(newseqid), description=""))        

        # Write new row to new db file
        outputdb.write(dbsep.join([row[i] for i in seqidfields])+"\n") 
#
outputdb.close()

# Check if sequences are aligned
aligned = True
seqlen = len(sequences[0].seq)
for seq in sequences:
    if (len(seq.seq) != seqlen):
        aligned = False
#
if (aligned):
    SeqIO.write(sequences,outputpath+prefix+".fas","fasta")
    #SeqIO.write(sequences,outputpath+prefix+".stockholm", "stockholm")
else:
    sys.stdout.write("WARNING: Sequences not aligned! Alignment will not be written...\n")
#SeqIO.write(sequences_ungapped,outputpath+prefix+".fa","fasta")

# Write skipped sequences
if (len(skipped) > 0): 
    with open(outputpath+prefix+".skipped.csv","w") as f:
        f.write(dbsep.join(seqidfields)+",Reason\n")
        for skip in skipped:
            f.write(skip)

# Write sequences not found
if (alignfile != "" and len(notfound) > 0):
    with open(outputpath+prefix+".notfound.csv","w") as f:
        f.write("\n".join(notfound)+"\n")


# Print some statistics
sys.stdout.write("\nRESULTS SUMMARY\n")
sys.stdout.write("%10d sequence(s) evaluated\n" % i)
sys.stdout.write("%10d sequence(s) successfully extracted\n" % len(sequences))
sys.stdout.write("%10d sequence(s) skipped\n" % len(skipped))
if (alignfile != ""): 
    sys.stdout.write("%10d sequence(s) in the alignment not found in the database\n" % len(msaseqids))

end = time.time()
sys.stdout.write("Total time taken: "+str(end-start)+" seconds\n")
