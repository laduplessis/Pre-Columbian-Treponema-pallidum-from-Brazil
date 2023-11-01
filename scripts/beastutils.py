import os, sys, io
import numpy as np
from subprocess import call
from Bio import Phylo, Nexus, SeqIO

def writeDateTrait(dates,  output):

	traits = []	
	for i in dates.keys():
		traits.append('\n\t\t\t\t\t%s=%.10f' %(i, dates[i]))
	output.write(','.join(traits)+'\n')
#

def writeBlankAlignment(newicktree, output):
	"""
		Writes alignment with empty sequences and taxa in the tree
	"""

	tree   = Nexus.Trees.Tree(newicktree)	
	taxa   = tree.get_taxa()

	for i in range(0,len(taxa)):
		output.write('\n\t\t\t\t\t<sequence spec="Sequence" taxon="%d" value="?" />' % (i+1))
#

def writeAlignment(filename, output, fasep="|"):
      """
      	Writes alignment in BEAST format to output buffer
      	Return dictionary with sequence ids and dates
      """

      datedict = dict()
      tipdates = dict()
      for seq in SeqIO.parse(filename, "fasta"):		

            # Process date
            datestring = seq.id.split(fasep)[-1]
            #(seqdate, sequpper, seqlower)  = getYearDate(datestring, datefmt = "%d-%m-%Y")            
            #print("%s\t%f\t%f\t%f\n" % (datestring, seqdate, sequpper, seqlower))

            #if (sequpper - seqlower != 0): 
            #    sys.stdout.write("Estimating sampling date for: %s\n" % seq.id)
            #    tipdates[seq.id] = [sequpper, seqlower]
            #

            datedict[seq.id] = float(datestring)

            # Write sequence
            output.write('\t\t\t<sequence id="%s:msa" taxon="%s" totalcount="4" value="%s"/>\n' % (seq.id, seq.id, seq.seq))

      #

      #return({'datetrait' : datedict, 'tipdates' : tipdates})
      return(datedict)
#



def getSamplingTimesFromTree(newicktree, forwards=False):
	"""
		Read sampling times from tree
		 
		If forwards = True  gives times as time from tMRCA
		If forwards = False gives times as time from most recent sample		
	"""

	# Sampling times from tree
	tree    = Nexus.Trees.Tree(newicktree)
	leaves  = tree.get_terminals()
	heights = np.zeros(len(leaves))	
	for i in range(0,len(leaves)):
		heights[i] = tree.sum_branchlength(node=leaves[i])

	# Forwards or backwards
	if (forwards == True):
		treetimes = heights
	else:
		treetimes = max(heights) - heights

	# Associate with label
	times = dict()
	for i in range(0,len(treetimes)):
		label = tree.get_taxa(node_id=leaves[i])[0]
		times[label] = treetimes[i]		
	#	

	return times	
#


def makeXMLFile(pars, template, outputfile="", outputpath=""):

	if (outputpath == ""):
		outputpath = pars['outputpath']

	if (not os.path.exists(outputpath)):
		os.makedirs(outputpath)

	if (outputfile == ""):
		outputfile = pars["name"]

	sys.stdout.write(outputfile+"...\n")

	formatpars = dict()
	for par in pars:
		formatpars['$'+par] = pars[par]
	output = template.format(**formatpars)

	outfile = open(outputpath+"/"+outputfile+".xml", 'w')
	outfile.write(output)
	outfile.close()
#


def formatPars(pars, inputdates=None, fasep="|"):

	# Set inputtree, sampling times and dummy alignment
	output_align = io.StringIO()
	output_dates = io.StringIO()
	output_model = io.StringIO()

	dates = writeAlignment(pars["alnfile"], output_align, fasep=fasep)	

	if (inputdates != None):	
		for seq in dates:
			seqid = seq[:seq.find(fasep)]
			dates[seq] = inputdates[seqid]

	#
	writeDateTrait(dates, output_dates)

	pars["alignment"] = output_align.getvalue()
	pars["dateTrait"] = output_dates.getvalue()

#


def get_beast_call(beast, filename, seed):	
	return ["java", "-jar", beast, "-seed", str(seed), "-overwrite", filename]