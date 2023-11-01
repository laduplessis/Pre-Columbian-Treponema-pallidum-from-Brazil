import os, sys, io, yaml, random
import numpy as np
from fnmatch import fnmatch
from operator import methodcaller
from subprocess import call
from optparse import OptionParser
from Bio import Phylo, Nexus

from beastutils import *

usage = "usage: %prog [option]"
parser = OptionParser(usage=usage)

parser.add_option("-i","--inputpath",
                  dest = "inputpath",
                  default = "",
                  metavar = "path",
                  help = "Path to input files [required]")

parser.add_option("-c","--config",
                  dest = "config",
                  default = "*.cfg",
                  metavar = "",
                  help = "Pattern to match for config files [default = %default]")

parser.add_option("-x","--template",
                  dest = "template",
                  default = "",
                  metavar = "path",
                  help = "Path to template XML file [required]")

parser.add_option("-o","--outputpath",
                  dest = "outputpath",
                  default = "",
                  metavar = "path",
                  help = "Path to save output file in [required]")

parser.add_option("-d","--dates",
                  dest = "dates",
                  default = "",
                  metavar = "path",
                  help = "csv file with sequence ids in first line, clades in second and dates in subsequent lines [required]")

parser.add_option("-n","--name",
                  dest = "name",
                  default = "",
                  metavar = "path",
                  help = "Name of the runs [required]")

parser.add_option("-s","--seed",
                  dest = "seed",
                  default = "127",
                  metavar = "integer",
                  help = "Seed or comma separated list of seeds to use [required]")

(options,args) = parser.parse_args()

if (options.inputpath != ""):
	config         = options.config
	inputpath      = os.path.abspath(options.inputpath)+"/"
else:
	config         = options.config[options.config.rfind("/")+1:]
	inputpath      = os.path.abspath(options.config[:options.config.rfind("/")])+"/"


################################################################################################################################  


def getDateTraitBlock(taxa, dates): 

	traits = []	

	for i in range(0,len(taxa)):
		traits.append('\t\t\t\t\t%s=%s' % (taxa[i], dates[i]))

	return(',\n'.join(traits))
#


################################################################################################################################  

for filename in sorted(os.listdir(inputpath)):
	if (fnmatch(filename,config)):

		sys.stdout.write(filename+"\t"+config+"\n")

		# Load config file
		configfile = open(inputpath+filename, 'r').read().replace("\t"," ")
		pars 	   = yaml.load(configfile, Loader=yaml.FullLoader)
	
		# Set BEAST specific parameters	
		seeds      = list(map(int, options.seed.split(',')))				
		basename   = pars["name"] if options.name == '' else pars["name"]+"_"+options.name
		outputpath = os.path.abspath(pars["outputpath"]    if options.outputpath == '' else options.outputpath)
		template   = open(os.path.abspath(pars["template"] if options.template == '' else options.template), 'r').read()		
		datesFile  = os.path.abspath(pars["dates"]    if options.dates == "" else options.dates)


		# Get dates and taxa
		dates_raw   = open(datesFile,'r').readlines()
		date_tuples = list(map(methodcaller("split", ","), dates_raw))

		taxa  = list(map( lambda x : x[0].strip(), date_tuples[1:]))
		dates = list(map( lambda x : float(x[1].strip()), date_tuples[1:]))

		# Output scripts
		if (not os.path.exists(outputpath)):
			os.makedirs(outputpath)
		scriptfile = open(outputpath+"/"+basename+".sh",'w')

		# Euler script header
		scriptfile.write("#!/bin/bash\n\nmodule load java phylo\n\n\n")
		#scriptfile.write("#SBATCH --ntasks=%d\n")
		#scriptfile.write("#SBATCH --time=%s\n")
		#scriptfile.write('#SBATCH --job-name="%s"\n')
		#scriptfile.write('#SBATCH --output="%s.out"\n')
		#scriptfile.write('#SBATCH --error="%s.err"")\n')
		#scriptfile.write("\n\nmodule load java phylo\n\n\n")

		#formatPars(pars)
		#makeXMLFile(pars, template, outputfile=pars["name"], outputpath=outputpath)

		# Create XML file with unshuffled dates
		pars["name"] = basename+".truth"			
		pars['dateTrait'] = getDateTraitBlock(taxa, dates)
		makeXMLFile(pars, template, outputfile=pars["name"], outputpath=outputpath)

		# Write command to scripts
		for seed in seeds:
			beastcmd = '%s -seed %s -overwrite -threads 4 %s.xml' % (pars["beast"], pars["seed"], pars["name"])
			eulercmd = 'sbatch --time=%s:00:00 --ntasks=%s --job-name="%s" --output="%s.out" --error="%s.err" --wrap="%s"' % (pars["hours"], pars["ntasks"], pars["name"], pars["name"], pars["name"], beastcmd)
			scriptfile.write("%s\n" % (eulercmd))


		# Replicates with shuffled dates
		for i in range(0,pars["replicates"]):
			
			# Shuffle dates
			random.shuffle(dates)

			# Create XML file			
			pars["name"] = basename+".R"+str(i)			
			pars['dateTrait'] = getDateTraitBlock(taxa, dates)
			makeXMLFile(pars, template, outputfile=pars["name"], outputpath=outputpath)

			# Write command to scripts
			for seed in seeds:
				beastcmd = '%s -seed %s -overwrite %s.xml' % (pars["beast"], pars["seed"], pars["name"])
				eulercmd = 'sbatch --time=%s:00:00 --ntasks=%s --job-name="%s" --output="%s.out" --error="%s.err" --wrap="%s"' % (pars["hours"], pars["ntasks"], pars["name"], pars["name"], pars["name"], beastcmd)
				scriptfile.write("%s\n" % (eulercmd))
					
		#
		scriptfile.close()
	#
#
