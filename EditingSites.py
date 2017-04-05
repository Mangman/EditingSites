# -*- coding: utf-8 -*-

import sys

import pandas as pd
pd.set_option('expand_frame_repr', False)

import numpy as np

import time 
from timeit import timeit

from itertools import islice
import argparse
from tqdm import tqdm

from Utilities import count_file_lines, progress_bar

#  TSV mismatch information container.
class Mismatch :
	# seqname = None #  str
	# pos = None #  int
	# strand = None #  char
	# reference = None #  char
	# ACGT = [] #  [float] 
	# coverage = None # float

	def __init__ (self, seqname, pos, strand, reference, ACGT, coverage) :
		self.seqname = seqname
		self.pos = pos
		self.strand = strand
		self.reference = reference
		self.ACGT = ACGT
		self.coverage = coverage

	def __str__(self):
		return ' {}  {}  {}  {}  {}  {}'\
		.format(self.seqname, self.pos, self.strand, self.reference, self.ACGT, self.coverage)

class  AdarMismatchBase :
	#  Linked_list
	#  СПРОСИТЬ ИРУ
	__base = []
	base_len = 0

	def __init__(self, file_path) :
		print "\n--------------------------------------------------------------\n\
	    	   \rStarted possible ADAR editing sites base initialization.\n"

		tic = time.clock()
		# TODO tqdm
		with open(file_path) as file:
			file_len = count_file_lines(file)

			i = 0
			for line in file:
				#  Get site information 
				line = line.strip().split()

				#  Checking can_be_ADAR_editing parameter value
				if line[-1] == "TRUE":
					#  TSV parameters order:
					#  seqnames, pos, strand, reference, A, C, G, T, coverage, can_be_ADAR_editing
					# TODO
					self.__base.append(Mismatch(line[0], line[1], line[2], line[3], line[4:8], line[8]))
					
				if i%7000 == 0 : progress_bar(i, file_len)
				i+=1

			self.base_len = len(self.__base)


			toc = time.clock()
			execution_duration = toc-tic
			
			progress_bar (1,1) 
			print "\n"

		print "Initialization completed.\nTime estimated -- {0}.\n\
		   	 \r--------------------------------------------------------------\n"\
		   	   .format(execution_duration)

	def __getitem__ (self, key) :
		if abs(key) > self.base_len : raise IndexError
		
		return self.base[key]

	def __iter__ (self) :
		return iter(self.base)

class CommonSitesFinder :

	#  Finds common potentional Adar editing sites between two sources.
	@staticmethod
	def find_common (f1, f2) :

		print "\n--------------------------------------------------------------\n\
	    	   \rStarted common possible ADAR editing sites search.\n"

		tic = time.clock()

		common_sites = 0

		#  {chr_name<str>: positions<set>}
		sites = {}

		#  Fills dictionary with sites from first file and uses it for comparison to second file
		CommonSitesFinder._fill_dict(sites, f1)

		print "  Searching common sites..."

		second_file_len = count_file_lines (f2)
		i = 0
		for line in f2 :
			line = line.strip().split()
			if line[-1] == "TRUE" :
				if line[1] in sites[line[0]] : 
					common_sites += 1 

			if i%7000 == 0: progress_bar (i, second_file_len)
			i += 1

		progress_bar (1,1)
		print "\n\n"


		toc = time.clock()

		execution_duration = toc-tic

		print "Search completed.\nTime estimated -- {0}.\n\
		   	 \r--------------------------------------------------------------\n"\
		   	   .format(execution_duration)

		return common_sites


	@staticmethod
	def _fill_dict (sites, f) :
		print "\n  Filling dictionary..."

		file_len = count_file_lines(f)

		i = 0
		for line in f :
			#  Getting site information.
			line = line.strip().split()
			
			#  Checking can_be_ADAR_editing parameter value.
			if line[-1] == "TRUE" :
				#  Узнать про стрэнд
				if line[0] not in sites.keys() :
					sites[line[0]] = set()
				sites[line[0]].add(line[1])		

			if i%7000 == 0 : progress_bar(i, file_len)
			i += 1

		progress_bar(1,1)	
		print "\n\n"		

class DecentTsvFiltering:

	#  Reads TSV file into dataFrame
	@staticmethod
	def read_tsv (file_path) :
		print "\n--------------------------------------------------------------\n\
		       \rStarted TSV file reading.\n"

		tic = time.clock()

		positions = []
		parameters = ["seqname", "pos", "strand", "reference", "A", "C", "G", "T", "coverage"]

		with open(file_path) as file:
			file_len = count_file_lines(file)

			i = 0

			progress_bar (0,1)

			for line in file:
				#  Get site information 
				line = np.array(line.strip().split())

				#  Checking can_be_ADAR_editing parameter value
				if line[-1] == "TRUE":
					# TODO: очень медленно 
					seqname = line[0]
					pos = int(line[1])
					strand = line[2]
					reference = line[3]
					the_rest = [float(line[x]) for x in range(4,9)]
					positions.append([seqname, pos, strand, reference]+the_rest)

				if i%7000 == 0 : progress_bar(i, file_len)
				i+=1

			toc = time.clock()
			execution_duration = toc-tic
			
			progress_bar (1,1) 
		print "\n"

		site_df = pd.DataFrame(positions, columns = parameters)

		print site_df[:5]
		print ("...")

		print "\nEndend reading process. Time estimated - {0}\n\
		   	   \r--------------------------------------------------------------\n"\
		   	   .format(execution_duration)
		
		return site_df

	#  Filters read coverage
	@staticmethod
	def filter_by_coverage (site_df, threshold) :
		return site_df[site_df["coverage"] > threshold]


	#  Filters absolute edit value
	@staticmethod 
	def filter_by_absolute_TG_replacement(site_df, threshold) :
		return site_df[site_df["G"] > threshold]\
	   .append(site_df[site_df["C"] > threshold])

	#  Applies all the filters on sites dataFrame
	@staticmethod 
	def filter (site_df, coverage_threshold, absolute_substitution_threshold) :
		print "\n--------------------------------------------------------------\n\
		       \rStarted DataFrame filtering.\n"

		tic = time.clock()

		#  TODO: тут грязновато
		filtered = DecentTsvFiltering.filter_by_absolute_TG_replacement\
			  	  (DecentTsvFiltering.filter_by_coverage(site_df, coverage_threshold), absolute_substitution_threshold)

		toc = time.clock()
		execution_duration = toc-tic

		return filtered

		print "\nEndend filtering process. Time estimated - {0}\n\
		   	   \r--------------------------------------------------------------\n"\
		   	   .format(execution_duration)

#---------------------------------
#  Argument parser initialization
#---------------------------------
parser = argparse.ArgumentParser(description='Count mismathches on gene alignment.')

parser.add_argument('-f1', '--firstFilePath',        type=str, help='tsv file with mismatch information')
parser.add_argument('-f2', '--secondFilePath',       type=str, help='another tsv file with mismatch information', default = None)

parser.add_argument('-c',  '--coverage', 	         type=int, help='minimal coverage for site',                  default = 0)
parser.add_argument('-as', '--absoluteSubstitution', type=int, help='minimal coverage for edited site',           default = 0)

args = parser.parse_args()
#---------------------------------

#---------------------------------
#  Filtering sites
#---------------------------------
#all_sites_qty = count_file_lines(open(args.firstFilePath))
#sites = AdarMismatchBase (args.firstFilePath)
#print "--------------------------------------------------------------\n\
#	 \rProcess ended.\nPotentional sites - {0}.\nAll sites - {1}.\nPercentage - {2:f}%.\n\
#  	 \r--------------------------------------------------------------\n"\
#	 .format(sites.base_len, all_sites_qty, 100*sites.base_len/float(all_sites_qty))
#---------------------------------

#---------------------------------
#  Filtering sites in  a right way
#---------------------------------
site_df = DecentTsvFiltering.read_tsv(args.firstFilePath)
filtered = DecentTsvFiltering.filter(site_df)

print filtered
#---------------------------------

#---------------------------------
#  Finding common sites
#---------------------------------
#
#with open(args.firstFilePath) as f1:
#	with open(args.secondFilePath) as f2:
#		print CommonSitesFinder.find_common(f1,f2)
#---------------------------------







