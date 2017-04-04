# -*- coding: utf-8 -*-

import sys
import time 
from itertools import islice
import argparse

from Utilities import count_file_lines, progress_bar

#Tsv mismatch information container
class Mismatch :
	seqname = None #  str
	pos = None #  int
	strand = None #  char
	reference = None #  char
	ACGT = [] #  [float] 
	coverage = None # float

	def __init__ (self, seqname, pos, strand, reference, ACGT, coverage) :
		self.seqname = seqname
		self.pos = pos
		self.strand = strand
		self.reference = reference
		self.ACGT = ACGT
		self.coverage = coverage

	def __str__(self):
		return ' '+self.seqname+'  '+str(self.pos)+'  '+self.strand+'  '+self.reference+'  '+str(self.ACGT)+'  '+str(self.coverage)

class  ADAR_Mismatch_Base :
	#  Linked_list
	#  СПРОСИТЬ ИРУ
	__base = []
	base_len = 0

	def __init__ (self, file_path) :
		print "--------------------------------------------------------------\n\
	    	 \rStarted possible ADAR editing sites base initialization.\n"

		tic = time.clock()

		with open(file_path) as file :
			file_len = count_file_lines(file)

			i = 0
			for line in file :
				#  Get site information 
				line = line.strip().split()

				#  Checking can_be_ADAR_editing parameter value
				if line[-1] == "TRUE" :
					#  CSV parameters order:
					#  seqnames, pos, strand, reference, A, C, G, T, coverage, can_be_ADAR_editing
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

print "\n"
#--------------------------------
#Argument parser initialization 
parser = argparse.ArgumentParser(description='Count mismathches on gene alignment.')
parser.add_argument('-f1', '--firstFilePath',  type=str, help='tsv file with mismatch information')
parser.add_argument('-f2', '--secondFilePath', type=str, help='another tsv file with mismatch information', default = None)

args = parser.parse_args()
#--------------------------------

#--------------------------------
#Filtering sites
all_sites_qty = count_file_lines(open(args.firstFilePath))

sites = ADAR_Mismatch_Base (args.firstFilePath)
#--------------------------------

print "--------------------------------------------------------------\n\
	 \rProcess ended.\nPotentional sites - {0}.\nAll sites - {1}.\nPercentage - {2:f}%.\n\
  	 \r--------------------------------------------------------------\n"\
	 .format(sites.base_len, all_sites_qty, 100*sites.base_len/float(all_sites_qty))









