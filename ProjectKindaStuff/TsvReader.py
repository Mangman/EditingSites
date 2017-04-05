# -*- coding: utf-8 -*-

import sys

import pandas as pd
pd.set_option('expand_frame_repr', False)

import time 

from Utilities import progress_bar, count_file_lines

#  Reads TSV file into dataFrame
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
			line = line.strip().split()

			# TODO: очень медленно 
			seqname = line[0]
			pos = int(line[1])
			strand = line[2]
			reference = line[3]
			the_rest = [float(line[x]) for x in range(4,9)]
			positions.append([seqname, pos, strand, reference]+the_rest)

			if i%2000 == 0 : progress_bar(i, file_len)
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
