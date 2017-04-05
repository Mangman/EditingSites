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

	site_df = pd.read_table(file_path) 
	site_df = site_df.drop('can_be_APOBEC_editing', axis=1)

	print site_df[:5]
	print ("...")

	print "\nEndend reading process.\n\
	   	   \r--------------------------------------------------------------\n"\

	
	return site_df
