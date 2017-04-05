# -*- coding: utf-8 -*-

import sys

import pandas as pd
pd.set_option('expand_frame_repr', False)

import time 

from Utilities import  progress_bar

#  Filters read coverage
def filter_by_coverage (site_df, threshold) :
	return site_df[site_df["coverage"] > threshold]


#  Filters absolute edit value
def filter_by_absolute_TG_replacement(site_df, threshold) :
	straight_strand_filtered   = site_df[site_df["strand"] == '+']
	reverse_strand_filtered  = site_df[site_df["strand"] == '-']

	return straight_strand_filtered[ straight_strand_filtered['G'] > threshold ]\
    .append(reverse_strand_filtered[  reverse_strand_filtered['C'] > threshold ])

#  Filters by ratio of substituion coverage to whole coverage
def filter_by_relative_substitution (site_df, threshold) :
	straight_strand_filtered = site_df[site_df["strand"] == '+']
	reverse_strand_filtered  = site_df[site_df["strand"] == '-']


	return straight_strand_filtered[ straight_strand_filtered['G']/straight_strand_filtered["coverage" ]  > threshold ]\
    .append(reverse_strand_filtered[  reverse_strand_filtered['C']/ reverse_strand_filtered["coverage" ]  > threshold ])


#  Applies all the filters on sites dataFrame
def filter (site_df, coverage_threshold, absolute_substitution_threshold) :
	print "\n--------------------------------------------------------------\n\
	       \rStarted DataFrame filtering.\n"

	tic = time.clock()

	#  TODO: тут грязновато
	filtered = filter_by_absolute_TG_replacement\
		  	  (filter_by_coverage(site_df, coverage_threshold), absolute_substitution_threshold)

	toc = time.clock()
	execution_duration = toc-tic

	print "\nEndend filtering process. Time estimated - {0}\n\
	   	   \r--------------------------------------------------------------\n"\
	   	   .format(execution_duration)

	return filtered
