# -*- coding: utf-8 -*-

import sys

import pandas as pd
pd.set_option('expand_frame_repr', False)

import time 

import argparse

import TsvFiltering 
import TsvReader
import CommonSitesFinder

#-------------------------------------
#  Argument parser initialization
#-------------------------------------
parser = argparse.ArgumentParser(description='Count mismathches on gene alignment.')

parser.add_argument('-f1', '--firstFilePath',        type=str, 
		 										     help='tsv file with mismatch information')

parser.add_argument('-f2', '--secondFilePath',       type=str, 
													 help='tsv file with mismatch information', 
													 default = None)

parser.add_argument('-f3', '--thirdFilePath',        type=str, 
													 help='tsv file with mismatch information',
													 default = None)


parser.add_argument('-fc',  '--fullCoverage', 	     type=int, 
													 help='minimal coverage for site',                  
													 default = 0)

parser.add_argument('-as', '--absoluteSubstitution', type=int, 
													 help='minimal coverage for edited site',        
													 default = 0)

args = parser.parse_args()
#-------------------------------------

#-------------------------------------
#  Initialising Dataframes
#-------------------------------------
first_site_df  = TsvReader.read_tsv(args.firstFilePath)
second_site_df = TsvReader.read_tsv(args.secondFilePath)
#-------------------------------------

#-------------------------------------
#  Filtering sites
#-------------------------------------
first_filtered  = TsvFiltering.filter( first_site_df, args.fullCoverage, args.absoluteSubstitution)
second_filtered = TsvFiltering.filter(second_site_df, args.fullCoverage, args.absoluteSubstitution)
#-------------------------------------

#-------------------------------------
#  Finding common
#-------------------------------------
common = CommonSitesFinder.find_common_between_data_frames(first_filtered, second_filtered)
print common
#-------------------------------------