# -*- coding: utf-8 -*-

import sys

import pandas as pd
pd.set_option('expand_frame_repr', False)

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import numpy as np

import time 

import argparse

import TsvFiltering  as ft
import TsvReader
import MismatchDataFrameTools as Tools

#-------------------------------------
#  Argument parser initialization
#-------------------------------------
parser = argparse.ArgumentParser(description='Count mismathches on gene alignment.')

parser.add_argument('-f1', '--firstFilePath',        type=str, 
		 										     help='tsv file with mismatch information')

parser.add_argument('-f2', '--secondFilePath',       type=str, 
													 help='tsv file with mismatch information')

parser.add_argument('-f3', '--thirdFilePath',        type=str, 
													 help='tsv file with mismatch information',
													 default = None)


parser.add_argument('-fc',  '--fullCoverage', 	     type=int, 
													 help='minimal coverage for site',                  
													 default = 0)

parser.add_argument('-as', '--absoluteSubstitution', type=int, 
													 help='minimal coverage for edited site',        
													 default = 0)

parser.add_argument('-rs', '--relativeSubstitution', type=float, 
													 help='minimal coverage for relative\
													       substitution/coverage ratio',        
													 default = 0.05)


args = parser.parse_args()
#-------------------------------------

#-------------------------------------
#  Initialising Dataframes
#-------------------------------------
first_site_df  = TsvReader.read_tsv(args.firstFilePath)
second_site_df = TsvReader.read_tsv(args.secondFilePath)

if (args.thirdFilePath == None) :
	print ("oh no!")
	sys.exit()

third_site_df = TsvReader.read_tsv(args.thirdFilePath)

#-------------------------------------

#-------------------------------------
#  Filtering sites
#-------------------------------------
# first_filtered  = TsvFiltering.filter( first_site_df, args.fullCoverage, args.absoluteSubstitution)
# second_filtered = TsvFiltering.filter(second_site_df, args.fullCoverage, args.absoluteSubstitution)
# first_filtered  = TsvFiltering.filter_by_relative_substitution( first_site_df, args.relativeSubstitution)
# second_filtered = TsvFiltering.filter_by_relative_substitution(second_site_df, args.relativeSubstitution)

def filter_two (first_site_df, second_site_df, plot_color) :
	#  TODO: сделать абстрактно
	phases = [1,2,3]
	sites = []

	numbers = pd.DataFrame([ [0.]*4, [0.]*4, [0.]*4 ],
							columns= ["original", "first_phase", "second_phase", "third_phase"], 
							index =  ["first", "second", "merged"])

	original_common = Tools.find_common_between_data_frames(first_site_df, second_site_df)

	numbers["original"][0] = (len(first_site_df))
	numbers["original"][1] = (len(second_site_df))
	numbers["original"][2] = (len(original_common))


	first_phase_result1  = ft.filter_by_absolute_TG_replacement(first_site_df,  args.absoluteSubstitution)
	first_phase_result2  = ft.filter_by_absolute_TG_replacement(second_site_df, args.absoluteSubstitution)
	first_phase_common = Tools.find_common_between_data_frames(first_phase_result1, first_phase_result2)

	numbers["first_phase"][0] = (100*len(first_phase_result1)/float(numbers["original"][0]))
	numbers["first_phase"][1] = (100*len(first_phase_result2)/float(numbers["original"][1]))
	numbers["first_phase"][2] = (100*len(first_phase_common) /float(numbers["original"][2]))

	second_phase_result1  = ft.filter_by_coverage(first_phase_result1,  args.fullCoverage)
	second_phase_result2  = ft.filter_by_coverage(first_phase_result2,  args.fullCoverage)
	second_phase_common = Tools.find_common_between_data_frames(second_phase_result1, second_phase_result2)

	numbers["second_phase"][0] = (100*len(second_phase_result1)/float(numbers["original"][0]))
	numbers["second_phase"][1] = (100*len(second_phase_result2)/float(numbers["original"][1]))
	numbers["second_phase"][2] = (100*len(second_phase_common)/ float(numbers["original"][2]))

	third_phase_result1  = ft.filter_by_relative_substitution(second_phase_result1,  args.relativeSubstitution)
	third_phase_result2  = ft.filter_by_relative_substitution(second_phase_result2,  args.relativeSubstitution)
	third_phase_common = Tools.find_common_between_data_frames(third_phase_result1, third_phase_result2)

	numbers["third_phase"][0] = (100*len(third_phase_result1)/float(numbers["original"][0]))
	numbers["third_phase"][1] = (100*len(third_phase_result2)/float(numbers["original"][1]))
	numbers["third_phase"][2] = (100*len(third_phase_common) /float(numbers["original"][2]))

	numbers["original"][0] = 100
	numbers["original"][1] = 100
	numbers["original"][2] = 100

	print numbers

	plt.plot (np.array([0, 1, 2, 3]), np.array(numbers.loc["merged"][:]), color=plot_color)
	
def filter_three (first_site_df, second_site_df, third_site_df, plot_color) :
	#  TODO: сделать абстрактно
	phases = [1,2,3]
	sites = []

	numbers = pd.DataFrame([ [0.]*4, [0.]*4, [0.]*4, [0.]*4],
							columns= ["original", "first_phase", "second_phase", "third_phase"], 
							index =  ["first", "second", "third", "merged"])

	original_common = Tools.find_common_between_data_frames(first_site_df, second_site_df)
	original_common = Tools.find_common_between_data_frames(original_common, third_site_df)


	numbers["original"][0] = (len(first_site_df))
	numbers["original"][1] = (len(second_site_df))
	numbers["original"][2] = (len(third_site_df))
	numbers["original"][3] = (len(original_common))


	first_phase_result1  = ft.filter_by_absolute_TG_replacement(first_site_df,  args.absoluteSubstitution)
	first_phase_result2  = ft.filter_by_absolute_TG_replacement(second_site_df, args.absoluteSubstitution)
	first_phase_result3  = ft.filter_by_absolute_TG_replacement(third_site_df, args.absoluteSubstitution)
	
	first_phase_common = Tools.find_common_between_data_frames(first_phase_result1, first_phase_result2)
	first_phase_common = Tools.find_common_between_data_frames(first_phase_common, first_phase_result3)


	numbers["first_phase"][0] = (100*len(first_phase_result1)/float(numbers["original"][0]))
	numbers["first_phase"][1] = (100*len(first_phase_result2)/float(numbers["original"][1]))
	numbers["first_phase"][2] = (100*len(first_phase_result3)/float(numbers["original"][2]))
	numbers["first_phase"][3] = (100*len(first_phase_common) /float(numbers["original"][3]))

	second_phase_result1  = ft.filter_by_coverage(first_phase_result1,  args.fullCoverage)
	second_phase_result2  = ft.filter_by_coverage(first_phase_result2,  args.fullCoverage)
	second_phase_result3  = ft.filter_by_coverage(first_phase_result3,  args.fullCoverage)
	
	second_phase_common = Tools.find_common_between_data_frames(second_phase_result1, second_phase_result2)
	second_phase_common = Tools.find_common_between_data_frames(second_phase_common, second_phase_result3)


	numbers["second_phase"][0] = (100*len(second_phase_result1)/float(numbers["original"][0]))
	numbers["second_phase"][1] = (100*len(second_phase_result2)/float(numbers["original"][1]))
	numbers["second_phase"][2] = (100*len(second_phase_result2)/float(numbers["original"][2]))
	numbers["second_phase"][3] = (100*len(second_phase_common)/ float(numbers["original"][3]))

	third_phase_result1  = ft.filter_by_relative_substitution(second_phase_result1,  args.relativeSubstitution)
	third_phase_result2  = ft.filter_by_relative_substitution(second_phase_result2,  args.relativeSubstitution)
	third_phase_result3  = ft.filter_by_relative_substitution(second_phase_result3,  args.relativeSubstitution)


	third_phase_common = Tools.find_common_between_data_frames(third_phase_result1, third_phase_result2)
	third_phase_common = Tools.find_common_between_data_frames(third_phase_common, third_phase_result3)


	numbers["third_phase"][0] = (100*len(third_phase_result1)/float(numbers["original"][0]))
	numbers["third_phase"][1] = (100*len(third_phase_result2)/float(numbers["original"][1]))
	numbers["third_phase"][2] = (100*len(third_phase_result3)/float(numbers["original"][2]))
	numbers["third_phase"][3] = (100*len(third_phase_common) /float(numbers["original"][3]))

	numbers["original"][0] = 100
	numbers["original"][1] = 100
	numbers["original"][2] = 100
	numbers["original"][3] = 100

	print numbers

	plt.plot (np.array([0, 1, 2, 3]), np.array(numbers.loc["merged"][:]), color=plot_color)
	


filter_two (first_site_df, second_site_df, "r")
filter_two (first_site_df,  third_site_df, "g")
filter_two (second_site_df, third_site_df, "b")

filter_three (first_site_df, second_site_df, third_site_df, "m")

red_patch = mpatches.Patch(color='red', label="TNOR1-TNOR2")
green_patch = mpatches.Patch(color='green', label="TNOR1-TNOR3")
blue_patch = mpatches.Patch(color='blue', label="TNOR2-TNOR3")
magenta_patch = mpatches.Patch(color='magenta', label="TNOR1-TNOR2-TNOR3")

plt.legend(handles=[red_patch, green_patch, blue_patch, magenta_patch])


plt.show()

#-------------------------------------

#-------------------------------------
#  Finding common
#-------------------------------------
# common = Tools.find_common_between_data_frames(first_filtered, second_filtered)
# print common
#-------------------------------------