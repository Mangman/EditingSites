# -*- coding: utf-8 -*-

import sys

import pandas as pd
pd.set_option('expand_frame_repr', False)

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from tqdm import tqdm

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

parser.add_argument('-d', '--clusterDistance',       type=int, 
													 help='maximal distance between cluster members',        
													 default = 0)

parser.add_argument('-cl', '--clusterLength',        type=int, 
													 help='minimal length of sites cluster',        
													 default = 0)


args = parser.parse_args()
#-------------------------------------

#-------------------------------------
#  Initialising Dataframes
#-------------------------------------
first_site_df  = TsvReader.read_tsv(args.firstFilePath)
# second_site_df = TsvReader.read_tsv(args.secondFilePath)

# if (args.thirdFilePath == None) :
# 	print ("oh no!")
# 	sys.exit()

third_site_df = TsvReader.read_tsv(args.thirdFilePath)

#-------------------------------------

#-------------------------------------
#  Filtering sites
#-------------------------------------
# first_filtered  = TsvFiltering.filter( first_site_df, args.fullCoverage, args.absoluteSubstitution)
# second_filtered = TsvFiltering.filter(second_site_df, args.fullCoverage, args.absoluteSubstitution)
# first_filtered  = TsvFiltering.filter_by_relative_substitution( first_site_df, args.relativeSubstitution)
# second_filtered = TsvFiltering.filter_by_relative_substitution(second_site_df, args.relativeSubstitution)



# ft.filter_two (first_site_df, second_site_df, "r", args)
# ft.filter_two (first_site_df,  third_site_df, "g", args)
# ft.filter_two (second_site_df, third_site_df, "b", args)

# ft.filter_three (first_site_df, second_site_df, third_site_df, "m", args)

# red_patch = mpatches.Patch(color='red', label="TNOR1-TNOR2")
# green_patch = mpatches.Patch(color='green', label="TNOR1-TNOR3")
# blue_patch = mpatches.Patch(color='blue', label="TNOR2-TNOR3")
# magenta_patch = mpatches.Patch(color='magenta', label="TNOR1-TNOR2-TNOR3")

# plt.legend(handles=[red_patch, green_patch, blue_patch, magenta_patch])

# plt.show()

#-------------------------------------


#-------------------------------------
#  Finding clusters
#-------------------------------------


#------------------------
# cluster_params = []
# cluster_lengths = []

# for cluster_length in tqdm(range (5,100,5)) :
# 	cluster_params.append(cluster_length)
# 	cluster_lengths.append(ft.find_clusters (first_site_df, 500, cluster_length)[0])

# plt.plot (cluster_params, cluster_lengths, color='m')
# plt.xlabel('min #\n(max site distance - 500)') 
# plt.ylabel('number of clusters')


# cluster_params = []
# cluster_lengths = []

# for cluster_length in tqdm(range (5,100,5)) :
# 	cluster_params.append(cluster_length)
# 	cluster_lengths.append(ft.find_clusters (first_site_df, 900, cluster_length)[0])

#plt.plot (cluster_params, cluster_lengths, color='b')
#------------------------


#------------------------
# cluster_params = []
# cluster_lengths = []

# for distance in tqdm(range (300,1500,100)) :
# 	cluster_params.append(distance)
# 	cluster_lengths.append(ft.find_clusters (first_site_df, distance, 5)[0])

# plt.plot (cluster_params, cluster_lengths, color='b')
# plt.xlabel('max site distance\n(min # - 5)') 
# plt.ylabel('number of clusters')
#------------------------



#------------------------
# data = ft.find_clusters (first_site_df, args.clusterDistance, args.clusterLength)

# plt.hist(data[2], 150)
# plt.xlabel('length of cluster\n(min len - 10, max dist - 500)') 
# plt.ylabel('number of clusters')
#------------------------

#------------------------
# data = ft.find_clusters (first_site_df, args.clusterDistance, args.clusterLength)

# plt.hist(data[1], 150)
# plt.xlabel('number of sites\n(min len - 10, max dist - 500)') 
# plt.ylabel('number of clusters')
#------------------------

#------------------------
data = ft.find_clusters (first_site_df, args.clusterDistance, args.clusterLength)

plt.plot(data[2], data[1], '.')
plt.xlabel('length of cluster\n(min len - 10, max dist - 500)') 
plt.ylabel('number of sites')
#------------------------



plt.show()


# with open("result.lol", 'w') as f :
# 	# for length in data[2] :
# 	# 	f.write(length)
# 	sys.stdout = f 
# 	for length in data[2] :
# 		print length


#-------------------------------------

#-------------------------------------
#  Finding common
#-------------------------------------
# common = Tools.find_common_between_data_frames(first_filtered, second_filtered)
# print common
#-------------------------------------