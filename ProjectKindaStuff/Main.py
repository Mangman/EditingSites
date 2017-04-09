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

parser.add_argument('-f4', '--fourthFilePath',        type=str, 
													 help='tsv file with mismatch information',
													 default = None)

parser.add_argument('-f5', '--fifthFilePath',        type=str, 
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
second_site_df = TsvReader.read_tsv(args.secondFilePath)
third_site_df = TsvReader.read_tsv(args.thirdFilePath)
fourth_site_df = TsvReader.read_tsv(args.fourthFilePath)
fifth_site_df = TsvReader.read_tsv(args.fifthFilePath)
#-------------------------------------

#-------------------------------------
#  Filtering sites

	# ft.filter_two (first_site_df, second_site_df, "r", args)
	# ft.filter_two (first_site_df,  third_site_df, "g", args)
	# ft.filter_two (second_site_df, third_site_df, "b", args)

	# ft.filter_three (first_site_df, second_site_df, third_site_df, "m", args)

	# red_patch = mpatches.Patch(color='red', label="TNOR1-TNOR2")
	# green_patch = mpatches.Patch(color='green', label="TNOR1-TNOR3")
	# blue_patch = mpatches.Patch(color='blue', label="TNOR2-TNOR3")
	# magenta_patch = mpatches.Patch(color='magenta', label="TNOR1-TNOR2-TNOR3")

	# plt.legend(handles=[red_patch, green_patch, blue_patch, magenta_patch])

	# plt.xlabel('filtering stages') 
	# plt.ylabel('number of clusters')

	# plt.show()

	#----------------------------------------
	# ft.filter_two (first_site_df, second_site_df, "r", args)

	# red_patch = mpatches.Patch(color='red', label="THYP2-THYP3")

	# plt.legend(handles=[red_patch])

	# plt.xlabel('filtering stages') 
	# plt.ylabel('number of clusters')

	# plt.show()
	#----------------------------------------
	# ft.filter_three (first_site_df, second_site_df, third_site_df, "m", args)


	# red_patch = mpatches.Patch(color='red', label="TNOR1")
	# green_patch = mpatches.Patch(color='green', label="TNOR2")
	# blue_patch = mpatches.Patch(color='blue', label="TNOR3")

	# plt.legend(handles=[red_patch, green_patch, blue_patch])

	# plt.xlabel('filtering stages') 
	# plt.ylabel('number of clusters')

	# plt.show()

	#-------------------------------------

	# ft.filter_two (first_site_df, second_site_df, "m", args)

	# red_patch = mpatches.Patch(color='red', label="THYP2")
	# green_patch = mpatches.Patch(color='green', label="THYP3")

	# plt.legend(handles=[red_patch, green_patch])

	# plt.xlabel('filtering stages') 
	# plt.ylabel('number of clusters')

	# plt.show()
#-------------------------------------

#-------------------------------------
#  Cluster Filtering
	#-------------------------------------

	# filtered_first = ft.filter_df(first_site_df, args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)
	# filtered_second= ft.filter_df(second_site_df, args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)
	# filtered_third = ft.filter_df(third_site_df, args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)

	#-------------------------------------

	# cluster_params = []
	# cluster_lengths = []

	# for cluster_length in tqdm(range (5,100,5)) :
	# 	cluster_params.append(cluster_length)
	# 	cluster_lengths.append(ft.get_cluster_quantitive_info  (filtered_first, 500, cluster_length)[0])

	# plt.plot (cluster_params, cluster_lengths, color='m')
	# plt.xlabel('min #\n(max site distance - 500)') 
	# plt.ylabel('number of clusters')

	# cluster_params = []
	# cluster_lengths = []

	# for cluster_length in tqdm(range (5,100,5)) :
	# 	cluster_params.append(cluster_length)
	# 	cluster_lengths.append(ft.get_cluster_quantitive_info  (filtered_first, 900, cluster_length)[0])

	# plt.plot (cluster_params, cluster_lengths, color='b')
	#------------------------

	#------------------------
	# cluster_params = []
	# cluster_lengths = []

	# for distance in tqdm(range (300,1500,100)) :
	# 	cluster_params.append(distance)
	# 	cluster_lengths.append(ft.get_cluster_quantitive_info  (filtered_first, distance, 5)[0])

	# plt.plot (cluster_params, cluster_lengths, color='b')
	# plt.xlabel('max site distance\n(min # - 5)') 
	# plt.ylabel('number of clusters')
	# #------------------------

	#------------------------
	# data = ft.get_cluster_quantitive_info  (filtered_first, args.clusterDistance, args.clusterLength)

	# plt.hist(data[2], 150)
	# plt.xlabel('length of cluster\n(min len - 10, max dist - 500)') 
	# plt.ylabel('number of clusters')
	#------------------------

	#------------------------
	# data = ft.get_cluster_quantitive_info  (filtered_first, args.clusterDistance, args.clusterLength)

	# plt.hist(data[1], 150)
	# plt.xlabel('number of sites\n(min len - 10, max dist - 500)') 
	# plt.ylabel('number of clusters')
	#------------------------

	#------------------------
	# data = ft.get_cluster_quantitive_info  (filtered_first, args.clusterDistance, args.clusterLength)

	# plt.plot(data[2], data[1], '.')
	# plt.xlabel('length of cluster\n(min len - 10, max dist - 500)') 
	# plt.ylabel('number of sites')
	#------------------------
	# plt.show()
	#-------------------------------------
#-------------------------------------

#-------------------------------------
#  Cluster Reproducibility


	# print "TNOR1-TNOR2 -- {0}".format(ft.find_intersecting_clusters_between_two (first_site_df, second_site_df, args.clusterDistance, args.clusterLength,
	#  							   								  args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)) 


	# print "TNOR2-TNOR3 -- {0}".format(ft.find_intersecting_clusters_between_two (second_site_df, third_site_df, args.clusterDistance, args.clusterLength,
	#  							   								  args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)) 

	# print "TNOR1-TNOR3 -- {0}".format(ft.find_intersecting_clusters_between_two (first_site_df, third_site_df, args.clusterDistance, args.clusterLength,
	#  							   								  args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)) 

	# ft.find_intersecting_clusters_between_three (first_site_df, second_site_df, third_site_df, args.clusterDistance, args.clusterLength,
	# 							                 args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution) 

	# ft.check_reproducibility_between_three (first_site_df, second_site_df, third_site_df, args.clusterDistance, args.clusterLength,
	#											 args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution) 
#-------------------------------------

#-------------------------------------
 # Finding common
	# common = Tools.find_common_between_data_frames(first_filtered, second_filtered)
	# print common
#-------------------------------------

#-------------------------------------
#  Common THYP TNOR 

# filtered_first_THYP  = ft.filter_df(first_site_df, args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)
# filtered_second_THYP = ft.filter_df(second_site_df, args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)

# filtered_first_TNOR  = ft.filter_df(third_site_df, args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)
# filtered_second_TNOR = ft.filter_df(fourth_site_df, args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)
# filtered_third_TNOR  = ft.filter_df(fifth_site_df, args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)

# #  Clusters
# THYP_intersection   = ft.find_intersecting_clusters_between_two(filtered_first_THYP, filtered_second_THYP, args.clusterDistance, args.clusterLength,
# 												 args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)

# TNOR12_intersection = ft.find_intersecting_clusters_between_two(filtered_first_TNOR, filtered_second_TNOR, args.clusterDistance, args.clusterLength,
# 												 args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)

# TNOR13_intersection = ft.find_intersecting_clusters_between_two(filtered_first_TNOR, filtered_third_TNOR, args.clusterDistance, args.clusterLength,
# 												 args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)

# TNOR23_intersection = ft.find_intersecting_clusters_between_two(filtered_second_TNOR, filtered_third_TNOR, args.clusterDistance, args.clusterLength,
# 	  											 args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)

# THYP_TNOR12 = ft.find_common(THYP_intersection[0], TNOR12_intersection[0])

# THYP_TNOR13 = ft.find_common(THYP_intersection[0], TNOR13_intersection[0])

# THYP_TNOR23 = ft.find_common(THYP_intersection[0], TNOR23_intersection[0])

# print "ab - {0}\nac - {1}\nbc - {2}".format(THYP_TNOR12[1], THYP_TNOR13[1], THYP_TNOR23[1])

# Rrr
# filtered_first_THYP  = ft.filter_df(first_site_df, args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)
# filtered_second_THYP = ft.filter_df(second_site_df, args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)

# filtered_first_TNOR  = ft.filter_df(third_site_df, args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)
# filtered_second_TNOR = ft.filter_df(fourth_site_df, args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)
# filtered_third_TNOR  = ft.filter_df(fifth_site_df, args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)

# #  Clusters
# THYP_intersection   = ft.find_intersecting_clusters_between_two(filtered_first_THYP, filtered_second_THYP, 
# 	args.clusterDistance, args.clusterLength, args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)[0]

# TNOR123_intersection = ft.check_reproducibility_between_three(filtered_first_TNOR, filtered_second_TNOR, filtered_third_TNOR, 
# 	args.clusterDistance, args.clusterLength, args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)

# THYP_TNOR123 = ft.find_common(THYP_intersection, TNOR123_intersection)[0]
# print(len(THYP_TNOR123[0])+len(THYP_TNOR123[1]))

#  THYP \ TNOR 

filtered_first_THYP  = ft.filter_df(first_site_df, args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)
filtered_second_THYP = ft.filter_df(second_site_df, args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)

filtered_first_TNOR  = ft.filter_df(third_site_df, args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)
filtered_second_TNOR = ft.filter_df(fourth_site_df, args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)
filtered_third_TNOR  = ft.filter_df(fifth_site_df, args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)

#  Clusters
THYP_intersection   = ft.find_intersecting_clusters_between_two(filtered_first_THYP, filtered_second_THYP, args.clusterDistance, args.clusterLength,
												 args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)

TNOR12_intersection = ft.find_intersecting_clusters_between_two(filtered_first_TNOR, filtered_second_TNOR, args.clusterDistance, args.clusterLength,
												 args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)

TNOR13_intersection = ft.find_intersecting_clusters_between_two(filtered_first_TNOR, filtered_third_TNOR, args.clusterDistance, args.clusterLength,
												 args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)

TNOR23_intersection = ft.find_intersecting_clusters_between_two(filtered_second_TNOR, filtered_third_TNOR, args.clusterDistance, args.clusterLength,
	  											 args.absoluteSubstitution, args.fullCoverage, args.relativeSubstitution)


THYP_TNOR12 = ft.find_common(THYP_intersection[0], TNOR12_intersection[0])[0]

THYP_TNOR13 = ft.find_common(THYP_intersection[0], TNOR13_intersection[0])[0]

THYP_TNOR23 = ft.find_common(THYP_intersection[0], TNOR23_intersection[0])[0]

rest = ft.exclude (THYP_intersection[0], THYP_TNOR12)
rest = ft.exclude (rest, THYP_TNOR13)
rest = ft.exclude (rest, THYP_TNOR23)

with open('pos_res.format', 'w') as pos_f, open("rev_res.format", 'w') as rev_f:
	old_out = sys.stdout
	

	# print(map(lambda x: x.pos, rest[0]))
	# print(map(lambda x: x.pos, rest[1]))


	pos_positions = [  (rest[0][i][0].seqnames, rest[0][i][0].pos, rest[0][i][-1].pos) for i in range(len(rest[0])) ]

	rev_positions = [  (rest[1][i][0].seqnames, rest[1][i][0].pos, rest[1][i][-1].pos) for i in range(len(rest[1])) ]

	sys.stdout = pos_f

	for cl in pos_positions :
		print cl[0], cl[1], cl[2]
	
	sys.stdout = rev_f

	for cl in rev_positions :
		print cl[0], cl[1], cl[2]


	sys.stdout = old_out

#print "ab - {0}\nac - {1}\nbc - {2}".format(THYP_TNOR12[1], THYP_TNOR13[1], THYP_TNOR23[1])
print "Num of THYP-only sites -- {0}".format(len(rest[0])+len(rest[1]))


#-------------------------------------