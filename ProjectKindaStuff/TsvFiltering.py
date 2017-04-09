# -*- coding: utf-8 -*-

import sys

import pandas as pd
pd.set_option('expand_frame_repr', False)

import matplotlib.pyplot as plt
import numpy as np
import time 
import collections
import operator

from Utilities import  progress_bar
import MismatchDataFrameTools as Tools
from Graph import Graph, Searcher


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


	return straight_strand_filtered[ straight_strand_filtered['G']/straight_strand_filtered["coverage" ].astype('float')  > threshold ]\
    .append(reverse_strand_filtered[  reverse_strand_filtered['C']/reverse_strand_filtered["coverage" ].astype('float') > threshold ])

#  Applies all the filters on sites dataFrame
def filter_df (site_df, absolute_substitution_threshold, coverage_threshold, relative_substitution_threshold) :
	# print "\n--------------------------------------------------------------\n\
	#        \rStarted DataFrame filtering.\n"

	tic = time.clock()

	abs_filtered = filter_by_absolute_TG_replacement(site_df, absolute_substitution_threshold)

	coverage_filtered = filter_by_coverage(abs_filtered, coverage_threshold)

	relative_filtered = filter_by_relative_substitution(coverage_filtered, relative_substitution_threshold)

	toc = time.clock()
	execution_duration = toc-tic

	# print "\nEndend filtering process. Time estimated - {0}\n\
	#    	   \r--------------------------------------------------------------\n"\
	#    	   .format(execution_duration)

	return relative_filtered

def filter_two (first_site_df, second_site_df, plot_color, args) :
	#  TODO: сделать абстрактно
	numbers = pd.DataFrame([ [0.]*4, [0.]*4, [0.]*4 ],
							columns= ["original", "first_phase", "second_phase", "third_phase"], 
							index =  ["first", "second", "merged"])

	original_common = Tools.find_common_between_data_frames(first_site_df, second_site_df)

	numbers["original"][0] = len(first_site_df)
	numbers["original"][1] = len(second_site_df)
	numbers["original"][2] = len(original_common)


	first_phase_result1  = filter_by_absolute_TG_replacement(first_site_df,  args.absoluteSubstitution)
	first_phase_result2  = filter_by_absolute_TG_replacement(second_site_df, args.absoluteSubstitution)
	first_phase_common = Tools.find_common_between_data_frames(first_phase_result1, first_phase_result2)

	numbers["first_phase"][0] = len(first_phase_result1)
	numbers["first_phase"][1] = len(first_phase_result2)
	numbers["first_phase"][2] = len(first_phase_common)

	second_phase_result1  = filter_by_coverage(first_phase_result1,  args.fullCoverage)
	second_phase_result2  = filter_by_coverage(first_phase_result2,  args.fullCoverage)
	second_phase_common = Tools.find_common_between_data_frames(second_phase_result1, second_phase_result2)

	numbers["second_phase"][0] = len(second_phase_result1)
	numbers["second_phase"][1] = len(second_phase_result2)
	numbers["second_phase"][2] = len(second_phase_common)

	third_phase_result1  = filter_by_relative_substitution(second_phase_result1,  args.relativeSubstitution)
	third_phase_result2  = filter_by_relative_substitution(second_phase_result2,  args.relativeSubstitution)
	third_phase_common = Tools.find_common_between_data_frames(third_phase_result1, third_phase_result2)

	numbers["third_phase"][0] = len(third_phase_result1)
	numbers["third_phase"][1] = len(third_phase_result2)
	numbers["third_phase"][2] = len(third_phase_common)

	print numbers

	#plt.plot (np.array([0, 1, 2, 3]), np.array(numbers.loc["merged"][:]), color=plot_color)
	plt.plot (np.array([0, 1, 2, 3]), np.array(numbers.loc["first"][:]), color='r')
	plt.plot (np.array([0, 1, 2, 3]), np.array(numbers.loc["second"][:]), color='g')
	
def filter_three (first_site_df, second_site_df, third_site_df, plot_color, args) :
	#  TODO: сделать абстрактно
	numbers = pd.DataFrame([ [0.]*4, [0.]*4, [0.]*4, [0.]*4],
							columns= ["original", "first_phase", "second_phase", "third_phase"], 
							index =  ["first", "second", "third", "merged"])

	original_common = Tools.find_common_between_data_frames(first_site_df, second_site_df)
	original_common = Tools.find_common_between_data_frames(original_common, third_site_df)


	numbers["original"][0] = len(first_site_df)
	numbers["original"][1] = len(second_site_df)
	numbers["original"][2] = len(third_site_df)
	numbers["original"][3] = len(original_common)


	first_phase_result1  = filter_by_absolute_TG_replacement(first_site_df,  args.absoluteSubstitution)
	first_phase_result2  = filter_by_absolute_TG_replacement(second_site_df, args.absoluteSubstitution)
	first_phase_result3  = filter_by_absolute_TG_replacement(third_site_df, args.absoluteSubstitution)
	
	first_phase_common = Tools.find_common_between_data_frames(first_phase_result1, first_phase_result2)
	first_phase_common = Tools.find_common_between_data_frames(first_phase_common, first_phase_result3)


	numbers["first_phase"][0] = len(first_phase_result1)
	numbers["first_phase"][1] = len(first_phase_result2)
	numbers["first_phase"][2] = len(first_phase_result3)
	numbers["first_phase"][3] = len(first_phase_common)

	second_phase_result1  = filter_by_coverage(first_phase_result1,  args.fullCoverage)
	second_phase_result2  = filter_by_coverage(first_phase_result2,  args.fullCoverage)
	second_phase_result3  = filter_by_coverage(first_phase_result3,  args.fullCoverage)
	
	second_phase_common = Tools.find_common_between_data_frames(second_phase_result1, second_phase_result2)
	second_phase_common = Tools.find_common_between_data_frames(second_phase_common, second_phase_result3)


	numbers["second_phase"][0] = len(second_phase_result1)
	numbers["second_phase"][1] = len(second_phase_result2)
	numbers["second_phase"][2] = len(second_phase_result3)
	numbers["second_phase"][3] = len(second_phase_common)

	third_phase_result1  = filter_by_relative_substitution(second_phase_result1,  args.relativeSubstitution)
	third_phase_result2  = filter_by_relative_substitution(second_phase_result2,  args.relativeSubstitution)
	third_phase_result3  = filter_by_relative_substitution(second_phase_result3,  args.relativeSubstitution)


	third_phase_common = Tools.find_common_between_data_frames(third_phase_result1, third_phase_result2)
	third_phase_common = Tools.find_common_between_data_frames(third_phase_common, third_phase_result3)


	numbers["third_phase"][0] = len(third_phase_result1)
	numbers["third_phase"][1] = len(third_phase_result2)
	numbers["third_phase"][2] = len(third_phase_result3)
	numbers["third_phase"][3] = len(third_phase_common)

	print numbers

	plt.plot (np.array([0, 1, 2, 3]), np.array(numbers.loc["merged"][:]), color=plot_color)
	# plt.plot (np.array([0, 1, 2, 3]), np.array(numbers.loc["first"][:]), color='r')
	# plt.plot (np.array([0, 1, 2, 3]), np.array(numbers.loc["second"][:]), color='g')
	# plt.plot (np.array([0, 1, 2, 3]), np.array(numbers.loc["third"][:]), color='b')

def mean_coverage (value, full1, full2, full3) :
	return (value/float(full1) + value/float(full2) + value/float(full3))/3


#  Finds sites situated close to each other based on threholds
def get_clusters (site_df, distance, c_len) :
	straight_strand_filtered = site_df[site_df["strand"] == '+']
	reverse_strand_filtered  = site_df[site_df["strand"] == '-']

	clusters = [fill_clusters(straight_strand_filtered, distance, c_len), 
				fill_clusters(reverse_strand_filtered, distance, c_len) ]

	return clusters

def fill_clusters (site_df, distance, c_len) :
	#print site_df

	clusters = []
	
	current_cluster = []
	last_pos = site_df["pos"][site_df.index[0]]
	last = site_df.ix[site_df.index[0]]

	current_cluster.append(last_pos)
	
	i = 0
	for current in site_df.itertuples() :
		if i == 0 :
			i = 1
			continue

		current_pos = current.pos
		current_dist = current_pos-last_pos
		
		if current_dist <= distance and current_dist > 0 :
			current_cluster.append(current)
		else :
			if len(current_cluster) >= c_len : 
				clusters.append(current_cluster)
			
			current_cluster = []
			current_cluster.append(current)

		last_pos = current_pos
		last = current

	return clusters

def quantitive_info (clusters) :
	num_of_clusters = len(clusters[0])+len(clusters[1])
	
	df_clusters = []
	qty_of_each = []

	len_of_each = []

	for strand in clusters :
		for cl in strand :
			qty_of_each.append(len(cl))
			len_of_each.append(cl[-1].pos-cl[0].pos)

	return (num_of_clusters, qty_of_each, len_of_each)

def get_cluster_quantitive_info (site_df, distance, c_len) :
	
	clusters = get_clusters(site_df, distance, c_len)

	return quantitive_info(clusters)

def find_intersecting_clusters_between_two (site_df1, site_df2, clusterDistance, clusterLength, 
								absolute_substitution_threshold, coverage_threshold, relative_substitution_threshold) :
	filtered1 = filter_df (site_df1, absolute_substitution_threshold, coverage_threshold, relative_substitution_threshold)
	filtered2 = filter_df (site_df2, absolute_substitution_threshold, coverage_threshold, relative_substitution_threshold)
	
	clusters1 = get_clusters(filtered1, clusterDistance, clusterLength)
	clusters2 = get_clusters(filtered2, clusterDistance, clusterLength)

	filtered_clusters, number_of_intersections, percents = find_common (clusters1, clusters2)
	
	# print len(clusters1[0]) + len(clusters1[1])
	# print len(clusters2[0]) + len(clusters2[1])
	
	print "Number of intersections -- {0}".format(number_of_intersections)

	site_reproduction = np.array(percents).mean()

	return filtered_clusters, number_of_intersections, site_reproduction

def find_common (clusters1, clusters2) :
	filtered_clusters = []
	number_of_intersections = 0
	percents = []

	#  TODO: слишком просто
	for strand1, strand2 in zip(clusters1, clusters2) :
		filtered_clusters.append([])
		cl1_pos = 0

		for cl1 in strand1 :
			
			cl2_pos = 0
			iter_strand = iter(strand2)
			
			for cl2 in iter_strand :
				if cl1[0].seqnames == cl2[0].seqnames :
					intersects, merged, percentage = check_cluster_intersection(cl1, cl2, 0)
					
					# val = [Bool, merged_cluster]
					if intersects :
						
						number_of_intersections += 1
						filtered_clusters[-1].append(merged)
						percents.append(percentage)

						number_of_intersections += find_consecutive(cl1, iter_strand, cl1_pos, cl2_pos, filtered_clusters, percents, 0)
						break
				cl2_pos += 1
			cl1_pos += 1
	
	return filtered_clusters, number_of_intersections, percents

def exclude (clusters1, clusters2) :
	filtered_clusters = []

	for strand1, strand2 in zip(clusters1, clusters2) :
		filtered_clusters.append([])
		for cl in strand1 :
			if find_cluster(cl, strand2) == False :
				filtered_clusters[-1].append(cl)
	return filtered_clusters

def find_cluster(cluster, clusters) :
	for cl in clusters :
		if cl[0].seqnames == cluster[0].seqnames :
			if cluster[0] in cl:
				return True 
	return False

def find_consecutive (cl, iter_strand, cl1_pos, cl2_pos, filtered_clusters, percents, intersections) :
	try :
		next_cl = next(iter_strand) 
		if cl[0].seqnames == next_cl[0].seqnames :
			intersects, merged, percentage = check_cluster_intersection(cl, next_cl, 0)
			if intersects :
				
				intersections += 1
				filtered_clusters[-1].append(merged)
				percents.append(percentage)

				return find_consecutive(cl, iter_strand, cl1_pos, cl2_pos+1, filtered_clusters, percents, intersections)
		return intersections
	except :
		return intersections							


def find_common_with_graph (clusters1, clusters2, graph, name1, name2) :
	#  TODO: слишком просто
	number_of_intersections = 0

	for strand1, strand2, graph in zip(clusters1, clusters2, graph) :
		cl1_pos = 0
	
		for cl1 in strand1 :
			
			cl2_pos = 0
			iter_strand2 = iter(strand2)

			for cl2 in iter_strand2 :
				if cl1[0].seqnames == cl2[0].seqnames :
					val = check_cluster_intersection(cl1, cl2, 0)
					# val = [Bool, merged_cluster]
					if val[0] == True :
						number_of_intersections += 1
						graph.add_edge (name1+"{}".format(cl1_pos), name2+"{}".format(cl2_pos))
						number_of_intersections += find_consecutive_g(cl1, iter_strand2, cl1_pos, cl2_pos, graph, name1, name2, 0)
						break
				cl2_pos += 1
			cl1_pos += 1
	return number_of_intersections

def find_consecutive_g (cl, iter_strand, cl1_pos, cl2_pos, graph, name1, name2, intersections) :
	try :
		next_cl = next(iter_strand) 
		if cl[0].seqnames == next_cl[0].seqnames :
			val = check_cluster_intersection(cl, next_cl, 0)
			if val[0] == True :
				intersections += 1
				graph.add_edge (name1+"{}".format(cl1_pos), name2+"{}".format(cl2_pos+1))
				find_consecutive(cl, iter_strand, cl1_pos, cl2_pos+1, graph, name1, name2)
		return intersections
	except :
		return intersections							

def find_intersecting_clusters_between_three (site_df1, site_df2, site_df3, clusterDistance, clusterLength, 
								absolute_substitution_threshold, coverage_threshold, relative_substitution_threshold) :
	graph = [Graph(), Graph()]

	filtered1 = filter_df (site_df1, absolute_substitution_threshold, coverage_threshold, relative_substitution_threshold)
	filtered2 = filter_df (site_df2, absolute_substitution_threshold, coverage_threshold, relative_substitution_threshold)
	filtered3 = filter_df (site_df3, absolute_substitution_threshold, coverage_threshold, relative_substitution_threshold)
	
	clusters1 = get_clusters(filtered1, clusterDistance, clusterLength)
	clusters2 = get_clusters(filtered2, clusterDistance, clusterLength)
	clusters3 = get_clusters(filtered3, clusterDistance, clusterLength)

	print ("First clusters -- {0}\nSecond clusters -- {1}\nThird clusters -- {2}").format(
															len(clusters1[0])+len(clusters1[1]),
															len(clusters2[0])+len(clusters2[1]),
															len(clusters3[0])+len(clusters3[1]))
	
	ab_len = find_common_with_graph(clusters1, clusters2, graph, "a", "b")
	bc_len = find_common_with_graph(clusters2, clusters3, graph, "b", "c")
	ac_len = find_common_with_graph(clusters1, clusters3, graph, "a", "c")

	print ("ab -- {0}\nbc -- {1}\nac -- {2}".format(ab_len, bc_len, ac_len))

	searcher = Searcher()

	searcher.calculateClusters(graph[0], "res1.txt")

	searcher.calculateClusters(graph[1], "res2.txt")

	#print graph[0].edges
	#print graph[1].edges

def check_cluster_intersection (cl1, cl2, intersection) :
	# left bound, right bound
	lb1 = cl1[0 ].pos
	rb1 = cl1[-1].pos

	lb2 = cl2[0 ].pos
	rb2 = cl2[-1].pos

	#Sort
	if (rb1 >= rb2 and lb1 <= rb2-intersection):
		new_cluster = sorted(cl1+cl2, key=operator.attrgetter('pos'))	

		non_common_sites = set(map(lambda x: x.pos, cl1)+map(lambda x: x.pos, cl2))
	

		return True, new_cluster, mean_percentage(len(cl1), len(cl2), len(non_common_sites))

	elif (lb1 <= lb2 and rb1 >= lb2+intersection):
		new_cluster = sorted(cl1+cl2, key=operator.attrgetter('pos'))		
		
		non_common_sites = set(map(lambda x: x.pos, cl1)+map(lambda x: x.pos, cl2))
	
		return True, new_cluster, mean_percentage(len(cl1), len(cl2), len(non_common_sites))
	
	elif (lb1 >= lb2 and rb1 <= rb2) :
		new_cluster = sorted(cl1+cl2, key=operator.attrgetter('pos'))		
		
		non_common_sites = set(map(lambda x: x.pos, cl1)+map(lambda x: x.pos, cl2))
	
		return True, new_cluster, mean_percentage(len(cl1), len(cl2), len(non_common_sites))
	else : 
		return False, [], -1

def mean_percentage(len1, len2, len_intersection) :
	diff = (len1+len2-len_intersection)
	percentage = 100*(diff/float(len1)+diff/float(len2))/2
	return percentage

def check_reproducibility_between_three (site_df1, site_df2, site_df3, clusterDistance, clusterLength, 
								absolute_substitution_threshold, coverage_threshold, relative_substitution_threshold) :
	
	filtered1 = filter_df (site_df1, absolute_substitution_threshold, coverage_threshold, relative_substitution_threshold)
	filtered2 = filter_df (site_df2, absolute_substitution_threshold, coverage_threshold, relative_substitution_threshold)
	filtered3 = filter_df (site_df3, absolute_substitution_threshold, coverage_threshold, relative_substitution_threshold)
	
	clusters1 = get_clusters(filtered1, clusterDistance, clusterLength)
	clusters2 = get_clusters(filtered2, clusterDistance, clusterLength)
	clusters3 = get_clusters(filtered3, clusterDistance, clusterLength)

	common_clusters = []
	with open("clusters.txt") as f:
		common_clusters = eval(f.readlines()[0])
	
	percents = []

	sites = []

	all_the_stuff = []

	for strand1, strand2, strand3, common_clusters in zip(clusters1, clusters2, clusters3, common_clusters) :
		sites.append([])
		all_the_stuff.append([])

		
		for cl in common_clusters :
			sites1 = []
			sites2 = []
			sites3 = []

			all_the_stuff1 = []
			all_the_stuff2 = []
			all_the_stuff3 = []


			for little_cl in cl :
				if little_cl[0] == 'a' :
					sites1 += get_sites(strand1, int(little_cl[1:]))
					all_the_stuff1 += (strand1[int(little_cl[1:])])
				if little_cl[0] == 'b' :
					sites2 += get_sites(strand2, int(little_cl[1:]))
					all_the_stuff2 +=(strand2[int(little_cl[1:])])
				if little_cl[0] == 'c' :
					sites3 += get_sites(strand3, int(little_cl[1:]))
					all_the_stuff3 += (strand3[int(little_cl[1:])])	
			
			common_sites = common_for_three(sites1,sites2,sites3)

			all_the_stuff[-1].append(sorted(all_the_stuff1+all_the_stuff2+all_the_stuff3, key=operator.attrgetter('pos')))
			sites[-1] += list(set(sites1+sites2+sites3))

			percents.append(mean_percentage3(len(sites1), len(sites2), len(sites3), len(common_sites)))
	
		
	mean = np.array(percents).mean()
	print "Mean common sites percentage -- {0}".format(mean)
	print "Number of intersections -- {0}".format(len(all_the_stuff[0])+len(all_the_stuff[1]))



	return all_the_stuff


def common_for_three (sites1, sites2, sites3) :
	return set(sites1) & set(sites2) & set(sites3)

def mean_percentage3(len1, len2, len3, len_intersection) :
	percentage = 100*(len_intersection/float(len1)+len_intersection/float(len2)+len_intersection/float(len3))/3
	return percentage


def get_sites(clusters, cluster_num) :
	return map(lambda x: x.pos, clusters[cluster_num])

