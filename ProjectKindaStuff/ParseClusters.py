from itertools import islice
import sys

with open ("res1.txt") as f1, open ("res2.txt") as f2, open ("clusters.txt", 'w') as result:
	clusters = [[],[]]
	for line in islice(f1, None, None, 3):
		current_cluster = []

		line = line.strip().split('\t')

		sets = set(line[i][0] for i in range(len(line)))
		if len(sets) == 3:
			clusters[0].append(line)

	for line in islice(f2, None, None, 3):
	 	current_cluster = []
		
		line = line.strip().split('\t')
		
		sets = set(line[i][0] for i in range(len(line)))
		if len(sets) == 3:
			clusters[1].append(line)

	print "Number of clusters -- {0}".format(len(clusters[0])+len(clusters[1]))

	sys.stdout = result

	print clusters

