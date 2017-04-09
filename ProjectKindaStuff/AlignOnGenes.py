import sys

def check_cluster_intersection (lb1, rb1, lb2, rb2) :
	#Sort

	# 154206766 154207544
	# 154220179 154271510 
	if (rb1 >= rb2 and lb1 <= rb2):
		return True
	elif (lb1 <= lb2 and rb1 >= lb2):
		return True
	elif (lb1 >= lb2 and rb1 <= rb2) :
		return True
	else : 
		return False

with open("Results/pos_res.format") as pos_res,	     \
	 open("Results/rev_res.format") as rev_res,  \
	 open("gencode.v26.basic.annotation.gtf") as genes, \
	 open("result.genes", 'w') as result :

	old_out = sys.stdout;
	sys.stdout = result
	i = 0
 	for line in genes:
 		if i < 5 :
 			i+=1
 			continue
	 	line = line.strip().split()

	 	if line[2] == 'gene':
	 		strand    = line[ 6]
	 		seqname   = line[ 0]
	 		gene_name = line[13]
	 		lb = int(line[3])
	 		rb = int(line[4])

	 		if strand == '+':
		 		for cluster in pos_res :
		 			
		 			cluster = cluster.strip().split()
		 			cl_seqname = cluster[0]
		 	
		 			if cl_seqname == seqname :
		 				cl_lb  = int(cluster[1])
		 				cl_rb  = int(cluster[2])

		 				if check_cluster_intersection(lb, rb, cl_lb, cl_rb) == True:
		 					#print "strand -- '+'"
		 					print gene_name[1:-2]

		 		pos_res.seek(0,0)
		 	elif strand == '-':
		 		for cluster in pos_res.readlines() :
		 			
		 			cluster = cluster.strip().split()
		 			cl_seqname = cluster[0]
		 	
		 			if cl_seqname == seqname :
		 				cl_lb  = int(cluster[1])
		 				cl_rb  = int(cluster[2])

		 				if check_cluster_intersection(lb, rb, cl_lb, cl_rb) == True:
							#print "strand -- '-'"	
							#print cluster	 					
		 					print gene_name[1:-2]


		 		rev_res.seek(0,0)