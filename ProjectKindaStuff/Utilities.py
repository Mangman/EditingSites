# -*- coding: utf-8 -*-

import sys
import time

def count_file_lines (file) :
	#tic = time.clock()
		
	file_len = sum(1 for line in file)
	file.seek(0,0)
	
	#toc = time.clock()
	#print toc-tic

	return file_len

def progress_bar (current_iteration, total) :
	whole_len = 25
	
	fill = '█'
	percentage = (current_iteration/float(total))*100

	fill_len = whole_len*current_iteration//total

	sys.stdout.write ("\rProgress [{0:s}{1:s}] {2:f}%".format(fill*fill_len, '-'*(whole_len-fill_len), percentage))
