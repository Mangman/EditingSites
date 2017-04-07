# -*- coding: utf-8 -*-

import numpy as np
import sys

class Graph :
	
	def  __init__(self):
		self.edges = {}
		self.names = set()

	def add_edge (self, first, second) :
		# TODO: сделать проверку на повторное добавление
		if first not in self.edges :
			self.edges[first] = []
		self.edges[first].append(second)

		if second not in self.edges :
			self.edges[second] = []
		self.edges[second].append(first)

		self.names.add(first)
		self.names.add(second)

class Searcher :
	
	def __init__(self) :
		self.visitedVerticies = []
		self.names = []

    #  Uses depth-first search in order to find connected component of the graph
	def calculateClusters (self, graph, filePath) :
		#   TODO: НОРМАЛЬНО РАЗБИТЬ НА ФУНКЦИИ
		self.names = list(graph.names)
		i = len(self.names)-1

		with open(filePath, "w") as file :
			component = 0
			while i >= 0 :
				self.__search (graph, self.names[i], file)

				component += 1
				
				clusterLen = i - len(self.names) + 1
				file.write("\nPrevious cluster ({})^\n\n".format(clusterLen))
				
				i = len(self.names)-1

			self.__init__()
			return component

	#  Depth-first search
	def search (self, graph, startVertex, outputFile) :
		self.names = list(graph.names)
		self.__search(graph, startVertex, outputFile)

	# It's private, dude
	def __search (self, graph, startVertex, outputFile) :
		index = self.names.index(startVertex)
		self.names.pop(index)
		self.visitedVerticies.append(startVertex)
		outputFile.write("{}\t".format(startVertex))
		if startVertex not in graph.edges :
			return
		for vertex in graph.edges[startVertex] :
			if vertex not in self.visitedVerticies :
				self.__search(graph, vertex, outputFile)