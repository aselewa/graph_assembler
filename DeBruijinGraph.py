#!/usr/local/bin/python
import os
import sys
import json
import numpy as np
import argparse
import networkx as rx
import matplotlib.pyplot as plt
import random as rn
from networkx.readwrite import json_graph

parser = argparse.ArgumentParser()
parser.add_argument("-f", metavar="--filename", help="FASTA filename")
parser.add_argument("-k", metavar="--kmer length",help="kmer length", type=int)
args = parser.parse_args()

class DeBruijnGraph:

    """A De Bruijin multigraph from a collection of strings. User supplies 
    sequence string and kmer length k. Nodes correspond to the left and right
    k-1mers of each kmer."""

    def chop(self,st, k):
        """ Convert string into a set of kmers """
        for i in xrange(0,len(st)-k+1): yield st[i:i+k]

    class Node:
        """ Node in a De Bruijin graph that represents k-1-mer """

        def __init__(self,km1mer):
            self.km1mer = km1mer
            self.nin = 0
            self.nout = 0

        def isSemiBalanced(self):
            return abs(self.nin - self.nout) > 0 

        def isBalanced(self):
            return self.nin == self.nout

        def __hash__(self):
            return hash(self.km1mer)

        def __len__(self):
            return len(self.km1mer)

    def __init__(self,reads,k):
        """ Builds the De Bruijn Graph """

        self.edges = {}         #Represents edges between neighbors
        self.nodes = {}     #Maps k-1mers to their Node object
        self.reduc_graph = {}
        self.reduc_nodes = []
        self.k = k
        self.semi, self.nbal, self.nn = 0, 0, 0
        self.tails, self.heads = [], []
        self.numReads = len(reads)

        for j in range(0,self.numReads):
            kmers = list(self.chop(reads[j],k))
            for i in range(0,len(kmers)):
                kmrL, kmrR = ''.join(kmers[i][:-1]), ''.join(kmers[i][1:]) #Split each k-mer into its left (k-1)mer and right (k-1)mer
                nodeL, nodeR = None, None
                if(self.nodes.setdefault(kmrL,None)): #if a node with kmrL exists already, nodeL points to this node
                    nodeL = self.nodes[kmrL]
                else:
                    nodeL = self.nodes[kmrL] = self.Node(kmrL) #if a node with kmrL doesn't exist, make a new node
                if(self.nodes.setdefault(kmrR,None)): #Same as above but for kmrR
                    nodeR = self.nodes[kmrR]
                else:
                    nodeR = self.nodes[kmrR] = self.Node(kmrR)
                if(nodeL.km1mer in self.edges): #if the left node exists in the edge space
                    if(self.edges[nodeL.km1mer].count(nodeR.km1mer) == 0): #if left is NOT already connected to the right
                        nodeL.nout += 1
                        nodeR.nin += 1
                        self.edges.setdefault(nodeL.km1mer,[]).append(nodeR.km1mer) #append right to the left
                else: 
                    nodeL.nout += 1 
                    nodeR.nin += 1
                    self.edges.setdefault(nodeL.km1mer,[]).append(nodeR.km1mer)

        for node in self.nodes.itervalues():
            if(node.isBalanced()):
                self.nbal += 1
            elif (node.isSemiBalanced()):
                if(node.nout == 0):
                    self.tails.append(node.km1mer)
                    self.edges.setdefault(node.km1mer,[])
                if(node.nin == 0):
                    self.heads.append(node.km1mer)
                self.semi += 1
            else:
                self.nn += 1

    def nnodes(self):
        return len(self.nodes)

    def nedges(self):
        return len(self.edges)
    
    def hasEularianPath(self):
        return self.nn == 0 and self.semi == 2

    def hasEularianCycle(self):
        return self.nn == 0 and self.semi == 0

    def findEularianPath(self):
        self.edges.setdefault(self.tail,[]).append(self.head)
        tour = []
        def visit(n):
            nxt = None
            while (len(self.edges[n]) > 0):
                nxt = self.edges[n].pop()
                visit(nxt)
            if(nxt):
                tour.append(nxt)
        visit(self.head)
        tour = tour[::-1][:-1]

        sprStr = tour[0].km1mer
        for i in range(1,len(tour)):
            sprStr += (tour[i].km1mer)[-1]
        return sprStr
    

    def getEdgeWeight(self, data, sample):
        counter = 0
        for k,v in data.items():
            if(k == sample[0]):
                for i in range(0,len(v)):
                    if(v[i] == sample[1]):
                        counter += 1
        return counter

    def printTopToBottom(self,current, parent):
        stack = []
        while(current):
            stack.append(current)
            current = parent[current]
        return stack

    def findLinearContig(self,start): #Builds a contig until a node has indegree!=1 (merge) or outdegree !=1 (split)
        nodeStack = []
        nodeStack.append(start)
        parent = {start:None}
        contig_list = []
        while(nodeStack):
            curr = nodeStack.pop()
            if(self.nodes[curr].nin == 1 and self.nodes[curr].nout == 1 or self.heads.count(curr) > 0):
                    parent.update({self.edges[curr][0]:curr})
                    nodeStack.append(self.edges[curr][0])
            else:
                contig_list += self.printTopToBottom(curr,parent)
                contig_list = contig_list[::-1]
                nodeStack = []
        contig = contig_list[0]
        for i in range(1,len(contig_list)):
            contig += (contig_list[i])[-1]
        return contig

    def reduceGraph(self): #Builds a reduced graph where each node is a contig defined by findLinearContig()
        for i in range(0,len(self.heads)):
            head = self.heads[i]
            callStack = []
            callStack.append(head)
            while(callStack):
                start = callStack.pop()
                curr = self.findLinearContig(start)

                if(self.reduc_nodes.count(curr) == 0):
                    self.reduc_nodes.append(curr)

                subtail = curr[-(self.k-1):]
                if(self.edges[subtail]):
                    for i in range(0,len(self.edges[subtail])):
                        self.reduc_graph.setdefault(curr,[]).append(self.findLinearContig(self.edges[subtail][i]))
                    for i in range(0,len(self.edges[subtail])):
                        callStack.append(self.edges[subtail][i])


if(__name__ == '__main__'):

    """
    filename = args.f
    k = args.k

    def read_fasta(fp):
        seq = []
        for line in fp:
            line = line.rstrip()
            if(line):
                if (line.startswith(">") == False):
                    seq.append(line)
        return seq

    fp = open(filename)
    seq = read_fasta(fp)

    G = DeBruijnGraph(seq,k)
    G.reduceGraph()
    """
    
    #test De Bruijin graph by synthesizing kmers
    L = 10000
    read_len = 150
    samples = 500
    k = 100
    basis = ['A','T','C','G']

    seq = []
    temp = ''.join(np.random.choice(basis, L))

    f = open('seq.txt','w')
    for i in range(samples):
        t = temp[:5000]+''.join(np.random.choice(basis,1000))+temp[5000:]
        f.write('>read %d\n%s\n' %(i,t))
        seq.append(t)
    f.close()

    G = DeBruijnGraph(seq,k)
    G.reduceGraph()

    #visualizing graph using networkx and force-graph
    gg = rx.DiGraph()

    f = open('seq_out.txt','w')
    c = 0
    for i in range(0,len(G.reduc_nodes)):
        if(len(G.reduc_nodes[i]) > 1000):
            c += 1
            f.write('>Contig %d, len: %d\n%s\n' % (c,len(G.reduc_nodes[i]),G.reduc_nodes[i]))
            gg.add_node(G.reduc_nodes[i], size = len(G.reduc_nodes[i]))
    f.close()

    for k,v in G.reduc_graph.items():
        for i in range(0,len(v)):
            gg.add_edge(k,v[i], name=k+(v[i])[-1])
    
    #rx.draw(gg, with_labels=True)
    #plt.show()
    d = json_graph.node_link_data(gg)
    json.dump(d, open('Graphing/graph.json','w'))


