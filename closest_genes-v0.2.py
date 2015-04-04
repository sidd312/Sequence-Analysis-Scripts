#!/usr/bin/python
import numpy as np
import sys
import operator
import subprocess
import os

def findDistance(gene_start,gene_end,peak_start,peak_end,gene_location):
	#case 1: Upstream
	if(gene_location=='upstream'):
		return abs(peak_start-gene_end)
	
	#case 2: Downstream
	
	elif(gene_location=='downstream'):
		return abs(gene_start-peak_end)
	
	#case 3: Peak overlapping

	elif(gene_location=='peak_overlapping'):
		return 0 
		#according to closest features bedops
 		#return min(abs(gene_start-peak_start),abs(gene_end-peak_end))
	
	#case 4: Gene overlapping
	
	else:
		#assert gene_location =='gene_overlapping'
		return 0 
		#according to the closest features bedops
		#return min(abs(gene_start-peak_start),abs(gene_end-peak_end))
		#return abs(gene_end-peak_start)

def filter(peaks,genes,offset):

	#d=open("filtered"+str(offset)+".bed","wb")
	e=open("filtered"+str(offset)+"_explained"+".txt","wb+")
	print "Processing ...."
	for i in xrange(peaks.size):
		#i is the row and 0 is the column for the peaks file
		for j in xrange(genes.size):
		    # j is the row and 0 is the column for the genes file	
			if(peaks[i][0]==genes[j][0]):
				#Scaffold is matching and so now four cases
				#case 1: if there is a gene in the left side or upstream of the peak and not overlapping				
				#if(b[i][1]-offset>0) and (b[i][1]-offset)<=c[j][1]<=b[i][1] and (b[i][1]-offset)<=c[j][2]<=b[i][1]:
				if((peaks[i][1]-offset)<=genes[j][2] and genes[j][2]<=peaks[i][1]):
					e.write(peaks[i][0]+"\t"+str(peaks[i][1])+"\t"+str(peaks[i][2])+"\t"+str(peaks[i][3])+"\t"+str(peaks[i][4])+"\t"+str(genes[j][3])+"\t"+"upstream"+\
						"\t"+str(findDistance(genes[j][1],genes[j][2],peaks[i][1],peaks[i][2],"upstream"))+"\n")

				#case 2: if there is a gene in the right side or downstream of the peak

				elif(peaks[i][2]<=genes[j][1] and genes[j][1]<=(peaks[i][2]+offset)):
					e.write(peaks[i][0]+"\t"+str(peaks[i][1])+"\t"+str(peaks[i][2])+"\t"+str(peaks[i][3])+"\t"+str(peaks[i][4])+"\t"+str(genes[j][3])+"\t"+"downstream"+\
						"\t"+str(findDistance(genes[j][1],genes[j][2],peaks[i][1],peaks[i][2],"downstream"))+"\n")

				#case 3: if there is a gene overlapping with the peak
				#either the startpoint of the gene or the endpoint of the gene is in between the range of the peak 

				elif(peaks[i][1]<=genes[j][1] and genes[j][1]<=peaks[i][2]) or (peaks[i][1]<=genes[j][2] and genes[j][2]<=peaks[i][2]):
					e.write(peaks[i][0]+"\t"+str(peaks[i][1])+"\t"+str(peaks[i][2])+"\t"+str(peaks[i][3])+"\t"+str(peaks[i][4])+"\t"+str(genes[j][3])+"\t"+"gene_overlapping"+\
						"\t"+str(findDistance(genes[j][1],genes[j][2],peaks[i][1],peaks[i][2],"gene_overlapping"))+"\n")

				#case 4: peak lies in the gene
				elif(((peaks[i][1]-offset)<=genes[j][1]<=peaks[i][1]) and (peaks[i][2]<=genes[j][2]<=peaks[i][2]+offset)):
					e.write(peaks[i][0]+"\t"+str(peaks[i][1])+"\t"+str(peaks[i][2])+"\t"+str(peaks[i][3])+"\t"+str(peaks[i][4])+"\t"+str(genes[j][3])+"\t"+"peak_overlapping"+\
						"\t"+str(findDistance(genes[j][1],genes[j][2],peaks[i][1],peaks[i][2],"peak_overlapping"))+"\n")
				#case 5: peak lies in the gene1
				elif(genes[j][1]<=(peaks[i][1]-offset)<=genes[j][2] or genes[j][1]<=(peaks[i][2]+offset)<=genes[j][2]):
					e.write(peaks[i][0]+"\t"+str(peaks[i][1])+"\t"+str(peaks[i][2])+"\t"+str(peaks[i][3])+"\t"+str(peaks[i][4])+"\t"+str(genes[j][3])+"\t"+"peak_overlapping"+\
						"\t"+str(findDistance(genes[j][1],genes[j][2],peaks[i][1],peaks[i][2],"peak_overlapping"))+"\n")
	e.close()
	return 



def mycmp(item1,item2):
	if item1[1]>item2[1]:
		return 1
	elif item2[1]>item1[1]:
		return -1
	else:
		assert (item2[1]==item1[1])
		return 0
	return 0

"""
Not in use now!!!

def opposite(gene_location):
	if(gene_location=='upstream'):
		return 'downstream'
	elif(gene_location=='downstream'):
		return 'upstream'
	elif(gene_location=='peak_overlapping'):
		return 'gene_overlapping'
	else:
		assert gene_location=='gene_overlapping'
		return 'peak_overlapping'
"""
def closest_peaks(offset):
	peaks_dict={}
	file_open=open("filtered_final_"+str(offset)+".bed","r")
	lines=file_open.readlines()
	for line in lines:
		single_entry=line.split()
		peaks_dict[single_entry[5]]=[]
		if(len(single_entry)>7):
			peaks_dict[single_entry[7]]=[]

	for line in lines:
		single_entry=line.split()
		#opposite function was used for single entry 6
		peaks_dict[single_entry[5]].append([single_entry[3],single_entry[1],single_entry[2],single_entry[6]])
		if(len(single_entry)>7):
			#opposite function was used for single entry 9
			peaks_dict[single_entry[7]].append([single_entry[3],single_entry[1],single_entry[2],single_entry[8]])
	file_open.close()
	
	peaks_file=open("filtered_peaks_"+str(offset)+".txt","wb+")
	
	for key in peaks_dict:
		#print key,peaks_dict[key]
		peaks_file.write(key+"\t"+"\n\t".join(i[0]+"\t"+str(i[1])+"\t"+str(i[2])+"\t"+i[3] for i in peaks_dict[key])+"\n")
	peaks_file.close()
	
def closest_genes():
	if len(sys.argv) !=4:
		print "USAGE: python closest_genes.py <peaks file in bed format> <genes file in bed format> <offset>"
		exit()

	peaks=sys.argv[1]
	genes=sys.argv[2]
	offset=sys.argv[3]
	peaks_arr=np.genfromtxt(peaks, delimiter='\t',dtype=None)
	genes_arr=np.genfromtxt(genes, delimiter='\t',dtype=None)
	if not os.path.exists("filtered"+str(offset)+"_explained"+".txt") or os.path.getsize("filtered"+str(offset)+"_explained"+".txt") == 0:
		filter(peaks_arr,genes_arr,int(offset))
	
	filtered_peaks_arr=np.genfromtxt("filtered"+str(offset)+"_explained"+".txt",delimiter='\t',dtype=None)
	#gene start points and end points dictionary
	genes_dict={}
	
	for i in xrange(genes_arr.size):
		genes_dict[genes_arr[i][3]]=[]
	
	for i in xrange(genes_arr.size):
		genes_dict[genes_arr[i][3]]+=[genes_arr[i][1],genes_arr[i][2]]
	
	genes_dist_dict={}
	for i in xrange(len(filtered_peaks_arr)):
		genes_dist_dict[filtered_peaks_arr[i][3]]=[]

	for i in xrange(len(filtered_peaks_arr)):
		genes_dist_dict[filtered_peaks_arr[i][3]].append([filtered_peaks_arr[i][5],filtered_peaks_arr[i][7]])
	
	for key in genes_dist_dict:
		genes_dist_dict[key]=sorted(genes_dist_dict[key],cmp=mycmp)
	
	
	final=open("filtered_final_"+str(offset)+".bed","wb+")
	for i in xrange(filtered_peaks_arr.size):
		final.write(filtered_peaks_arr[i][0]+"\t"+str(filtered_peaks_arr[i][1])+"\t"+str(filtered_peaks_arr[i][2])+"\t"+str(filtered_peaks_arr[i][3])+"\t"+str(filtered_peaks_arr[i][4])\
			+"\t"+"\t".join(j[0]+"\t"+str(j[1]) for j in genes_dist_dict[filtered_peaks_arr[i][3]])+"\n")
	final.close()
	subprocess.call(["uniq","filtered_final_"+str(offset)+".bed","filtered_final1_"+str(offset)+".bed"])
	subprocess.call(["rm","filtered_final_"+str(offset)+".bed"])
	subprocess.call(["mv","filtered_final1_"+str(offset)+".bed","filtered_final_"+str(offset)+".bed"])
	value=raw_input("Do you want the peaks list also?  (Press y or n)\n")
	if value.lower() in {"y", "yes", "yea", "yeah" ,"si", "go", "aye", "sure"}:
		closest_peaks(offset)
closest_genes()