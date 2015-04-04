#!/usr/bin/python
import sys
import time
def count_reads(reads_list,peaks_list):
	count=0
	for i in xrange(len(reads_list)):
		# start point of the read is in between the peaks range
		if (peaks_list[0]<=reads_list[i][0]<=peaks_list[1]):
			count+=1
		# end point of the read is in between the peaks range
		elif(peaks_list[0]<=reads_list[i][1]<=peaks_list[1]):
			count+=1
		# read itself is in the range of the peak
		elif(reads_list[i][0]<=peaks_list[0]<=reads_list[i][1]):
			count+=1
	return count
	
def reads_mapped():
	if(len(sys.argv)!=3):
		print "USAGE: python read_mapped.py <peaks_bed_file_name> <filtered_bed_file_name>"
		exit()
	print "Processing peaks ..."
	peaks_bed_file=sys.argv[1]
	filtered_bed_file=sys.argv[2]
	peaks_dict={}
	file_open=open(filtered_bed_file,"r")
	lines=file_open.readlines()
	
	#filling the dictionary for filtered peaks file
	for line in lines:
		single_entry=line.split()
		peaks_dict[single_entry[3]]=[]
	for line in lines:
		single_entry=line.split()
		peaks_dict[single_entry[3]].append([int(single_entry[1]),int(single_entry[2])])
	
	for line in lines:
		single_entry=line.split()
		peaks_dict[single_entry[3]].append([0])
	
	file_open.close()
	#print peaks_dict['Scaffold143.475']
	#start_time=time.time()
	print "Processing reads ..."
	#filling the dictionary for peaks file
	file_open=open(peaks_bed_file,"r")
	reads_dict={}
	lines=file_open.readlines()
	count=len(lines)
	for line in lines:
		single_entry=line.split()
		reads_dict[single_entry[0]]=[]

	for line in lines:
		single_entry=line.split()
		reads_dict[single_entry[0]].append([int(single_entry[1]),int(single_entry[2])])
	file_open.close()
	start_time=time.time()/60
	print "Counting the reads mapped ..."
	for k in peaks_dict:
		peaks_dict[k][1][0]=count_reads(reads_dict['Scaffold'+str(int(float(k[8:])))],peaks_dict[k][0])
	#print peaks_dict['Scaffold10.274']
	print "Writing the count of the reads mapped into the file reads_count.txt..."
	file_open=open("reads_count"+"_"+filtered_bed_file[:-4]+".txt","wb+")
	file_open.write("Scafffold Name"+"\t"+"Peak ID"+"\t"+"Peak_start"+"\t"+"Peak_end"+"\t"+"reads_mapped"+"\n")
	total_count=0
	for key in peaks_dict:
		file_open.write('Scaffold'+str(int(float(key[8:])))+"\t"+key+"\t"+str(peaks_dict[key][0][0])+"\t"+str(peaks_dict[key][0][1])+"\t"+str(peaks_dict[key][1][0])+"\n")
		total_count+=peaks_dict[key][1][0]
	print "Total peaks matched:",total_count
	print "Total reads in the file "+peaks_bed_file+" are :",count
	print "The ratio of reads mapped to peaks and total reads:",float(total_count)/count
	print "Total time for counting reads:",time.time()/60-start_time
	file_open.close()
reads_mapped()