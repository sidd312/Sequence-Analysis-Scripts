#!/usr/bin/python
import sys
def npftobed():
	if len(sys.argv)!=2:
		print "USAGE: python filter.py <file_name>"
		exit()
	file_name=sys.argv[1]
	print "Creating bed files..."
	file_open=open(file_name,"r")
	bed_file=open(str(file_name[:-3])+"bed","wb+")
	lines=file_open.readlines()
	for line in lines:
		single_entry=line.split()
		bed_file.write(single_entry[0]+"\t"+single_entry[1]+"\t"+single_entry[2]+"\t"+single_entry[3]+"\t"+single_entry[6]+"\n")
    
npftobed()