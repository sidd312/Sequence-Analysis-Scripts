#!/usr/bin/python
import sys
from Bio import SeqIO

def extract():
	if(len(sys.argv)!=4):
		print "USAGE: python extract_sequence.py <bed_file_name> <fasta_file_name> <offset>"
		exit()
	filehandle=open(sys.argv[2],"rb")
	print "Extracting sequences and making an output file ...."
	offset=sys.argv[3]
	record_dict=SeqIO.to_dict(SeqIO.parse(filehandle,"fasta"))
	filehandle.close()
	#Now extracting the starting and ending positions of the peaks and then slice the record.seq array, write into a file with the key as well.
	bed_file=open(sys.argv[1],"rb")
	output_file=open("Extracted_Sequence"+"_"+offset+".fa","wb+")
	lines=bed_file.readlines()
	for line in lines:
		single_entry=line.split()
		output_file.write(">"+single_entry[3]+"\n"+str(record_dict[single_entry[0]].seq[int(single_entry[1]):int(single_entry[2])])+"\n")
	bed_file.close()

extract()


