# Sequence-Analysis-Scripts
Some Python scripts that I wrote for DNase-seq data analysis in summer for my research.

Closest_genes --> Resembles the BEDOPS 'closest-features' command, but works better than the BEDOPS command. closest-features provides all the features in that Scaffold whether or not it is close to the feature of interest. While my script asks the user a distance parameter, and finds all the closest features that lie within that distance from the border of the feature of interest.

extract_sequence --> Extracts the sequence from the user inputted fasta file at the positions determined by a user inputted bed file.

npf-to-bed --> Removes last 3 columns of the npf file and saves the resulting bed file.

read_mapped --> Maps a peak to the number of reads that lie within that peak and gives a confidence level for the peaks.