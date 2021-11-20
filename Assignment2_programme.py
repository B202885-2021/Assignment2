#!/usr/local/bin/python3

import os, subprocess

# First we retreive the relevant sequences protein sequences, given A) Protein name (we have used glucose-6-phosphatase, but change that portion of code to whatever protein you want to investigate) and B) Taxonomic group (we have used Aves, but change to whatever protein you want to investigate).

# Task 1: get the data

def input_info(prot, taxgroup) :
    import string
    print("You have provided the following details: \n \tProtein: ", prot, "\n\tTaxonomix group: ", taxgroup)

details = {}
details["Protein"] = input("Which protein family would you like to search? ")
details["Taxonomic group"] = input("Which taxonomic subset would you like to search? ")

input_info(*list(details.values()))
input_vals = list(details.values())

os.system("esearch -db protein -query '" + input_vals[0] + " [PROT]' | efilter -query '" + input_vals[1]+ "' | efetch -db protein -format fasta >> raw_seq_data.fasta")

# check everything is there:

os.system("cat raw_seq_data.fasta")

# Task 2: Check the sequence data and output 1) how many sequences there are, 2) what species/# sequences per species and 3) "Do you wish to continue?"

raw_seq_data = open("raw_seq_data.fasta").read()
raw_seq_data_array = raw_seq_data.split()

    # First count how many different sequences there are:

lines = subprocess.check_output("grep '>' raw_seq_data.fasta | wc -l", shell=True, text=True).rstrip()

print("There are " + lines + " different protein sequences in this set.")

# Create file with sequence info (FASTA header only)

os.system("grep '>' raw_seq_data.fasta >> sequence_info.txt")

print(" Removing partial sequences...")

partial_lines = subprocess.check_output("grep 'partial' sequence_info.txt | wc -l", shell=True, text=True).rstrip()

if int(partial_lines) > 0 :
    print("There are " + partial_lines + " partial sequences, which will be removed.")
    os.system("grep -v 'partial' sequence_info.txt >> new_seq_info.txt")
    lines = subprocess.check_output("cat new_seq_info.txt | wc -l", shell=True, text=True).rstrip() # We've now made 'lines' the number of FULL proteins.
    
print(" 'lines' has been redefined to be the number of full (non-partial) protein sequences. You have " + lines + " full sequences in your set.")

if int(lines) > 1000 :
    def check_nr(ans) : 
        import string
        print(" You have provided the following answer: ", ans)
        if ans == "y" :
            print("This might make the code very slow.")
        else:
            print("Please run the programme again, choosing a different protein or narrowing your taxonomix group.")
            exit()
    
check = {}
check["Answer"] = input("There are more than 1000 sequences in your set. Do you still wish to continue? (y/n)")
check_nr(*list(check.values())) 

# split predicted and existing

os.system("grep 'PREDICTED:' new_seq_info.txt >> predicted_seq_info.txt")
    
os.system("grep -v 'PREDICTED:' new_seq_info.txt >> existing_new_seq_info.txt")

    # Now output the species and # sequences per species
    
os.system("sed 's/>//' predicted_seq_info.txt >> new_predicted.txt")

os.system("sed 's/>//' existing_new_seq_info.txt >> new_existing.txt")

os.system('''awk 'BEGIN{FS=" "} {print $1}' new_predicted.txt >> accession_numbers.txt''') # list of accession numbers

os.system('''awk 'BEGIN{FS=" "} {print $4,$5}' new_predicted.txt >> species_list.txt''')

os.system('''awk 'BEGIN{FS=" "} {print $1}' new_existing.txt >> accession_numbers.txt''') # list of accession numbers

os.system('''awk 'BEGIN{FS=" "} {print $3,$4}' new_existing.txt >> species_list.txt''')

species_names = open("species_list.txt").read().split() # create vector of species names

# need to clean this up:

import re

names_clean = []

for name in species_names :
    if re.search(r'^\[', name) :
        names_clean.append(name[1:].upper())
    if re.search(r'\]$', name):
        names_clean.append(name[:-1].upper())

fullnames_clean = [ ' '.join(x) for x in zip(names_clean[0::2], names_clean[1::2]) ] # array of full names of species
        
names_dict = {i:fullnames_clean.count(i) for i in fullnames_clean} # count the number of names and output as a dictionary

import operator

sorted_d = sorted(names_dict.items(), key=operator.itemgetter(1), reverse=True)

print("Below are the species names and how many protein sequences exist for each species.") 

for i in sorted_d : 
    print(i)

print("There are " + str(len(sorted_d)) + " species in total.")

if len(sorted_d) < 15:
    def small_set(ans) : 
        if ans == "y" :
            print("This will result in a less diverse phylogeny, which may lead to less significant results.")
        else:
            print("Please run the programme again, choosing a different protein or expanding your taxonomix group.")
            exit()
    check = {}
    check["Answer"] = input("There are less than 15 species in your set. Do you still wish to continue? (y/n)")
    check_nr(*list(check.values())) 

# We now plot sequence conservation

os.system("clustalo -i raw_seq_data.fasta -o aligned.fa") 

os.system("prettyplot -blocksperline=5 -boxcol -consensus -ratio=0.59 -graph ps -sequence aligned.fa")

# Pretty graph of alignment, we can create a pdf file (though this is generally unnecessary, can just view the .ps file)

os.system("gs -dSAFER -dBATCH prettyplot.ps &")

# Now plot the sequence conservation using plotcon

os.system("plotcon -sequences aligned.fa -winsize 4 -graph cps")

os.system("gs -dSAFER -dBATCH plotcon.ps") 

# Task 3: Now find the conserved domains by comparing the full protein sequences against PROSITE database

# We first identify all the sequences that are full, and put them into a fasta file

raw_seq_data = open("raw_seq_data.fasta").read()
raw_seq_splitfasta = raw_seq_data.split('>')

partial_seqs_fasta = []
full_seqs_fasta = []

for seq in raw_seq_splitfasta :
    if re.search(r"partial", seq) :
        partial_seqs_fasta.append(seq)
    else :
        full_seqs_fasta.append(seq)

# Write each full fasta sequence to a separate file

for num, element in enumerate(full_seqs_fasta) :
    with open(f"{num}.fasta", 'w') as f:
        f.write(element)
    
# NOTE: the first file 0.fasta has no content. For future ref, start as 1.fasta

# We need to reinsert the '>'

for num, element in enumerate(full_seqs_fasta):
    subprocess.run(["sed -i '1s/^/>/' " + str(num) + ".fasta"], shell=True)

# Now run patmatmotifs on everything

for num in range(1, len(full_seqs_fasta)) :
    subprocess.run(["patmatmotifs " + str(num) + ".fasta output_" + str(num)], shell=True)

# This gives the conserved regions for each sequence

os.system("grep 'Motif' output_* >> present_motifs.txt")

os.system("cut -f 3- -d ' ' present_motifs.txt >> motifs_clean.txt")

motifs = open("motifs_clean.txt").read().split()

print("Below is a list of the different motifs and the number of occurences of each motif in the sequence set.")

motif_count = {i:motifs.count(i) for i in motifs}

motif_count

# Task 4: for the last portion, we identify double helices using helixturnhelix
# In order to do this, we need to get rid of special character X. 

for num in range(1, len(full_seqs_fasta)) :
    subprocess.call(["grep -v '>' " + str(num) + ".fasta >> " + str(num) + "_sequence_only.fa"], shell=True)

# identify and remove files with X in them (ie with unknow characters)

a_list = list(range(1, len(full_seqs_fasta)))
to_remove = []
to_keep = []

for index, item in enumerate(a_list) :
    item_open = open(str(item)+"_sequence_only.fa").read()
    if re.search(r"X", item_open) :
        to_remove.append(item)
    else :
        to_keep.append(item)
        
for item in to_remove :
    os.system("rm -f " + str(item) +"_sequence_only.fa")

# Do helixturnhelix on it and output hitcounts per sequence (output species)

for index, item in enumerate(to_keep) :
    os.system("helixturnhelix -sequence " + str(item) + ".fasta -outfile " + str(item) +"_helix.txt")
    
# Now we search for helix numbers

os.system('''grep 'HitCount:' *_helix.txt | awk 'BEGIN{FS = " "} {print $3}' >> number_of_hits.txt''') 

number_of_helices = open("number_of_hits.txt").read().split()

numbers = []

for index, item in enumerate(number_of_helices) :
    numbers.append(int(item))

print("There are " + str(sum(numbers)) + " double helices in this sequence set.")


    
    
