#!/usr/bin/env python
# coding: utf-8

from sourmash import signature
import pandas as pd
import glob
import os

sigs = glob.glob("outputs/gather_matches/*sig")

# read in all genome signatures that had gather matches for the var imp hashes
# create a dictionary, where the key is the genome and the values are the minhashes
genome_dict = {}
for sig in sigs:
    if os.path.getsize(sig) > 0:
        sigfp = open(sig, 'rt')
        siglist = list(signature.load_signatures(sigfp))
        loaded_sig = siglist[0] 
        mins = loaded_sig.minhash.get_mins() # Get the minhashes 
        genome_dict[sig] = mins

# read in vita variables
# CHANGE THIS TO READ IN ALL VITA VARS AND COMBINE 
#sigfp = open(str(input[1]), 'rt')
sigfp = open("outputs/vita_rf/vita_vars_merged.sig", 'rt')
vita_vars = sig = signature.load_one_signature(sigfp)
vita_vars = vita_vars.minhash.get_mins() 

# generate a list of all minhashes from all genomes
all_mins = []
for file in sigs:
    if os.path.getsize(file) > 0:
        sigfp = open(file, 'rt')
        siglist = list(signature.load_signatures(sigfp))
        loaded_sig = siglist[1]
        mins = loaded_sig.minhash.get_mins() # Get the minhashes 
        all_mins += mins

# get intersection size between hashes in vita and 
# hashes in genomes
# len(list(set(vita_vars) & set(all_mins)))

# define a function where if a hash is a value, return all key for which it is a value
def get_all_keys_if_value(dictionary, hash_query):
    genomes = list()
    for genome, v in dictionary.items():
        if hash_query in v:
            genomes.append(genome)
    return genomes

# create a dictionary where each vita_vars hash is a key, 
# and values are the genome signatures in which that hash
# is contained
vita_hash_dict = {}
for hashy in vita_vars:
    keys = get_all_keys_if_value(genome_dict, hashy)
    vita_hash_dict[hashy] = keys

# transform this dictionary into a dataframe and format the info nicely
df = pd.DataFrame(list(vita_hash_dict.values()), index = vita_hash_dict.keys())
df = df.reset_index()
df = pd.melt(df, id_vars=['index'], var_name= "drop", value_name='genome')
# remove tmp col drop
df = df.drop('drop', 1)
# drop duplicate rows in the df
df = df.drop_duplicates()
# write the dataframe to csv
df.to_csv("vita_hash_to_genome_mapping.csv", index = False)



