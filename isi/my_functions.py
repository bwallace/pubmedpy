import os
import pickle
from numpy import * 
from Bio import Entrez
# Tell Entrez who you are
Entrez.email = "ttrikalin@gmail.com"


def get_pmid(search_term, ret_number):
    """Downloads up to 2000 abstracts as an xml and saves it"""
    handle = Entrez.esearch(db="pubmed", term=search_term, retmode="xml", 
        RetMax=ret_number)
    
    #record the id_list and then fetch the file into a local xml
    summary_record = Entrez.read(handle)
    pmid_list = summary_record["IdList"]

    return pmid_list

    

def batch_fetch(article_ids, batch_size=100):
    all_records = []
    total = len(article_ids)
    fetched_so_far = 0
    while fetched_so_far < total:
        records = fetch_articles(article_ids[fetched_so_far:fetched_so_far+batch_size])
        fetched_so_far += batch_size
        all_records.extend([r for r in records])
    return all_records
  
  
def read_isi_file(filename):
    """Read an ISI filename and store the titles in a list"""
    infile = open(filename, "r")

    i=0
    titles=[]
    for line in infile:
        l=line.strip().split("\t")
        if (i==0):
            dict = {}
            for j, word in enumerate(l):
                dict[l[j]]=j
        else:
            titles.append(l[dict["TI"]])
        i+=1

    return titles 


def fetch_pmids_from_titles(titles):
    '''Fetch pmids from possible titles from ISI'''
    pmids=[]
    clean_titles=[]
    for ti in titles:
        id = get_pmid(ti,1)
        if (id!=[]):
            clean_titles.append(ti)
            pmids.extend(id)
    
    return [pmids, clean_titles]

