import os
import linecache, re
from numpy import * 
from Bio import Entrez
# Tell Entrez who you are
Entrez.email = "ttrikalin@gmail.com"


def download_from_pubmed(search_term, filename, ret_number):
    """Downloads up to 2000 abstracts as an xml and saves it"""
    handle = Entrez.esearch(db="pubmed", term=search_term, retmode="xml", 
        RetMax=ret_number)
    
    #record the id_list and then fetch the file into a local xml
    summary_record = Entrez.read(handle)
    id_list = summary_record["IdList"]

    if not os.path.isfile(filename):
        print "Downloading..."
        net_handle = Entrez.efetch(db="pubmed", id=id_list, retmode="xml" )
        out_handle = open(filename, "w")
        out_handle.write(net_handle.read())
        out_handle.close()
        net_handle.close()
        print "File \"" + filename + "\" saved!\n"
    elif os.path.isfile(filename):
        string= "File \"" + filename + "\" exists!\n"
        print string
    return



def fix_the_dtd(filename):
    """Read the second line of the pubmed xml and if the dtd is 'pubmed_100301.dtd' change it to 
    'pubmed_100101.dtd'"""
    
    line = linecache.getline(filename,2) 
    if re.search("pubmed_100301\.dtd", line)==None:
        return filename
    else: 
        correctedline = re.sub("pubmed_100301\.dtd", "pubmed_100101.dtd", line)
        newfilename = "corrected_"+filename
        outhandle = open(newfilename, "w")
        contents = linecache.getlines(filename)
        contents[1]=correctedline
        for line in contents:
            outhandle.write(line)
        outhandle.close()
        print "I corrected the DTD from 'pubmed_100301.dtd' to 'pubmed_100101.dtd'"
        return newfilename



def parse_pubmed_xml(filename):
    """parse_pubmed_xml(filename): give a Pubmed xml and get back a list with the PMIDs 
       and a list of the lists of authors """
    
    filename = fix_the_dtd(filename)
    
    in_handle = open(filename, "r")
    citations = Entrez.parse(in_handle)

    ids=[]
    authors = []
    i=0
    for citation in citations:
        ids.append(citation["MedlineCitation"]["PMID"])
        names = []
        i +=1 
        #print i
        try: 
            author_list = citation["MedlineCitation"]["Article"]["AuthorList"]
        except KeyError: 
            #print "****"
            errorStr = "no authors in citation " +  str(i)
            author_list = [{"CollectiveName":errorStr}]
        for who in range(len(author_list)): 
            if author_list[who].keys()[0]=="LastName":
                lastname = author_list[who]["LastName"]
                initials = author_list[who]["Initials"]
                new_name = lastname.lower() + "_" + initials.lower() 
            else: 
                newname = author_list[who]["CollectiveName"].lower()
            names.append(new_name)
        authors.append(names)

    # close the handle!!! 
    in_handle.close()
    return (ids, authors)


def is_this_the_same_author(author1, author2):
    """Returns true if the "Surname IM" and the initials are identical - handles 
    unicode"""
    
    if (author1==author2) & (author1!=""):
        return (True, "Identical")
    elif (author1==author2) & (author1==""):
        return (False, "Both authornames==''")
    else:
        if (len(author1)<len(author2)):
            if author2.count(author1)==1:
                return (True, "Ignoring (some) initials")
        else: 
            if author1.count(author2)==1:
                return (True, "Ignoring (some) initials")
    
    # if nothing matches
    return (False, "No match")


# This tells you how many common elements are between two author lists 
def how_many_common(list1, list2):
    """how_many_common(list1, list2): give it two lists and it will 
    return the number of common elements """
    l1 = sorted(list1)
    l2 = sorted(list2)
    same = 0
    for w1 in l1:
        for w2 in l2:
            same += is_this_the_same_author(w1,w2)[0]
    return same

## The comparison ids now handled by a small function
#        if w1==w2:
#            same +=1
#        else:
#            if (len(w1) < len(w2)):
#            if w2.count(w1)==1:    # is one a substring (initials problem) 
#                same +=1
#        else: 
#            if w2.count(w1)==1:
#                same +=1
#           
#    return same


#now write the graph-viz file
def author_network(ids, authors, out_filename):
    """author_network(ids, authors, filename): Give it a list of Pubmed ids and 
       a list of lists of authors, and it will write the GraphViz code for 
       the author network in a file named filename"""

    out_handle = open(out_filename, "w")
    out_handle.write("Graph {\n")
    for first in range(0, len(ids)-2,1):
        for second in range(first+1,len(ids)-1,1):
            n=how_many_common(authors[first], authors[second])
        line = ''
        if (n>=1):
            line = '    "' + ids[first] + '" -- "' + ids[second] + '"\n'
        elif (second==first+1):
            line = '    "' + ids[first] + '"\n' 
        
        if line != "":
            out_handle.write(line)

    out_handle.write("}\n")
    out_handle.close()
    print "Graph file \"" + out_filename + "\" saved/replaced!"
    return



def create_paper_adj_array(authors):
    """This function takes a list of lists of authors and 
    returns an adjacency array of all the PMIDs """
    n = len(authors)
    
    pmid_mat = zeros((n,n))
    for row in range(n):
        for col in range(row, n, 1):
            #print row, col
            #print authors[col]
            pmid_mat[row,col]=how_many_common(authors[row], authors[col])
            if (col>row):
                pmid_mat[col,row] = pmid_mat[row,col]
    
    return pmid_mat
    
        


def augment_unique_author_list(unique_list, new_author):
    """This function compares a UNIQUE list with a new name. If the new name is not 
    included in the list, it adds it and returns the augmented list. If it was included in the list
    it just returns the list"""
    
    if (how_many_common(unique_list, [new_author])==0):
        return unique_list + [new_author]
    else:
        return unique_list 

def create_unique_author_list(authors):
    """This function takes a list of lists of authors and 
    returns a list of unique author names """
    
    n = len(authors)
    unique_author_list = authors[0]
    
    for i in range(1, n, 1):
        paper_names = authors[i]
        #print paper_names
        if(how_many_common(unique_author_list, paper_names)==0):
            unique_author_list += paper_names 
        else: 
        #print i 
            for author in paper_names: 
            #print unique_author_list
                unique_author_list=augment_unique_author_list(unique_author_list, author)

    return unique_author_list
    

def create_author_adj_array(authors):
    """This function takes a list of lists of authors and 
    returns an adjacency array of all the authors """
    
    unique_list = create_unique_author_list(authors)
    n= len(unique_list)

    author_mat = zeros((n,n))
    for row in range(n):
        for col in range(row, n, 1):
            how_many= 0 
            pair = [unique_list[row], unique_list[col]]
            for paper_authors in authors:
                if (how_many_common(paper_authors, pair) == 2):
                    how_many += 1
            author_mat[row,col]=how_many
            if (col>row):
                author_mat[col,row] = author_mat[row,col]
    return (author_mat, unique_list)
    


# this is Bython 

def fetch_articles(article_ids):
    print "Fetching abstracts..."
    handle = Entrez.efetch(db="pubmed",id=article_ids,rettype="medline",retmode="text")
    records = Medline.parse(handle)
    print "ok"
    return records   
    
def batch_fetch(article_ids, batch_size=100):
    all_records = []
    total = len(article_ids)
    fetched_so_far = 0
    while fetched_so_far < total:
        records = fetch_articles(article_ids[fetched_so_far:fetched_so_far+batch_size])
        fetched_so_far += batch_size
        all_records.extend([r for r in records])
    return all_records
  
  
def article_search(db, search_str):
    handle = Entrez.esearch(db=db,term=search_str)
    record = Entrez.read(handle)
    article_ids = record["IdList"]
    print "Found %s articles" % len(article_ids)
    return fetch_articles(article_ids)
  

