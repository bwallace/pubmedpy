import os
import linecache, re
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
                try: 
                    initials = author_list[who]["Initials"]
                except KeyError: 
                	initials = "no initials for this person"
                new_name = lastname.lower() + "_" + initials.lower() 
        else: 
            newname = author_list[who]["CollectiveName"].lower()
        names.append(new_name)
    authors.append(names)

    # close the handle!!! 
    in_handle.close()
    return (ids, authors)


# this tells you if two authornames match 
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
    
# the comparison is now handled by a function 
#        if w1==w2:
#           same +=1
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

