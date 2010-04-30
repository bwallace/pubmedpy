from  author_net_functions import *
import pickle

#perform a dummy search to show how to make the author networks
#my_terms="multivariate meta-analysis AND Stat Med[so]"
my_terms="JAK2"
my_xml_file = "mini.xml"

my_graph_file = "meta2.dot"
#download_from_pubmed(my_terms, my_xml_file, 200)

parsed_xml = parse_pubmed_xml(my_xml_file)
#author_network(parsed_xml[0], parsed_xml[1], my_graph_file)

paper_mat = create_paper_adj_array(parsed_xml[1])
author_mat_tuple = create_author_adj_array(parsed_xml[1])

outfile = open("lala.raw", "wb")

pickle.dump([paper_mat, author_mat_tuple], outfile)

outfile.close()

