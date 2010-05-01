from my_functions import * 

my_search = 'Nomogram can help estimate risk of serious hyperbilirubinemia in healthy infants'
pmids = get_pmid(my_search, 2000)

titles = read_isi_file("test.txt")

ids = fetch_pmids_from_titles(titles)

f = open("pmids.pickle", "rb+")
pickle.dump(ids,f)

f.close()

