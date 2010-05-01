from my_functions import * 

#my_search = 'Nomogram can help estimate risk of serious hyperbilirubinemia in healthy infants'
#pmids = get_pmid(my_search, 2000)

titles = read_isi_file("test.txt")


f = open("pmids.pickle", "wb+")
ids = fetch_pmids_from_titles(titles)
pickle.dump(ids,f)
f.close()

