import pubmedpy 

#
# A demonstration of using the pubmedpy library. Here
# we have some pubmed ids and labels in the `cohen.tsv' file.
# For demonstrate purposes, the meaning of these labels
# is irrelevant.
#

# hard coded indices for accessing fields
DRUG = 0
PMID = 2
LABEL = 3

drug_name = "ACEInhibitors"

# example of using pubmedpy
def get_pmid_to_label_dict():
    pmid_to_labels = {}
    #
    # cohen contains a list of citations; if there's an
    # 'I' here, the correct label is '1'; else '-1'. (Again
    # the meaning of this is immaterial here)
    #
    cohen = [c.replace("\n", "").split("\t") for c in open("cohen.tsv", "r").readlines()]
    drugs = list(set([c[DRUG] for c in cohen]))
    for abstract in [c for c in cohen if c[DRUG] == drug_name]:
        pmid = abstract[PMID]
        pmid_to_labels[pmid] = 1.0 if abstract[LABEL] == 'I' else -1.0
    return pmid_to_labels
    
pubmedpy.set_email("byron.wallace@gmail.com")
lbl_dict = get_pmid_to_label_dict()
pubmedpy.fetch_and_encode(lbl_dict.keys(), "output", lbl_dict = lbl_dict)

    