import pubmedpy

DRUG  = 0
PMID = 2
LEVEL1 = 3
LEVEL2 = 4

def download_cohen():
    ''' Process the cohen dataset. '''
    pubmed_fetchr.set_email("bwallace@tuftsmedicalcenter.org")
    cohen = open("COHEN.tsv", "r").readlines()
    cohen = [c.replace("\n", "").split("\t") for c in cohen]
    drugs = list(set([c[DRUG] for c in cohen]))
    pmids = {}
    for d in drugs:
        cur_ids = [c[PMID] for c in cohen if c[DRUG] == d]
        pmids[d] = cur_ids

    for d in drugs:
        print "on drug: %s" % d
        records = pubmed_fetchr.batch_fetch(pmids[d])
        pubmed_fetchr.write_out_fields(["AB", "TI", "RN", "MH", "PT", "AU"], records, d)
    
    print "done!"
    

def get_labels():
    drugs_to_lbls = {}
    cohen = open("cohen.tsv", "r").readlines()
    cohen = [c.replace("\n", "").split("\t") for c in cohen]
    drugs = list(set([c[DRUG] for c in cohen]))
    for d in drugs:
        level1, level2 = {}, {}
        for abstract in [c for c in cohen if c[DRUG] == d]:
            pmid = abstract[PMID]
            level1[pmid] = 1.0 if abstract[LEVEL1] == 'I' else -1.0
            level2[pmid] = 1.0 if abstract[LEVEL2] == 'I' else -1.0
        drugs_to_lbls[d] = (level1, level2)
    return drugs_to_lbls