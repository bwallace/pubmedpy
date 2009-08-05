#!/usr/bin/env python
# encoding: utf-8
#
# 
# The MIT License
# 
# Copyright (c) 2009 Byron C Wallace
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# http://www.opensource.org/licenses/mit-license.php

"""
Byron C Wallace
Tufts Medical Center
pubmed_fetchr.py
--

Example use:

>python pubmed_fetchr.py -e bwallace@tuftsmedicalcenter.org -s biopython

(NCBI wants your email address for the web services stuff). -s is the search string, here we search 
for abstracts related to "biopython". 

"""

# std libraries
import sys
import getopt
import pdb
import os
from optparse import OptionParser
import Bio
from Bio import Entrez
from Bio import Medline

# home-grown
import tfidf2 

def main(argv=None):
  parser = OptionParser()
  parser.add_option("-e", "--email", action="store", type="string", dest="email")
  parser.add_option("-s", "--search_string", action="store", type="string", dest="search_str")
  (options, args) = parser.parse_args()
  
  Entrez.email = options.email
  records = article_search("pubmed", options.search_str)
  print "Writing out title and abstract data..."
  write_out_fields(["TI", "AB"], records)
  print("fin.")
    
def set_email(email):
    Entrez.email = email

def fetch_and_encode(article_ids, out_dir, binary_features=False, 
                                    labels=None, fields = ["AB", "TI"], out_f_name = ""):
    '''
    First fetches from the web, then encodes them.
    '''
    # first, fetch the articles
    fetch_and_write_out(article_ids, out_dir, fields = fields)
    
    for field in fields:   
        print "encoding %s..." %field
        # now, clean and encode them
        out_for_field = os.path.join(out_dir, field)
        tfidf2.encode_docs(out_for_field, os.path.join(out_for_field, "encoded"), out_f_name + field)
    print "finito."
    
def fetch_and_write_out(article_ids, base_out_dir, fields = ["AB", "TI", "RN", "MH", "PT", "AU"]):
    records = batch_fetch(article_ids)
    write_out_fields(fields, records, base_out_dir)
    
def fetch_articles(article_ids):
    print "Fetching abstracts..."
    handle = Entrez.efetch(db="pubmed",id=article_ids,rettype="medline",retmode="text")
    records = Medline.parse(handle)
    print "Done." 
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
  
 
def write_out_fields(field_keys, records, base_out_dir):
    # make a directory for each field
    os.mkdir(base_out_dir)
    for field in field_keys:
      os.mkdir(os.path.join(base_out_dir, field))
    # write out 
    for record in records:
      for field in field_keys:
        out_f = open(os.path.join(base_out_dir, field, record["PMID"]), "w")
        cur_field = None
        try:
          cur_field = record[field]
        except Exception, inst:
          cur_field = ""
          
        if isinstance(cur_field, list):
            cur_field = ", ".join(cur_field)
        out_f.write(cur_field)
        out_f.close()



if __name__ == "__main__":
  sys.exit(main())
