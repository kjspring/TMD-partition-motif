# PubMed Crawler
# Author: Kevin Spring
# Date: 20 March 2014
# Language: Python
# Description: This script searches/crawls the PubMed database
#              for publications that have protein and lipid
#              raft associations.
##############################################################

from scrapy.item import Item, Field
# http://doc.scrapy.org/en/latest/intro/tutorial.html
from Bio.Medline import PubMed

# Can use XML to parse the abstract
# ex. http://www.ncbi.nlm.nih.gov/pubmed/24643062?report=xml&format=text
# <ABSTRACT>
# <PMID>
# can't get the list of references

# Set the seed URLs
# 1. Using the list of raft motif containing proteins to get a 
#    set of references.
# 2. Go through each of those references and save the citations
#    as future links to search.
# 3. Search for articles that have cited the seed references. 
#    Determine that they have term "lipid raft" and a protein
#    name.

# Define the data you want to scrape

# Write a spider to extract the data
# 1. Define the start URL
# 2. Define the rules for following the links
# 3. Define the rules to extract the data from the pages.
#    Data we want.
#    Number of references that have both search terms
#    Reference PMID numbers
#    Total citations on that seed paper

# Store the data in an information:
#   Protein Name/ID
#   PMID ID
#   

