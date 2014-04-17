#!/usr/bin/env python
import os
'''
This script assumes there is a list of files obtained with find_dupl.sh script.
See this twiki for details
https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideCrabFaq#Duplicated_jobs_and_duplicate_ou
'''

duplFile = "duplicates.txt"

f = open(duplFile)
for l in f:
    print l,
    os.system("rm "+l)
