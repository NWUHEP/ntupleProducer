#!/usr/bin/env python
import os,sys
'''
This script assumes there is a list of files obtained with find_dupl.sh script.
See this twiki for details
https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideCrabFaq#Duplicated_jobs_and_duplicate_ou
'''

if len(sys.argv) == 3:
    path     = sys.argv[1]
    duplFile = sys.argv[2]
else:
    print 'Need to specify path and file with list of unwanted files'

file = open(duplFile)
for line in file:
    target = path + line
    print 'Deleting file {0}'.format(target)
    os.system('rm {0}'.format(target))
