#!/usr/bin/env python

import re, os
import getopt, sys
import FjrParser as fp

def usage (prog_name):
    print ">>> usage: " + prog_name + " -c <crab directory> [--help | -h] [--quite | -q]"

def help (prog_name):
    usage(prog_name)

    print """
    -c\t\t\t (Mandatory) CRAB project directory
    --quiet, -q\t\t Wehn turned on, works silently
    --help, -h\t\t Print this message
    """
    sys.exit(1)

if __name__ == '__main__':
  (prog_path, prog_name) = os.path.split(sys.argv[0])
  if len(sys.argv) < 2:
      help(prog_name)

  try:
      opts, argv = getopt.getopt(sys.argv[1:], "hc:q", ["help","quiet"])
  except getopt.GetoptError, err:
      print ">>> " + str(err)
      help(prog_name)

  quiet = False
  directory = ""
  for o, a in opts:
    if o == "-c":
      directory = a
    elif o in ("-q", "--quiet"):
      quiet = True
    elif o in ("-h", "--help"):
      help(prog_name)
    else:
      print ">>> ", "unhandled option", o
      help(prog_name)

  if directory == "" or not os.path.exists(directory):
    print ">>> path <" +  directory + "> may not exist"
    help(prog_name)

  fjrs = fp.get_fjrs(directory)
  localPath = None
  expList = []
  if not quiet:
    print 'This is the expected list:'
  for f in fjrs:
    #if '74' not in f: continue
    if not quiet:
      print ">>> ", "processing fjr", f

    # parse fjr
    try:
      doc = fp.parse_xml(f)
    except:
      if not quiet:
        print ">>> ", 'FrameworkJobReport', f, 'could not be parsed and is skipped'
      continue

    # if the job finished successfully get the PFN
    # check to make sure that file actually exists in your storage space (set for EOS)
    #try:
    if fp.is_goodfile(doc) :
      (lfn, pfn, surl) = fp.get_filenames(doc)
      localFile = '/'+'/'.join(pfn.split('/')[6:])
      if not quiet:
        print localFile
      if localPath == None:
        localPath = '/'+'/'.join(pfn.split('/')[6:-1])
      if not os.path.isfile(localFile):
        print pfn, 'does not actually exist'
        raw_input('resubmit this job, something is wrong')
      expList.append(pfn.split('/')[-1])

  #except:
    #if not quiet:
      #print ">>> ", "skipping fjr", f


  if not quiet:
    print 'This is the observed list:'
  obsList = None
  for dirname, dirnames, filenames in os.walk(localPath):
    obsList = set(filenames)
    # print path to all subdirectories first.
    #for subdirname in dirnames:
        #print os.path.join(dirname, subdirname)

    # print path to all filenames.
    for filename in filenames:
      if not quiet:
        print os.path.join(dirname, filename)


  expList = set(expList)
  #print obsList, expList
  # is obsList in expList?
  #if obsList <= expList: print 'obs is a subset of exp'
  # is expList in obsList?
  #if obsList >= expList: print 'exp is a subset of obs'

  #symmetric diff
  if len(obsList ^ expList) != 0:
    print 'these are missing:',expList-obsList
    print 'these are extra:',obsList-expList

