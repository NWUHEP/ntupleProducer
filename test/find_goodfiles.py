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
  for f in fjrs:
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
    try:
      if fp.is_goodfile(doc) : 
        (lfn, pfn, surl) = fp.get_filenames(doc)
        print pfn
    except:
      if not quiet:
        print ">>> ", "skipping fjr", f
