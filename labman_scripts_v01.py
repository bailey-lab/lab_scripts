#!/usr/bin/python3
VERSION=0.10_20210926

###standard modules
import argparse #parse initial arguements from command line
import configparser
import re    
import os
import os.path
#import pickle
import sys
#import subprocess
import textwrap
#import glob
#import copy
import csv

import collections

#these can potentially move into the subroutine
import configparser
import pyfiglet
import termcolor

#import random
#from collections import defaultdict, Counter, OrderedDict
#from collections import namedtuple


###PACKAGES REQUIRING INSTALLATION
#key science modules
#import numpy
#import scipy
#import numpy.string_
#import csv 
## key bio python modules 

###################################################################
###################################################################
#### FEATURE REQEUST     ##########################################
###################################################################
 
#make major subroutines contain specialized packages
#   particularlly those requireing installs
# rename bins --> bins/bags ### 

###################################################################
###################################################################
#### DEVELOPMENT HISTORY ##########################################
###################################################################

#2021-08-08 program started
#2021-09-26 stable version with all key features.
#2021-09-27 added color and better formatting
#2021-09-28 added bag # for storage and bag # for labeling instructions
#           fixed error for blank scan /enter only --> being a badcode

###################################################################
################################################################### 
### INITIALIZE GLOBAL VARIABLES ################################### 
################################################################### 

subcommands={} 

#add other global commands here. 

###################################################################
###################################################################
###################################################################  
###################################################################
###  MAIN        ##################################################
###################################################################


def main( args ):
  """Main handles general arguments and executing subcommands.
  
  Each subcommand launches a separate function. The pydoc subcommand 
  launches pydoc on this overall program file. 
  
  :param args: the main command line arguments passed minus subcommand
  """
  #print globals().keys()
  #print "ARGUMENTS", args
  
  if len(args)  == 0 or args[0] in ["h", "help", "-h", "--h", "--help","-help"] :
    verbosity= 'shortDesc'
    if args[0] in ["help" , "--help", "-help"]:
      verbosity = 'longDesc' 
    program_name=os.path.basename(__file__)
    print ("VERSION: ", VERSION)
    print ("USAGE:",program_name, "[-h] subcommand [suboptions]")
    print ("DESCRIPTION: various scripts complementing MIPTools pipelines")
    print ("SUBCOMMANDS:")
    #tw=TextWrap()
    for k in subcommands.keys():
      text=subcommands[k][verbosity]
      text= textwrap.dedent(text)
      if text:
        text =  "%s:   %s " %(k, text )
        print (textwrap.fill(text,77, initial_indent='', subsequent_indent='         ') )
    print ("HELP:")
    print ("pydoc      detailed documentation of program structure")
    print ("-h/-help   short / long  subcommand descriptions")
    print ("For specific options:",program_name,"[subcommand] --help")
  elif args[0] == '--pydoc':
    os.system( "pydoc " + os.path.abspath(__file__) )
  elif args[0] in subcommands.keys():
    globals()[args[0]](args[1:])
  else:
    print ("unknown subcommand (" + args[0] + ") use -h for list of subcommands!")
    sys.exit(-1)
  sys.exit(0) #normal exit

#------------------------------------------------------------------------------
###############################################################################
####  SCAN SAMPLE BARCODES AND BIN INTO BAGS   ################################
###############################################################################
#------------------------------------------------------------------------------

shortDescText="scan barcoded samples  into bins  -- ordering chaos"
longDescText="""loads a tab-delimitde and config file controlling scanning and sorting of barcoded samples into bins """
subcommands['barcodescanbinner'] = { 'shortDesc':shortDescText, 'longDesc': longDescText }

 
def barcodescanbinner(args):
  """commandline routine for scanning barcodes and binning samples appropriate
  
  This command will take a list of samples/items (one per row) that contain a barcode and designated bin.   
  
  """
  ##### add global variables imports here ####
  #config parser 
  #https://stackoverflow.com/questions/11990556/how-to-make-global-imports-from-a-function 
  
  aparser=argparse.ArgumentParser( prog = os.path.basename(__file__) + " barcode_binner", 
      description= "set maximum number of reads at given position within the bam",
      epilog= "Note: to avoid errors most data is read in from a fixed sheet",
      formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )

  aparser.add_argument("--samplefile", required=True,  help='tab delimited table with samples and bins [BARCODE,BIN,SAMPLEID]')
  aparser.add_argument("--configfile", required=True, help='file with bins and key settings')
  aparser.add_argument("--skipbinningcheck", required=False, action='store_true',  help='skip checking if sample bins map to config bins')
  aparser.add_argument("--skipscancheck", required=False, action='store_true',  help='scan all bag label barcodes to make sure all function')
  aparser.add_argument("--validateonly", required=False, action='store_true',  help='stop execution at the validation stage')
  args=aparser.parse_args(args=args)

 
  ##############################################
  ### SETUP AND LOADING 
  ##############################################
  
  ### load the configuration file ###
  reverseheader ("LOAD CONFIGFILE", "digital")
  config= configparser.ConfigParser()
  config.optionxform = str #this forces mixed case in variable names
  config.read(args.configfile)
  print ( config.sections() )
  for sectionname in ["BINS", "REQUIRED"]:
    if not  sectionname in config.sections():
      print ( f'ERROR: CONFIGFILE missing SECTION: {sectionname}', file=sys.stderr  )
      exit(1)
  ###### check for required fields that will be free text##
  for requiredfield in ["FILEOUTSAMPLE","FILEOUTBAD","BARCODEMATCH", 'BAGSETNAME' ] :
    if not requiredfield in config['REQUIRED']:
      print ( f"ERROR: Missing [REQUIRED] field variable {requiredfield}= ????", file = sys.stderr)
      exit(1)
  ##### check for required fiels that will be free text and 
  for requiredfield in ["BAGNUM" , 'BAGSIZE'] :
    if not requiredfield in config['REQUIRED']:
      print ( f"ERROR: Missing [REQUIRED] field variable {requiredfield}= ????", file = sys.stderr)
      exit(1)
    if  config['REQUIRED'][requiredfield].isdigit()  ==False :
      print ( f"ERROR: [REQUIRED]  {requiredfield} =", config['REQUIRED'][requiredfield]
           , " must be an integer!", file = sys.stderr)
      exit(1)
  
  ######################################## 
  ##### create bag status dictionary  ###
  bagstatus={}
  binlist= [x for x in config['BINS'] ]
  print ("# Creating bags for BADBIN BADCODE  and user bins: ", binlist)
  for bag in   binlist + ["BADBIN", "BADCODE"]:
    bagstatus[bag]={}
    bagstatus[bag]['BAGSIZE']= config['REQUIRED'].getint("BAGSIZE")
    bagstatus[bag]['BAGNUM']=config['REQUIRED'].getint("BAGNUM")
    bagstatus[bag]['BAGSAMPLENUM']=0
  for section in ['BAGSIZE','BAGNUM']:
    if section in config.sections():
      for bag in config[section]:
        if  config[section][bag].isdigit()  == False :
          print ( f"ERROR: {section} {bag} =", config[section][bag]
             , " must be an integer!", file = sys.stderr)
          exit(1)
        bagstatus[bag][section]= config[section].getint(bag)
        
  for key in bagstatus:
    print (key, bagstatus[key])
  
  ###########################################3
  ### load the sample file ###
  reverseheader ("LOAD SAMPLEFILE", "digital")
  print ("#Loading and validating sample file format", args.samplefile)
  print ("## required column headers are BARCODE BIN and SAMPLEID")
  print ("## output file will have additional BAG column")
  with open(args.samplefile) as samplefile:
    ### load the table into memory for fast retrieval ###
    reader= csv.reader( samplefile, delimiter ="\t")
    samples = list ( reader)
  header={}
  print ("## Loaded ", len(samples)-1, "rows of samples ")
  for pos, name in enumerate( samples[0] ):
    if name == "":
      print ( f"ERROR: Sample file column pos is header is empty", file = sys.stderr)
      exit(1)
    if name in header:
      print (f"ERROR: Sample file header name {name} exists twice", file=sys.stderr)
      exit(1)
    header[name]=pos 
  for colname in ["BARCODE",  "BIN", "SAMPLEID"]:
    if not colname in header:
      print( f"SAMPLEFILE missing required column {colname})", file=sys.stderr)
      exit(1)   

  ##################################################
  ###### generate fullbarcodes and matched barcodes dictionaries
  ### full  = sample sheet column barcodes
  ### matched = pattern found using BARCODEMATCH variable in config file
  reverseheader ("EXTRACT BARCODEMATCH", "digital")
  fullbarcodes={}
  matchedbarcodes={}
  matchedbarcodecounts=collections.Counter()
  for  pos, row in enumerate (samples):
    if pos == 0:
      continue
    if len (row) !=  len(samples[0]):
      print ("ERROR: Length of row",  pos+1,"is", len(row)
          ,"and differs from header row", len(samples[0]), file=sys.stderr)
      exit(1)
    bc=row[ header['BARCODE'] ]
    thematch= re.search( config['REQUIRED']['BARCODEMATCH'], bc)
    if thematch:
      #print (thematch)
      fullbarcodes[bc]=pos
      matchedbarcodecounts.update({ thematch.group(1) :1})
      if thematch.group(1) in matchedbarcodes:
        matchedbarcodes[thematch.group(1)].append(bc)
      else:
        matchedbarcodes[thematch.group(1)]=[bc]
    else:
      print( "BARCODEMATCH ", config['REQUIRED']['BARCODEMATCH'], " has no match in barcode ", bc, "in row ", pos,  file=sys.stderr)
      exit(1)
  print ("MOST COMMON BARCODES (UNIQUE=1)\n", matchedbarcodecounts.most_common(50) )
  
  ########################################################################
  ### check of bins and what might fall in them in terms of barcodes/bins    
  reverseheader ("LOAD BINS", "digital")   
  print ("User defined bins",binlist)
 
  reverseheader ("CHECK SAMPLE MAPPING", "digital")   
  if args.skipbinningcheck:
    print ("#SKIPPING binning check by cross tabulation of samples and bins")
  else:
    print ("#VALIDATION crosstabulating and counting samples vs bins")
    samplebincounts=collections.Counter()
    finalbincounts=collections.Counter()
    #looop through all the samples except the header
    for samp in samples[1:]:
      bincode=samp[header['BIN']]
      samplebincounts.update({bincode:1})
      finalbin=lookup_bin(config['BINS'], bincode)
      finalbincounts.update({finalbin:1})
    print ("##SAMPLE BIN CATEGORIES from sample file:\n", samplebincounts, 
            " ===>   TOTAL SUM ", sum(samplebincounts.values())  )
    print ("#FINAL BAGGING BASED ON [BIN]\n", finalbincounts
              , " ===>   TTOTAL SUM ", sum(finalbincounts.values()) )

  ### check of barcode bin labels by scanning
  reverseheader ("CHECK LABEL SCAN", "digital")   
  if args.skipscancheck:
    print ("#SKIPPING SCAN CHECK of bin labels")
  else:
    print ("#INTERACTIVE  SCAN  CHECK OF BIN LABEL BARCODES")
    for label in binlist + [ "BADBIN","BADCODE",  "SKIP/REDO" ]:
      validinput=[label, '!N', '!NEXT', '!EXIT'] 
      if label in config['BINS']:   validinput[0] =  label  + label 
      response= get_validresponse (validinput , f'Scan {label}'
           , f"Scan barcode label({label}) or type !N(ext) or !E(end) =>",'cyan')
      if  response in ['!NEXT', '!N']:
        print ("Proceeding to next label...")
        continue
      if  response == '!EXIT':
        print ("Exiting label check ...")
        break
      print ("... good!")
        
        
  ########################################
  #setup output files ##
  ### fixed additional fields in output defined above
  reverseheader ("OUTPUT FILES", "digital")   
  outputcolumns=["OUTBINBAG","BAGNUM","BAGSAMPLENUM", "SCANNEDBARCODE", "BAGSETNAME"]
  if not os.path.isfile ( os.getcwd()+"/"+config['REQUIRED']['FILEOUTSAMPLE']) :
     print ("#Creating output file ",config['REQUIRED']['FILEOUTSAMPLE'] )
     os.system('echo "'+ "\t".join(samples[0]+outputcolumns) 
              +'" > ' + config['REQUIRED']['FILEOUTSAMPLE'] )
  else:
    print   ("#Existing output file ",config['REQUIRED']['FILEOUTSAMPLE'] )
  if not os.path.isfile ( os.getcwd()+"/"+config['REQUIRED']['FILEOUTBAD']) :
     print ("#Creating output file ",config['REQUIRED']['FILEOUTBAD'] )
     os.system('echo "'+ "\t".join(["ERRORTYPE"]+outputcolumns) 
              +'" > ' + config['REQUIRED']['FILEOUTBAD'] )
  else:
    print ("#Existing output file ",config['REQUIRED']['FILEOUTBAD'] )
  reverseheader ("LOAD BAGS", "digital")  
  for ofile in [config['REQUIRED']['FILEOUTSAMPLE'], config['REQUIRED']['FILEOUTBAD'] ] : 
    with open( ofile ) as outfile:
      ## load the table into memory for fast retrieval ###
      dictreader= csv.DictReader( outfile, delimiter ="\t")
      for row in dictreader:
         bag=row['OUTBINBAG']
         bagnum=int (row['BAGNUM'])
         samplenum=int (row['BAGSAMPLENUM'])
         if bagstatus[bag]['BAGNUM']<= bagnum:
           bagstatus[bag]['BAGNUM']= bagnum
           bagstatus[bag]['BAGSAMPLENUM']= samplenum
  
  print ("#CHECK  Current BAG NUM and SAMPLENUM")
  for key in bagstatus:
    print (key, bagstatus[key])
  
  ###################################################
  ##avoid messing with files if just validating 
  if args.validateonly:
    exit(0) 

  ################################################################
  ################################################################
  #mainloop scanning again and again and again until user aborts#
  ################################################################
    ## sophesticated https://pypi.org/project/pynput/
  print (termcolor.colored( pyfiglet.figlet_format( "---------\nSTARTING\n---------\n"
        ,font='banner'),'grey', attrs=['reverse'] ))
  while ( True ):
    ##SCAN 
    # banner, big
    print ( termcolor.colored( pyfiglet.figlet_format ("SCAN BARCODE", font='slant'),attrs=['bold'] ))
    rawinputcommand=input("SCAN BARCODE or !QUIT =>")
    fullbarcode=None
    if rawinputcommand == "!QUIT" :
      break 
    if rawinputcommand == '':
      print ("No input -- please rescan.")
      continue
    fout = "INPUT:"+rawinputcommand+" BARCODEMATCH:"+ config['REQUIRED']['BARCODEMATCH']
    #check for the barcode pattern 
    thematch=re.search(config['REQUIRED']['BARCODEMATCH'], rawinputcommand)
    if  not  thematch:
      ####### BADCODE -- no pattern
      print ("BAD: NO  PATTERN MATCHED IN BARCODE! [BARCODEMATCH:" , config['REQUIRED']['BARCODEMATCH'] )
      #failed match
      #SAVE OR SKIP
      response=get_validresponse( ["BADCODE", "SKIP/REDO"], 
            "Bag BADCODE ", "Scan 'BADCODE' as you bag OR 'SKIP/REDO' to not bag ->", 'magenta')
      if response=="BADCODE":
        write_outfilebad(rawinputcommand, 'FAILEDPATTERNMATCH', config, bagstatus )
      continue
    mbc = thematch.group(1)
    ### check if failed match in sample barcodes 
    if not mbc in matchedbarcodes:
      print ("BAD: BC PATTERN DOES NOT MATCH EXISTING BARCODE [BARCODEMATCH:"
          + config['REQUIRED']['BARCODEMATCH']+"]" )
      response=get_validresponse( ["BADCODE", "SKIP/REDO"], 
            "Bag BADCODE", "Scan 'BADCODE' as you bag OR 'SKIP/REDO' to not bag ->",'magenta')
      print (response)
      if response=="BADCODE":
        write_outfilebad(rawinputcommand, 'FAILEDSAMPLEMATCH', config, bagstatus )
      continue
    ### I have a good barcode match in the sample sheet but can have more than one
    if len(matchedbarcodes[mbc])==1:
      fullbarcode=matchedbarcodes[mbc][0]
    else:
      print ( "MULTIPLE MATCHES (", len(matchedbarcodes[mbc]), ") IN SAMPLES FROM [BARCODEMATCH:"
        + config['REQUIRED']['BARCODEMATCH']+"]"  ) 
      validinput=[]
      for pos, fbc in enumerate( matchedbarcodes[mbc] ):
         print ("  [",pos, "]  ", fbc)
         validinput.append(str(pos))
      response= get_validresponse (validinput , None, "Type # to choose or SKIP/REDO to abandon ->",'grey')
      if response== "SKIP/REDO":
        break
      if response in validinput:
        fullbarcode=matchedbarcodes[mbc][int(response)] 
    ### we should  have one and only one bar code
    print (fullbarcode)
    finalsample=samples[fullbarcodes[fullbarcode]]
    mappedbin=lookup_bin(config['BINS'], finalsample[header["BIN"]])
    response=get_validresponse( [mappedbin+mappedbin, "SKIP/REDO"], "Bag "+mappedbin,
          "Scan '"+mappedbin+"' as you bag OR  'SKIP/REDO to not bag ->",'green')
    if (response==mappedbin+mappedbin): 
      write_outfilesample(finalsample, mappedbin , rawinputcommand , config, bagstatus )
 
 ###########################################
 # functions for barcodescanbinner
 #########################################

def write_outfilebad(bc, errortype, config, bagstatus ):
#  outputcolumns=["OUTBINBAG","BAGNUM","BAGSAMPLENUM", "SCANNEDBARCODE", "BAGSETNAME"]dsadfasf
  bagstatus["BADCODE"]['BAGSAMPLENUM']+=1
  outtext=[  errortype
           ,"BADCODE"
           , str(bagstatus["BADCODE"]['BAGNUM'])
           , str(bagstatus["BADCODE"]['BAGSAMPLENUM'])
           , bc
           , config['REQUIRED']['BAGSETNAME']
          ]
  os.system('echo "'+ "\t".join(outtext) 
              +'" >> ' + config['REQUIRED']['FILEOUTBAD'] )
  bagcheck("BADCODE",bagstatus)
              
              
def write_outfilesample(samplerow, outbag , bc , config, bagstatus ):
  bagstatus[outbag]['BAGSAMPLENUM']+=1
  outtext= samplerow +[  outbag
           , str(bagstatus[outbag]['BAGNUM'])
           , str(bagstatus[outbag]['BAGSAMPLENUM'])
           , bc
           , config['REQUIRED']['BAGSETNAME']
          ]
  os.system('echo "'+ "\t".join(outtext) 
              +'" >> ' + config['REQUIRED']['FILEOUTSAMPLE'] )
  bagcheck (outbag, bagstatus)

def bagcheck ( bag, bagstatus):
  print (f"*** SAVED****  {bag} BAGNUMBER:", bagstatus[bag]['BAGNUM']
          , "SAMPLENUM:", bagstatus[bag]['BAGSAMPLENUM'] )
  if bagstatus[bag]['BAGSAMPLENUM'] >= bagstatus[bag]['BAGSIZE']:
     print (f"\nCurrent bag is labeled {bag} #"+str(bagstatus[bag]['BAGNUM']))
     print (f"PLEASE LABEL NEW BAG {bag} #"+str(bagstatus[bag]['BAGNUM']+1 ))
 
     response=get_validresponse(['BAG_CHANGED'], "#BAG SWAP# #NOW#", f"Bag {bag} " 
          + str(bagstatus[bag]['BAGNUM']) +" is full! SCAN BAG_CHANGED after swapping", 'blue',)
     bagstatus[bag]['BAGNUM']+=1
     bagstatus[bag]['BAGSAMPLENUM']=0
     print (" Thank you for swapping and properly storing old and labeling new bag !!!")
def get_validresponse(valid_responses, bigtext, smallquerytext,textcolor):
    while (True) :  
      if bigtext !=None:
        print ( termcolor.colored(pyfiglet.figlet_format (bigtext,font='standard') 
                  ,textcolor,attrs=['bold'])) 
      answer=input(smallquerytext )
      textcolor='red'
      if answer in valid_responses:
        return answer
      print (termcolor.colored("Not a valid response.",textcolor))
  
def lookup_bin (binlist, bincode ):
  #requires the config bin list and a bin code
  for key in binlist:
    textvalue=binlist[key]
    if textvalue.startswith("regex:"):
       pattern=textvalue.replace("regex:","")
       if re.search(pattern, bincode):
         return key
    else:
      #treat as regular string and test if strings match
      if (  textvalue==bincode):
        return key
  return "FAILEDBIN"
  
def reverseheader (text, font):
  print (termcolor.colored( pyfiglet.figlet_format(text,font=font),'grey', attrs=['reverse'] ))

 
#-----------------------------------------------------------------
#-----------------------------------------------------------------
#-----------------------------------------------------------------
#main is left to the bottom so that global variables can be populated

if __name__ == "__main__":
  print ("SYS VERSION", sys.version)
  print ("WORKING DIRECTORY", os.getcwd() )
  print ("COMMAND:"," " . join(sys.argv) )
  if len (sys.argv)==1:
    sys.argv.append("--help")  #if no command then it is a cry for help
  main(sys.argv[1:])

