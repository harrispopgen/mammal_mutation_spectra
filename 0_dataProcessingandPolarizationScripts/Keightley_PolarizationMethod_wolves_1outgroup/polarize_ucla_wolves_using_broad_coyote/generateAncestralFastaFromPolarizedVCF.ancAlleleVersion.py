#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 09:32:10 2020

@author: annabelbeichman

goal: make 'ancestral' fasta based on AA calls in a vcf file's INFO field

need to make a copy of reference fasta *do NOT modify original in place*

and have a vcf file with ancestral alleles in the AA= field

and then want to update the fasta with the coordinates from the vcf

*be careful of coordinates. cyvcf2 and pyfaidx are consistent with each other.*

If a site was lower case (soft masked) in the original fasta, the replacement will be soft-masked too.
!! note this script assumes the vcf and fasta are both on positive strand. no reverse complementing will occur !!

 # TODO: decide on behavior if no AA field is present. Could just break since that is not the purpose of the script (every site should have AA even if its AA=N)
 # or could sub in N or the ref allele. Ask Will for preference. For now it just exits.
## NOTE: currently am inputting a vcf that has ambiguous anc state sites REMOVED # so they don't get updated to N or anything. So the anc fasta will just have ref fasta bp at those sites. NOT IDEAL need to be cautious (warned Tom)
"""
import os
import sys
import cyvcf2
import pyfaidx
import shutil # for file copying
import argparse
######### set up args #############
parser = argparse.ArgumentParser(description='This script will make a copy of your reference genome fasta and sub in ancestral allele calls from the vcf file (in the ancAllele='' info field).')
parser.add_argument("--vcffile",required=True,help="path to input vcf that has had ancestral alleles assigned in the AA field. Sites where no AA was called should have ancAllele=N rather than no AA info field. BUT NOTE THAT IF A SITE WAS REMOVED FROM VCF FOR BEING AMBIGUOUS IT WONT GET UPDATED BY THIS SCRIPT")
parser.add_argument("--originalReferenceFasta",required=True,help="path to original reference fasta file. This file will be copied and the copy will be modified.")
parser.add_argument("--refPrefix",required=True,help="reference fasta prefix for modified fasta file, e.g. brown_bear")
parser.add_argument("--outdir",required=True,help="path to where you want modified fasta to be written out (for safety's sake, a different directory from the input fasta is recommended)")

args = parser.parse_args()
vcffilename=str(args.vcffile)
originalfastafilename= str(args.originalReferenceFasta) #DON'T CHANGE THE ORIGINAL FASTA!
refPrefix=str(args.refPrefix)
outdir=str(args.outdir)
mutablefastafilename=outdir+"/MODIFIED.ANCESTRAL."+refPrefix+".fasta" ### THIS IS THE FASTA YOU CAN CHANGE.


# for testing purposes
#outdir="/Users/annabelbeichman/Documents/UW/BearProject/scripts/sandbox/makingAncestralFasta/"
#vcffilename="/Users/annabelbeichman/Documents/UW/BearProject/scripts/sandbox/testAA.out.vcf.gz"
##must have .fasta file
#originalfastafilename="/Users/annabelbeichman/Documents/UW/BearProject/DataInformation/reference_genomes/brown_bear/brown_bear.fasta" # THIS IS THE ORIGINAL FASTA~ be very careful and do not mess with this file!!
#refPrefix="brown_bear"
#mutablefastafilename=outdir+"MODIFIED.ANCESTRAL."+refPrefix+".fasta" ### THIS IS THE FASTA YOU CAN CHANGE.


#%% functions
def make_ancestral_fasta_fromVCFAA(vcffilename,mutablefastaname):
    '''


    Parameters
    ----------
    vcffilename : str
        Path to input vcf file.
    mutablefastaname : str
        path to mutable fasta file that is NOT the original fasta because will get modified in place.

    Returns
    -------
    None.

    '''

    #### open the mutable fasta copy. this also will generate a new .fai file for the mutable fasta file
    print("opening fasta (mutable)")
    mutableFasta = pyfaidx.Fasta(mutablefastaname, mutable=True) ### THIS FILE WILL BE MODIFIED IN PLACE. Key function will split the fasta header
    # note if you have a weird header you can use key_function = lambda x: x.split(':')[0]
    # >NW_020656121.1:0-92727749 at the ":" and only use the first part since that's the scaffold name. From the pyfaidx docs.
    # read through the vcf file and update the fasta file. (# note the vcf file will not be modified.)
    for variant in cyvcf2.VCF(vcffilename):
        vcfAA=None
        vcfref = variant.REF
        vcfchrom = variant.CHROM
        vcfsitestart = variant.start # cyvcf converts to 0 based. So if vcf coord is 1492 then this will be 1491 and end will be 1492 (like a bed file)
        vcfsiteend = variant.end
        if not variant.INFO.get('ancAllele'):
            print("no AA")
            print("A site has no AA info field. Did you polarize this vcf? ...") # if a site has no AA field, for now, break the script. but could instead do either of the following
            os.remove(mutablefastafilename)
            sys.exit(1)
            # TODO: choose an alternative to exiting when you run into a site without an AA field.
            # alternative 1, set to N
            #vcfAA='N'
            # or alternative 2, set to ref
            #vcfAA=vcfref
        # make sure the scaffold is in the fasta file
        elif not vcfchrom in mutableFasta.keys():
            print("scaffold not in fasta")
            print("VCF scaffold/chromosome not found in fasta file, exiting ...")
            os.remove(mutablefastafilename)
            sys.exit(1)
        else:
            vcfAA =  variant.INFO.get('ancAllele') # get from AA = info field (must be pre-polarized.)
            # the coordinates of pyfaidx are start+1 and end =end (like bed file) so if you do [200:230] you'll get back 201:230. The manual states this in an unclear manner
            # but I have checked it. Also there is an internal check in this script that makes sure the reference base matches between vcf and fasta.
            fastaReferenceBase=str(mutableFasta[vcfchrom][vcfsitestart:vcfsiteend]) # convert to string and capitalize. these coordinates work.
            # if the fasta ref base is lower case, that means it's in a soft-masked repeat region
            # and I want the replacement to stay lower case in case I decide to use that softmasking in mutyper
            # note that fasta base may be lower case. so capitalize it for comparison to the vcfRef base
            if fastaReferenceBase.capitalize() != vcfref:
                print("something is wrong. your fasta position doesn't match the vcf's ref allele.\nAre your positions wrong?\n0 vs 1 based issues?\nDid you somehow alter the ref fasta?")
                os.remove(mutablefastafilename)
                sys.exit(1)
            else:
                # check if fastaRefBase is upper or lower case:
                # update fasta with AA in place. **Note** the fasta is changed in place so be very careful with this!
                # if lower case
                if fastaReferenceBase.islower():
                    mutableFasta[vcfchrom][vcfsitestart:vcfsiteend] = vcfAA.lower() # # make replacement base lower case if it was lower case in the original fasta.
                # if upper case
                else:
                    mutableFasta[vcfchrom][vcfsitestart:vcfsiteend] = vcfAA

                # note if you have AA = N it will replace with an 'N'.
    print("closing fasta")
    mutableFasta.close()
    return()


def checkAncestralFasta(vcffilename,fastaToCheck):
    '''
    Helper function that you can use to check if a fasta file has already been udpated with AA alleles
    Note this doesn't check for correct capitalization.

    Parameters
    ----------
    vcffilename : str
         path to vcf file that contains the AA info fields.
    fastaToCheck : str
        path of fasta you want to check, ie one that has had alleles changed, to be checked to make sure that previous function completed.

    Returns
    -------
    None.

    '''
    # open up fasta for checking NOT in mutable state (unchangeable)
    print("opening fasta (unmutable)")
    fastaToCheck=pyfaidx.Fasta(fastaToCheck, mutable=False) # this is the fasta you just updated, but you're reopening it

    for variant in cyvcf2.VCF(vcffilename):
        vcfchrom = variant.CHROM
        vcfsitestart = variant.start # cyvcf converts to 0 based. So if vcf coord is 1492 then this will be 1491 and end will be 1492 (like a bed file)
        vcfsiteend = variant.end
        if not variant.INFO.get('ancAllele'):
        	print("A site has no AA info field. Did you polarize this vcf?")
        	sys.exit(1) # if a site has no AA field, for now, break the script. but could instead do either of the following
        elif not vcfchrom in fastaToCheck.keys():
            print("VCF scaffold/chromosome not found in fasta file")
            sys.exit(1)
        else:
            vcfAA =  variant.INFO.get('ancAllele') # get from AA = info field (must be pre-polarized.)
            # the coordinates of pyfaidx are start+1 and end =end (like bed file) so if you do [200:230] you'll get back 201:230. The manual states this in an unclear manner
            # but I have checked it. Also there is an internal check in this script that makes sure the reference base matches between vcf and fasta.
            fastaAncestralBase=str(fastaToCheck[vcfchrom][vcfsitestart:vcfsiteend]).capitalize() # convert to string and capitalize. these coordinates work.
            # if the fasta base already matches vcfAA, you are fine.
            if fastaAncestralBase==vcfAA:
                continue
            else:
                #print("fasta base = "+str(fastaAncestralBase)+"\nVCF AA = "+vcfAA)
                print("FAIL"+str(vcfchrom)+":"+str(vcfsitestart)+"-"+str(vcfsiteend)+"Checking failed. Fasta ancestral base doesn't match VCF AA base.")
                sys.exit(1)

    print("closing fasta")
    fastaToCheck.close()
    print("The ancestral fasta matches the vcf AA calls. This ancestral fasta has been made correctly.")
    return()




#%% main
def mainFunction(vcffilename,originalfastafilename,mutablefastafilename):
    #### make a copy of the fasta file labeled MODIFIED and only modify that!
    shutil.copy(originalfastafilename,mutablefastafilename)
    make_ancestral_fasta_fromVCFAA(vcffilename,mutablefastafilename)
    #### check to make sure that you actually updated all the sites (this isn't necessary to do every time but is useful for script testing)
    # if you aren't sure if the fasta file has been updated correctly, you could just run checkAncestralFasta to make sure the vcf and fasta are consistent
    checkAncestralFasta(vcffilename,mutablefastafilename)
    # rename it if script makes it this far (a way to make sure that job didn't end prematurely)
    os.rename(mutablefastafilename,outdir+"/FINAL.MODIFIED.ANCESTRAL."+refPrefix+".fasta")
    os.rename(mutablefastafilename+".fai",outdir+"/FINAL.MODIFIED.ANCESTRAL."+refPrefix+".fasta.fai")
#%% run it
mainFunction(vcffilename,originalfastafilename, mutablefastafilename)

# TODO institute a check for an excess of lowerUPPERlower sites in the new ancestral fasta indicating the lower/upper chagnes didn't work
