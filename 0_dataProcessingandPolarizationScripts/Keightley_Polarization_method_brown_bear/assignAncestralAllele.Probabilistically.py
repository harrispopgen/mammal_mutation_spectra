#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This script will take output from extsfs that has been combined with extsfs input tables in separate bash script

 and assign ancestral alleles based on a Bernoulli trial with the probability of success (1) 
 
 as the probability of the major allele being ancestral from extsfs output. 

infile should have end with *pvalues.WITHPOSITIONS.usethis.txt.gz

# maybe combine them all into one file? 
Created on Fri May  6 10:44:06 2022

@author: annabelbeichman
"""
import gzip
import numpy as np
import argparse
 ####################### arguments #################
parser = argparse.ArgumentParser(description='This script will assign ancestral allele status based on prob of major allele being ancestral from ext-sfs')
parser.add_argument("--indir",required=True,help="path to dir containing extsfs output that has been combined with positions (files that end in *.pvalues.WITHPOSITIONS.usethis.txt.gz)",type=str)
parser.add_argument("--inprefix",required=True,help="prefix of infile from filename like extsfs.output.${inprefix}.pvalues.WITHPOSITIONS.usethis.txt.gz",type=str)
parser.add_argument("--intervalCount",required=True,help="total intervals -- varies depending on your species!! ",type=int)
parser.add_argument("--outdir",required=True,help="path to output dir (note this script doesn't make the dir, must exist)",type=str)
parser.add_argument("--outprefix",required=True,help="prefix for outout",type=str)


np.random.seed(42)  # setting seed so that should get same assignments if re-run though are random draws


args = parser.parse_args()
indir=str(args.indir)
inprefix=str(args.inprefix)
intervalCount=int(args.intervalCount)
outdir=str(args.outdir)
outprefix=str(args.outprefix)
# going to output 2 output files ; one with all the output and one with just 3 OOOOOOOOcolumns to be used downstream

# for testing: 
#outdir="/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/scripts/multispecies_spectra/polarization_assignAllelesProbabilistically/test/"
#outprefix=""
#outprefix="test"
#### these outputs will concat all the infiles
outfileFull_Filename=outdir+"/"+outprefix+".allIntervals.allInfoForTroubleshooting.DONOTUSEFORBCFTOOLSANNOTATE.txt"
outfileFull=open(outfileFull_Filename,"w")
outfile3Columns_Filename=outdir+"/"+"ThreeColumnsOnly."+outprefix+".allIntervals.USETHIS.txt" ### this will get used with bcftools annotate 
outfile3Columns=open(outfile3Columns_Filename,"w")
totalSitesCounter=0
sitesWhereAncDoesMatchRefCounter=0
# these counters are across all infiles

##### get the header from file 1 (same for all and only want to write to outfile once): ####
    
infile1_forHeader_filename=indir+"extsfs.output."+inprefix+"_00001"+".pvalues.WITHPOSITIONS.usethis.txt.gz"
infile1=gzip.open(infile1_forHeader_filename,"rt")

infile1.seek(0) # make sure it's at beginning

for line0 in infile1:
    line=line0.strip().split('\t') # this splits line by tabs
    if line[0]=="chr": # header
        infile1_header=line
        newheaderFull="\t".join(line)+"\tancAllele\n"
        outfileFull.write(newheaderFull)
        newheader3Columns="chr\tpos\tancAllele\n" # 3 columns only for use with bcftools annotate *** this is critical that this only be 3 columns in this order * 
        outfile3Columns.write(newheader3Columns)
        break
        
infile1.close() # close infile 1 (going to reread in next loop)
# for testing: infilename="/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/scripts/multispecies_spectra/polarization_assignAllelesProbabilistically/test/extsfs.output.KimuraModel.10MLRuns.focal.Mmd.out1.Mmm.out2.Ms_00001.pvalues.WITHPOSITIONS.usethis.txt.gz"

for interval in range(1,intervalCount+1): # aha you have to be careful with python range() because is not inclusive at end so must be +1; but don't want to start at 0 here because I don't have a 0 interval
    print("starting"+str(interval)) 
    intervalText = str(interval).zfill(5) # prepend zeroes so that total size is 5 (aka 00001 00002 00003 etc) 
    infilename=indir+"/extsfs.output."+inprefix+"_"+intervalText+".pvalues.WITHPOSITIONS.usethis.txt.gz" # format of infile (from my bash script)
    infile=gzip.open(infilename,"rt")

    # get the header: 
    infile.seek(0)
    for line0 in infile:
        line=line0.strip().split('\t') # this splits line by tabs
        if line[0]=="chr": # header --- to be safe, am gettig indices for every file to catch if one has different formatting for some reason
            header=line
            # make sure it's consistent with header of infile1
            if header!=infile1_header:
                print("this file"+intervalText+"has a different header than infile1 does - what is goign on!")
                exit(1)
            # get indices : (note that will pull first instance of string, works fine for this because nothing is repeated)
            # am doing this for every file even though redundant just as backup safety check
            chrHeaderIndex=header.index("chr")
            posHeaderIndex=header.index("pos")
            ProbMajorAlleleIsAncestralHeaderIndex=header.index("ProbMajorAlleleIsAncestral") # so this should work if there are 1 or 2 outgroups (not setting position manually)
            MajorAlleleHeaderIndex=header.index("major")
            MinorAlleleHeaderIndex=header.index("minor")
            # don't need idnexes of anything else 
            # optionally also going to get ref allele so that I can count up when anc doesn't match ref:
            RefAlelleHeaderIndex=header.index("ref")
            # **not ** writing out header each time because don't want header to appear over and over again in concatted file 
    
        else:
            # reset each time just in case:
            ancAllele=None # resetting each line to be cautious
            # now go through regular lines
            totalSitesCounter+=1
            # get chr and pos:
            mychr=line[chrHeaderIndex]
            mypos=line[posHeaderIndex]
            
            # get major and minor alleles
            mymajor=line[MajorAlleleHeaderIndex]
            myminor=line[MinorAlleleHeaderIndex]
            
            # optional : get ref alelle just for keeping track of what % of sites aren't ref polarized 
            myref=line[RefAlelleHeaderIndex] # only using ref to keep count of times when anc alelle doesn't match ref
            
            # get prob major allele is ancestral from extsfs
            ProbMajorAlleleIsAncestral=float(line[ProbMajorAlleleIsAncestralHeaderIndex]) # get probability
            
            # check that it's between 0 and 1
            if ProbMajorAlleleIsAncestral > 1.0 or ProbMajorAlleleIsAncestral < 0:
                print("Something is wrong with this probability!"+line0)
                exit(1)
            else:
                ## NOTE: doing a bernoulli trial where the probability of success is the probability that major alleel is ancestral
                # therefore success (outcome=1) means that anc allele is the major alelle
                # and failure  (outcome=0) means that the minor alelel is the anc allele (so not actually a failure)
                # https://numpy.org/doc/stable/reference/random/generated/numpy.random.binomial.html 
                # format (n,p) where n is how many draws you want to do . we are just drawing 1 site at a time so n = 1
                # p is probability of success aka prob major allele is ancestral
                bernoulliTrialResult=np.random.binomial(n=1,p=ProbMajorAlleleIsAncestral)
                
                # if trial result is 1 then major allele is ancestral
                # if trial result is 0 then minor allele is ancestral 
                if bernoulliTrialResult==1:
                    ancAllele=mymajor
                elif bernoulliTrialResult==0:
                    ancAllele=myminor
                else:
                    print("Something has gone wrong with bernoulli trial - result isn't 0 or 1")
                    print(line0)
                    exit(1)
                # write it out:
                # check if anc matches ref or not: (optional)
                if ancAllele!=myref:
                    sitesWhereAncDoesMatchRefCounter+=1 # add to the counter 
                # write it out:  
                outfileFull.write("\t".join(line)+"\t"+ancAllele+"\n") # write it out to full file
                outfile3Columns.write('\t'.join([mychr,mypos,ancAllele])+"\n")
    infile.close()
print("total sites: " + str(totalSitesCounter))
print("total sites where anc!=ref: "+str(sitesWhereAncDoesMatchRefCounter))
outfileFull.close()
outfile3Columns.close()
infile.close()
