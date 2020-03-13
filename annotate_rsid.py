#!/usr/bin/python
#-*- coding:utf-8 -*-
################################################
# File Name: annotate_rsid.py
# Author: Mischa Lundberg
# Mail: m.lundberg@uq.net.au
################################################

## TODO ##
##use dbSNP to update database

import os,sys
import pandas as pd
import argparse
import dask.dataframe as dd
import numpy as np
import subprocess
from dask.multiprocessing import get

input_data = ""
reference_data = ""

def get_ldsnp_db(db_type):

    directory = os.path.dirname
    print "Downloading and cleaning database"
    if db_type == "common" :
        print "Creating database for common variants"
        command = "cd "+directory+"/dbsnp ; rm *.*; \
            wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-common_all.vcf.gz;  zcat 00-common_all.vcf.gz | \
                grep -v \"##\" | cut  -f1-5 | \
                awk \\'{ if (\$1 == \"#CHROM\") print \"RSID\\tKEY1\\tKEY2\"; else print \$3\"\\t\"\$1\":\"\$2\"_\"\$4\"_\"\$5\"\\t\"\$1\":\"\$2\"_\"\$5\"_\"\$4}\\' \
                > dbsnp_common.tsv; rm *.gz; gzip dbsnp_common.tsv"
    if db_type == "all":
        print "Creating database for common and rare variants"
        command = "cd "+directory+"/dbsnp ; rm *.*; \
            wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz;  zcat 00-All.vcf.gz | \
                grep -v \"##\" | cut  -f1-5 | \
                awk \\'{ if (\$1 == \"#CHROM\") print \"RSID\\tKEY1\\tKEY2\"; else print \$3\"\\t\"\$1\":\"\$2\"_\"\$4\"_\"\$5\"\\t\"\$1\":\"\$2\"_\"\$5\"_\"\$4}\\' \
                > dbsnp_all.tsv; rm *.gz; gzip dbsnp_all.tsv"
    subprocess.call(command, shell=True)
    print "database successfully created"

def load_data(r_sep, i_sep, input_file, reference, strandedness, rck, file_type):

    #if rck: strandedness = rck
    supported_compression_formats = ['gz', 'gzip', 'bz2', 'zip', 'xz']
    dtypes = { 'CHR' : np.unicode_, 'BP' : np.int32, 'A1' : np.unicode_, 'A2' : np.unicode_, 'RSID' : np.unicode_, 'LOCUS' : np.unicode_, 'SE' : np.float, 'P' : np.float, 'FREQ' : np.float, 'BETA' : np.float}
    ################TODO#####################
    ## maybe add chunksize=500 to read_csv ##
    ## this will decrease the ram needed   ##
    #########################################
    
    ## reference file 
    if file_type == "reference":
        try:
            global reference_data
            reference_data = ""
            if rck:
                if os.path.splitext(reference)[-1][1:] in supported_compression_formats:
                    print "Compression detected for %s. Reading data." %reference
                    reference_data = pd.read_csv(reference, sep=r_sep, compression='infer', header=0, dtype=dtypes, usecols=['RSID', 'KEY'])
                else:
                    print "No compression detected for %s. Reading data." %reference
                    reference_data = pd.read_csv(reference, sep=r_sep, header=0, dtype=dtypes, usecols=['RSID', 'KEY'])
            else:    
                if os.path.splitext(reference)[-1][1:] in supported_compression_formats:
                    print "Compression detected for %s. Reading data." %reference
                    reference_data = pd.read_csv(reference, sep=r_sep, compression='infer', header=0, dtype=dtypes)
                else:
                    print "No compression detected for %s. Reading data." %reference
                    reference_data = pd.read_csv(reference, sep=r_sep, header=0, dtype=dtypes)
            
            reference_data.columns = reference_data.columns.str.upper()   
            ref_header = []
            for col in reference_data.columns: ref_header.extend([col])
            if not rck: reference_data.columns = check_and_swap_header(ref_header, "ref")
            check_strandedness(strandedness, reference, file_type, rck)
            #print reference_data.head()
            if not rck and not strandedness and ( not ('KEY1' in reference_data.columns) and not ('KEY1' in reference_data.columns)): 
                reference_data['KEY1'] = reference_data['CHR'].astype(str) + ":" + reference_data['BP'].astype(str) + "_" + reference_data['A1'].astype(str) + "_" + reference_data['A2']
                reference_data['KEY2'] = reference_data['CHR'].astype(str) + ":" + reference_data['BP'].astype(str) + "_" + reference_data['A2'].astype(str) + "_" + reference_data['A1']
            #print reference_data.head()
        except:
            raise Warning("An exception occurred while loading the reference file into memory.") 

    ## input file (sumstats)
    if file_type == "input":
        try:
            global input_data 
            input_data = "" 
            if os.path.splitext(input_file)[-1][1:] in supported_compression_formats:
                print "Compression detected for %s. Reading data." %input_file
                input_data = pd.read_csv(input_file, sep=i_sep, compression='infer', header=0, dtype=dtypes)
            else:
                print "No compression detected for %s. Reading data." %input_file
                input_data = pd.read_csv(input_file, sep=i_sep, header=0)
            input_data.columns = input_data.columns.str.upper()
            input_header = []
            for col in input_data.columns: input_header.extend([col])
            input_data.columns = check_and_swap_header(input_header, "input")
            check_strandedness(strandedness, input_file, file_type, rck)
            #print input_data.head()
            input_data['KEY2'] = input_data['CHR'].astype(str) + ":" + input_data['BP'].astype(str) + "_" + input_data['A1'].astype(str) + "_" + input_data['A2']
            input_data['KEY1'] = input_data['CHR'].astype(str) + ":" + input_data['BP'].astype(str) + "_" + input_data['A2'].astype(str) + "_" + input_data['A1']
            #print input_data.head()   
        except:
            raise Warning("An exception occurred while loading the input file %s into memory. Please check this input file and rerun it." %input_file ) 
 
    

def check_strandedness(strandedness, check_data, file_type, rck):

    global input_data   
    global reference_data

    if strandedness:
        if "_strandedness" not in check_data and file_type == "input":
            print "Adapt for strandedness for input(summarystats) file"
            input_data = adapt_strandedness (input_data)
            file_name = os.path.splitext(check_data)[0]+"_strandedness"+os.path.splitext(check_data)[1]
            input_data.to_csv(path_or_buf=file_name, sep=args.i_sep, index=False, compression='infer') #save_df(input_data, file_name, args.i_sep)
            print "Strandedness for input is finished. New input(summarystats) file is saved here: %s" %file_name
        else:
            print "Strandedness adaption detected in input(summarystats) file"

        if "_strandedness" not in check_data and not rck and file_type == "reference":
            print "Adapt for strandedness for reference file"
            reference_data = adapt_strandedness (reference_data)
            file_name = os.path.splitext(check_data)[0]+"_strandedness"+os.path.splitext(check_data)[1]
            reference_data['KEY'] = reference_data['CHR'].astype(str) + reference_data['BP'].astype(str) + reference_data['A1'].astype(str) + reference_data['A2']
            reference_data.to_csv(path_or_buf=file_name, sep=args.i_sep, index=False, compression='infer') #save_df(input_data, file_name, args.i_sep)
            key_reference_data = reference_data.filter(['RSID', 'KEY'], axis=1)
            file_name = os.path.splitext(check_data)[0]+"_strandedness_key"+os.path.splitext(check_data)[1]
            key_reference_data.to_csv(path_or_buf=file_name, sep=args.i_sep, index=False, compression='infer')
            print "Strandedness for reference is finished. New reference file is saved here: %s" %file_name
        else:
            print "Strandedness adaption detected in reference file"
    else:
        raise Warning("There is NO strandedness check information supplied.")
        #raise Warning("There is NO strandedness check implemented yet.\nTODO - check for same strand, if not on same strand change.")


def check_and_swap_header(df_header, data_type):

    needed_ref_header = ['CHR', 'BP', 'A1', 'A2', 'RSID']
    needed_input_header = ['CHR', 'BP', 'A1', 'A2']

    if data_type == "ref":
        needed_header = needed_ref_header
    elif data_type == "input":
        needed_header = needed_input_header
    else:
        sys.exit("Something went terribly wrong. Go to your corner and cry.")

    counter = 0
    header_lib = {'CHR' : 'CHR', 'BP' : 'BP', 'A1' : 'A1', 'A2': 'A2', 'RSID' : 'RSID', 'CHROMOSOME': 'CHR', 'CHROM': 'CHR', 'POS': 'BP', 'POSITION': 'BP', 'MINOR_ALLELE': 'A2', 'MAJOR_ALLELE': 'A1', 'RS': 'RSID'}
    #print df_header
    for item in range(len(df_header)):
        if df_header[item] in header_lib:
            df_header[item] = header_lib.get(df_header[item])
            if df_header[item] in needed_header:
                counter += 1
    
    if counter < len(needed_header):
            sys.exit("%s header doesnt contain needed information. Header should contain at leas: %s" %(data_type, needed_header))
    
    return df_header


def adapt_strandedness (df):

    df['A1_initial'] = df['A1']
    df['A2_initial'] = df['A2']
    df['A1'] = df['A1_initial'].str.replace('T','P')
    df['A1'] = df['A1'].str.replace('A','T/A')
    df['A1'] = df['A1'].str.replace('P','T/A')
    df['A1'] = df['A1'].str.replace('C','Q')
    df['A1'] = df['A1'].str.replace('G','G/C')
    df['A1'] = df['A1'].str.replace('Q','G/C')

    df['A2'] = df['A2_initial'].str.replace('T','P')
    df['A2'] = df['A2'].str.replace('A','T/A')
    df['A2'] = df['A2'].str.replace('P','T/A')
    df['A2'] = df['A2'].str.replace('C','Q')
    df['A2'] = df['A2'].str.replace('G','G/C')
    df['A2'] = df['A2'].str.replace('Q','G/C')
        
    return df

def strandedness_translation(allele_df):

    allele = allele_df.iloc[0]
    dna_lib = {'A': 'T/A', 'C': 'G/C', 'G': 'G/C', 'T': 'T/A'}
    new_allele = dna_lib.get(allele)
    allele_df = allele_df.replace(to_replace=allele, value=new_allele)
    return allele_df

def revert_strandedness_translation(df, rck):

    try:
        if rck: return df.drop(['A1','A2'], axis=1).rename(index=str, columns={"A1_initial":"A1", "A2_initial":"A2"})
            
        else: return df.drop(['A1_y','A2_y','A1_initial_y','A2_initial_y','CHR_y','BP_y','KEY','A1_x','A2_x'], axis=1).rename(index=str, columns={"CHR_x": "CHR", "BP_x": "BP", "A1_initial_x":"A1", "A2_initial_x":"A2"})
    except AssertionError as error:
        raise Warning("An exception occurred while reverting strandedness information. Error: %s" %error) 

def merge_dfs(strandedness):

    try:
        global input_data
        global reference_data
        if strandedness:
            return pd.merge(input_data, reference_data, on='KEY1')
        else:
            tempdf = pd.merge(input_data, reference_data, left_on='KEY1', right_on='KEY2')
            tempdf2 = pd.merge(input_data, reference_data,  left_on='KEY2', right_on='KEY1')
            tempdf3 = pd.concat([tempdf, tempdf2], axis=1)
            tempdf = pd.merge(input_data, reference_data, on='KEY1')
            tempdf = pd.merge(input_data, reference_data, on='KEY2')
            return pd.concat([tempdf, tempdf2, tempdf3], axis=1)
    except AssertionError as error:
        raise Warning("An exception occurred while merging dataframes. Error: %s" %error) 

def save_df(df, output_loc, o_sep):

    df.to_csv(sep=o_sep, header=True, path_or_buf=output_loc, index=False)


def main(args):

    if args.sort: raise Warning("Not yet implemented! TODO")
    if args.i == 'none' or args.create_database:
        if not args.create_database:
            raise Exception("No input file has been selected.")
        else:
            raise Warning("No input file selected. Updating dbsnp database!")
            get_ldsnp_db(args.snps)
    else:
        load_data(args.r_sep,args.i_sep,args.i,args.r,args.strandedness, args.rck, "reference")
        if "," in args.i or "files.txt" in args.i or args.m: ## if several input files (sumstats) given.
            print "in files"
            if "files.txt" in args.i or args.m:
                if args.mis == 't':
                    with open(args.i, 'r') as f: input_file_list = f.read().splitlines()
                else:
                    with open(args.i, 'r') as f: input_file_list = f.readlines()[0].split(args.mis)[:-1]
                #print "files.txt"
            else:
                input_file_list = args.i.split(",")
                #print "list"
            print input_file_list
            o_was_none =  False
            for input_file in input_file_list:
                print "Running now %s" %input_file
                load_data(args.r_sep,args.i_sep,input_file,args.r,args.strandedness, args.rck, "input")
                merged_df = merge_dfs(args.strandedness)
                print merged_df.head()
                if args.strandedness:
                    merged_df = revert_strandedness_translation(merged_df, args.rck)     
                try:  
                    if args.o == "none":
                        args.o = os.path.splitext(input_file)[0]+"_sumstats_rsid.bed"
                        o_was_none = True
                    print "output is saved in %s" %args.o
                except AssertionError as error:
                    raise Warning("An exception occurred while generating output filename. Error: %s" %error)
                try:
                    save_df(merged_df, args.o, args.o_sep)
                except AssertionError as error:
                    raise Warning("An exception occurred while saving output dataframe. Error: %s" %error)
                if o_was_none:
                    args.o = "none"
                
        else: 
            load_data(args.r_sep,args.i_sep,args.i,args.r,args.strandedness, args.rck, "input")
            merged_df = merge_dfs(args.strandedness)
            print merged_df.head()
            if args.strandedness:
                merged_df = revert_strandedness_translation(merged_df, args.rck)
            if args.o == "none":
                args.o = os.path.splitext(args.i)[0]+"_sumstats_rsid.bed"
            save_df(merged_df, args.o, args.o_sep)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Annotates RSID to chromosmal positions.\nPlease remember to check if your reference and input file are both either onebased or zerobased.')
    parser.add_argument('-i', default='none', help='Target input summary statistics file or list of files separated by \",\" or a file containing all input files separated by \",\" (the file name of this has to be files.txt, e.g. generated by the following command:  ls *.fastGWA.txt.gz | tr \'\n\' \',\' > files.txt); has to contain CHR, BP, A1, A2 (these labels are needed for annotation)')
    parser.add_argument('-r', default="./dbsnp/dbsnp.tsv.gz", help='Reference file containing RSIDs and chromosomal positions with allele information; has to contain CHR, BP, RSID, A1, A2 (these labels are needed for annotation). Deflaut: \"./dbsnp/dbsnp.tsv.gz\"')
    parser.add_argument('-o', default='none', help='Target output File, for name of the file in argument i leave blank (will be extended by \"sumstats_rsid.bed\"')
    parser.add_argument('--nochr', default=False, action='store_true', help='Make both input files nochr (no \"chr\" prefix to the chromosomal location)')
    parser.add_argument('--sort', default=False, action='store_true', help='Sorts the BED File')
    parser.add_argument('--i_sep', default="\t", help='Set the separator used for the -i option (\\t is default)')
    parser.add_argument('--r_sep', default="\t", help='Set the separator used for the -r option (\\t is default)')
    parser.add_argument('--o_sep', default="\t", help='Set the separator used for the -o option (\\t is default)')
    parser.add_argument('--strandedness', default=False, action='store_true', help='Set if strandedness is unknown or different between the two input files')
    parser.add_argument('-rck', default=False, action='store_true', help='(reference contains key) - If reference file contains a column \"KEY\" which harbors (), then set this key. In this case, the needed KEY is not generated, but used from your input. This will speed-up the whole process by far (depending on the size of your reference file). Setting this key, activates the strandedness option!')
    parser.add_argument('-mis', default=',', help='mass input separator; separator of the files to use if using files.txt as input argument for -i option. Default is \',\'')
    parser.add_argument('-m', default=False, action='store_true', help='Set this key, if autodetection of a file containing filenames to run as input doesnt work.')
    parser.add_argument('--create_database', default=False, action='store_true', help='downloads and prepares the snpdb database that is used for annotating SNPs without rsid. Needs to be ran as first initialization step')
    parser.add_argument('--snps', default="common", help='select which database you want to use/align for; options: \"common\" (only common SNPs) \"all\" (common and rare SNPs);deflaut \"common\"')
    args = parser.parse_args()
    main(args)
#import os,sys
#import pandas as pd
#import argparse
#import dask.dataframe as dd
#import numpy as np
#from dask.multiprocessing import get
#reference_data = pd.read_csv('./snp151CodingDbSnp_onebased_nochr_header_short_stradedness.bed', sep='\t', compression='infer', header=0)
#input_data = pd.read_csv('./25690.v1.fastGWA.txt_stradedness.gz', sep='\t', compression='infer', header=0)
#reference_data['KEY'] = reference_data['CHR'].astype(str) + reference_data['BP'].astype(str) + reference_data['A1'].astype(str) + reference_data['A2']
#input_data['KEY'] = input_data['CHR'].astype(str) + input_data['BP'].astype(str) + input_data['A1'].astype(str) + input_data['A2']
#temp = pd.merge(input_data, reference_data, on='KEY')
#temp.drop(['A1_y','A2_y','A1_initial_y','A2_initial_y','CHR_y','BP_y','KEY','A1_x','A2_x'], axis=1).rename(index=str, columns={"CHR_x": "CHR", "BP_x": "BP", "A1_initial_x":"A1", "A2_initial_x":"A2"})
#echo -e "CHR\tPOS\tRSID\tA1\tA2" > snp151_onebased.bed ;awk '{print $2"\t"$4"\t"$5"\t"$10}' snp151.txt | sed 's/\//\t/g' | grep -v lengthTooLong >> snp151_onebased.bed
#sed 's/chr//g' snp151_onebased.bed > snp151_onebased_nochr.bed
#echo -e "CHR\tPOS\tRSID\tA1\tA2" > snp151_zerobased.bed ;awk '{print $2"\t"$3"\t"$5"\t"$10}' snp151.txt | sed 's/\//\t/g' | grep -v lengthTooLong >> snp151_zerobased.bed
#sed 's/chr//g' snp151_zerobased.bed > snp151_zerobased_nochr.bed