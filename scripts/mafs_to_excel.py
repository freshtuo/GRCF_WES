#!/usr/bin/env python
# mafs_to_excel.py
# calculate target enrichment rate
# 

import os.path
import logging
import gzip
import pandas as pd
from re import search

# MAF files to process
maf_files = snakemake.input["maf"]

# output summary file
out_file = str(snakemake.output["report"])

# log file
log_file = str(snakemake.log)

# MAF columns to keep
cols = snakemake.config["cols"]

# rename columns if needed
newcols = snakemake.config["newcols"]

# MAF freeze panes
mrow = int(snakemake.config["mrow"])
mcol = int(snakemake.config["mcol"])

# functions
def setup_logging(log_file, level=logging.INFO):
    """set up logging"""
    # prepare loggings
    log_formatter = logging.Formatter('%(levelname)s: %(message)s')
    root_logger = logging.getLogger()
    root_logger.setLevel(level)
    
    # logging to file
    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setFormatter(log_formatter)
    root_logger.addHandler(file_handler)
    
    return root_logger

def load_maf(maf_file):
    """load and clean a MAF file"""
    logging.info('processing {}'.format(maf_file))
    # confirm header lines in MAF file
    with gzip.open(maf_file, 'rb') as fin:
        skip = 0
        for line in fin:
            if search('^#', line.decode()):
                skip += 1
            else:
                break
    logging.info('{} header lines to skip.'.format(skip))
    # load MAF file
    maf = pd.read_table(maf_file, header=0, sep='\t', skiprows=skip)
    logging.info('{} variants detected.'.format(maf.shape[0]))
    # clean MAF columns
    return maf

def combine_maf(maf_files, out_file, columns_to_keep, mrow, mcol):
    """write MAFs to a single excel file"""
    with pd.ExcelWriter(out_file, engine='xlsxwriter') as writer:
        for maf_file in maf_files:
            # get tumor-normal pair name: {pair}.funcotated.maf.gz
            pair = os.path.basename(maf_file).replace('.funcotated.maf.gz','')
            # load MAF
            maf = load_maf(maf_file)
            # get available columns
            columns_to_keep = check_columns(maf, columns_to_keep)
            # subset columns for output
            mafout = maf[columns_to_keep]
            # rename columns if needed
            if newcols:
                # get position of each available column in the original colume list
                pos = [cols.index(x) for x in columns_to_keep]
                # build a column name dictionary
                namedic = dict(zip([cols[k] for k in pos], [newcols[k] for k in pos]))
                # rename columns
                if namedic:
                    mafout = mafout.rename(columns=namedic)
            # write to file
            mafout.to_excel(writer, sheet_name=pair, index=False)
            # apply format
            worksheet = writer.sheets[pair]
            nrows, ncols = maf.shape
            for k in range(0, ncols):
                worksheet.set_column(k, k, 18)
            worksheet.freeze_panes(mrow,mcol)

def check_columns(maf, columns_to_keep):
        missing = []
        matching = []
        for x in columns_to_keep:
            if x not in maf.columns:
                missing.append(x)
            else:
                matching.append(x)
        if missing:
            logging.warning('cannot find columns: {}'.format(','.join(missing)))
        return matching

# main
# set up logging
root_logger = setup_logging(log_file, level=logging.INFO)
#root_logger = setup_logging(log_file, level=logging.DEBUG)

combine_maf(maf_files, out_file, cols, mrow, mcol)

logging.info('Complete!')

