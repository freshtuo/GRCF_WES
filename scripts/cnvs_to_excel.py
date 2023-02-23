#!/usr/bin/env python
# cnvs_to_excel.py
# combine annotated CNVs and write to an excel file
# 

import os.path
import logging
import gzip
import pandas as pd
from re import search

# CNV files to process
cnv_files = snakemake.input["cnv"]

# output summary file
out_file = str(snakemake.output["report"])

# log file
log_file = str(snakemake.log)

# filter regions without cnvs?
filt = bool(snakemake.params["filt"])

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

def load_cnv(cnv_file):
    """load and clean a CNV file"""
    logging.info('processing {}'.format(cnv_file))
    # load CNV file
    cnv = pd.read_table(cnv_file, header=0, sep='\t', low_memory=False)
    logging.info('{} CNV entries detected.'.format(cnv.shape[0]))
    # filter regions without cnvs
    if filt:
        cnv = cnv[cnv['cn'] != 2]
        logging.info('filter out entries without cnvs, {} entries left.'.format(cnv.shape[0]))
    return cnv

def combine_cnv(cnv_files, out_file):
    """write CNVs to a single excel file"""
    with pd.ExcelWriter(out_file, engine='xlsxwriter') as writer:
        for cnv_file in cnv_files:
            # get tumor-normal pair name: {pair}.call.annotated.txt.gz
            pair = os.path.basename(cnv_file).replace('.call.annotated.txt.gz','')
            # load CNV
            cnv = load_cnv(cnv_file)
            # write to file
            cnv.to_excel(writer, sheet_name=pair, index=False)
            # apply format
            worksheet = writer.sheets[pair]
            nrows, ncols = cnv.shape
            for k in range(0, ncols):
                worksheet.set_column(k, k, 18)
            ##worksheet.freeze_panes(1,1)

# main
# set up logging
root_logger = setup_logging(log_file, level=logging.INFO)
#root_logger = setup_logging(log_file, level=logging.DEBUG)

combine_cnv(cnv_files, out_file)

logging.info('Complete!')

