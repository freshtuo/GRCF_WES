#!/usr/bin/env python
# review_base_coverage.py
# summarize per-base coverage results generated by samtools
# 

import os.path
import logging
import gzip
import pandas as pd
import pybedtools as bd
from re import search

# base coverage file
cov_file = snakemake.input["cov"]

# output summary file
out_file = str(snakemake.output["report"])

# temperary bed file
temp_file = str(snakemake.output["tempbed"])

# target bed file (including gene annotations)
bed_file = str(snakemake.config["bed"])

# low coverage regions file
lowcov_file = str(snakemake.output["lowcov"])

# log file
log_file = str(snakemake.log)

# sample ids
samples = list(snakemake.config["read1"].keys())

# cutoff for defining a low-coverage base
lowcut = int(snakemake.config["lowcut"])

# base coverage: %>?X
cstart = int(snakemake.config["covstart"])
cend = int(snakemake.config["covend"])
cstep = int(snakemake.config["covstep"])

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

def summarize_coverage(cov_file, out_file, temp_file, lowcov_file, lowcut=10):
    """load and summarize the per-base coverage results"""
    logging.info('processing {}'.format(cov_file))
    # read in per-base coverage results
    cov = pd.read_table(cov_file, header=None, names=['chrom','pos']+samples, sep='\t', low_memory=False)
    logging.info('{} bases and {} samples detected.'.format(cov.shape[0], cov.shape[1]))
    # dictionary to store summary results
    # for creating DataFrame
    results = { 'sample':samples }
    # review coverage
    # minimum coverage per sample
    results.update({'mincov': [cov[sid].min() for sid in results['sample']]})
    # maximum coverage per sample
    results.update({'maxcov': [cov[sid].max() for sid in results['sample']]})
    # mean coverage per sample
    results.update({'meancov': [cov[sid].mean().round(2) for sid in results['sample']]})
    # median coverage per sample
    results.update({'mediancov': [cov[sid].median() for sid in results['sample']]})
    # fold 80 base penalty
    results.update({'fold80': [cov[sid].mean()/cov[sid].quantile(0.2) for sid in results['sample']]})
    # % bases > ?X
    for c in range(cstart, cend+cstep, cstep):
        results.update({'%>{}X'.format(c): [(cov[sid] > c).sum()/len(cov) for sid in results['sample']]})
    # save coverage results
    with pd.ExcelWriter(out_file, engine='xlsxwriter') as writer:
        # write data
        pd.DataFrame(results).to_excel(writer, sheet_name='coverage', index=False)
        # get book/sheets object
        workbook = writer.book
        worksheet = writer.sheets['coverage']
        # set percentage format
        cell_format1 = workbook.add_format()
        cell_format1.set_num_format(10)# 0.00%
        # set two digits number format
        cell_format2 = workbook.add_format()
        cell_format2.set_num_format('0.00')
        # apply format
        worksheet.set_column(0, 0, 18)
        worksheet.set_column(3, 3, 16, cell_format2)
        worksheet.set_column(4, 4, 16)
        worksheet.set_column(5, 5, 12, cell_format2)
        for k in range(6, int((cend-cstart)/cstep)+7):
            worksheet.set_column(k, k, 12, cell_format1)
        worksheet.freeze_panes(1,0)
    # check low coverage regions
    detect_low_cov_regions(cov, temp_file, lowcov_file, lowcut)

def detect_low_cov_regions(cov, temp_file, lowcov_file, lowcut=10):
    """search for regions of low-coverage"""
    # select coverage columns
    columns_to_check = [x for x in cov.columns.tolist() if x not in ['chrom', 'pos']]
    # find bases of consistently low coverage
    cov_df = cov[(cov[columns_to_check] < lowcut).all(axis=1)][['chrom','pos']]
    # prepare a bed file
    cov_df['start'] = cov_df['pos'] - 1
    cov_df['end'] = cov_df['pos']
    cov_df.to_csv(temp_file, columns=['chrom','start','end'], sep='\t', index=False, header=False)
    # load bed file
    lowbed = bd.BedTool(temp_file)
    logging.info('{} low coverage bases.'.format(len(lowbed)))
    # sort and merge
    lowbed = lowbed.sort().merge()
    logging.info('{} regions after merging bases.'.format(len(lowbed)))
    # load target capture bed file (including gene annotations)
    refbed = bd.BedTool(bed_file)
    # annotate low coverage regions
    annotated_lowbed = lowbed.intersect(refbed, wa=True, wb=True).to_dataframe()
    # write to file
    annotated_lowbed.to_csv(lowcov_file, sep='\t', index=False)

# main
# set up logging
root_logger = setup_logging(log_file, level=logging.INFO)
#root_logger = setup_logging(log_file, level=logging.DEBUG)

summarize_coverage(cov_file, out_file, temp_file, lowcov_file, lowcut)

logging.info('Complete!')

