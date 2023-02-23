#!/usr/bin/env python
# annotate_cnv.py
# annotate copy number calls based on annotations from the target bed file
# 

import logging
import pandas as pd
import pybedtools as bd
from re import search

# copy number variations (cnv) file
cnv_file = str(snakemake.input["call"])

# output annotated cnv file
out_file = str(snakemake.output["call"])

# temporary bed file
temp_file = str(snakemake.output["tempbed"])

# target bed file (including gene annotations)
bed_file = str(snakemake.config["bed"])

# log file
log_file = str(snakemake.log)

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

def get_gene_name(tx):
    """given an annotation, extract gene name"""
    tglist = []
    for y in tx.split(','):
        tpat = search('entg\|(.*)', y)
        if tpat:
            tg = tpat.groups()[0]
            if tg not in tglist:
                tglist.append(tg)
    return ';'.join(tglist)

def collapse_genes(td):
    """merge gene annotations that target the same region"""
    tglist = []
    for y in td['gene'].tolist():
        for z in y.split(';'):
            if z != '' and z not in tglist:
                tglist.append(z)
    if tglist:
        return ';'.join(tglist)
    else:
        return '-'
    

# annotate cnv calls
def annotate_cnvs(cnv_file, out_file, temp_file, bed_file):
    """prepare a bed file and annotate with bedtools"""
    # load cnv calls
    cnv = pd.read_table(cnv_file, header=0, sep='\t', low_memory=False)
    logging.info('CNV calls: {}'.format(cnv.shape))
    # prepare a bed file
    cnv.to_csv(temp_file, sep='\t', header=False, index=False, columns=['chromosome','start','end'])
    # load bed file
    bed = bd.BedTool(temp_file)
    logging.info('{} cnv regions.'.format(len(bed)))
    # load target capture bed file (including gene annotations)
    refbed = bd.BedTool(bed_file)
    # annotate cnv regions
    annotated_bed = bed.intersect(refbed, wa=True, wb=True).to_dataframe(low_memory=False)
    # clean gene annotations
    annotated_bed['gene'] = annotated_bed['thickStart'].apply(get_gene_name)
    # merge annotations on the same region
    collapsed_bed = annotated_bed.groupby(by=['chrom','start','end']).apply(collapse_genes).rename(index='gene').reset_index().rename(columns={'chrom':'chromosome'})
    # combine annotations with the original cnv results, and write to file
    cnv.drop(labels='gene', axis=1).merge(collapsed_bed, how='left', on=['chromosome','start','end']).to_csv(out_file, sep='\t', index=False)

# main
# set up logging
root_logger = setup_logging(log_file, level=logging.INFO)
#root_logger = setup_logging(log_file, level=logging.DEBUG)

annotate_cnvs(cnv_file, out_file, temp_file, bed_file)

logging.info('Complete!')

