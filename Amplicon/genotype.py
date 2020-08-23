'''

dab-seq: single-cell dna genotyping and antibody sequencing
ben demaree 7.9.2019

the main script for pipeline execution

'''

import os
import subprocess
import argparse
import copy
import sys
import time
from multiprocessing.pool import ThreadPool
from multiprocessing import Process, Queue

# import functions from external files
import resources


def wait(processes):
    # waits for processes to finish
    return [process.communicate() for process in processes]


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='''
    
    dab-seq: single-cell dna genotyping and antibody sequencing
    ben demaree 2019
    
    input requirements:
    -config file defining file paths and variables (dabseq.cfg)
    -raw fastq files (panel and optionally antibody tags)
    -cell and ab barcode csvs
    -panel bed files
    
    requires the following programs in path:
    -gatk
    -bowtie2
    -itdseek (flt3-calling only)
    -samtools
    -bedtools
    -bcftools
    -cutadapt
    -bbmap
    -snpeff
    
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('sample_name', type=str, help='sample basename')
    parser.add_argument('cfg_file', type=str, help='config filename')
    parser.add_argument('--skip-flt3', action='store_true', default=False, help='option to skip FLT3-ITD calling')

    args = parser.parse_args()  # parse arguments

    sample_basename = args.sample_name
    cfg_f = args.cfg_file
    dna_only = True
    skip_flt3 = args.skip_flt3

    print '''Initializing pipeline...    '''

    # load config file variables
    # be careful about using exec
    if not os.path.isfile(cfg_f):
        print 'Config file not found! Please check the file name and path.'
        raise SystemExit

    else:
        with open(cfg_f, 'r') as cfg:
            for line in cfg:
                if line[0] == '#' or line[0] == '[' or line[0] == ' ':
                    continue
                else:
                    var = line.split("#", 1)[0].strip()  # to remove inline comments
                    exec(var)

    # check all files exist
    all_vars = copy.copy(globals())
    input_files = [all_vars[f] for f in all_vars if '_file' in f and f != '__file__']
    input_files.append(bt2_ref + '.1.bt2')
    missing_files = []
    for f in input_files:
        if not os.path.exists(f):
            missing_files.append(f)

    # print missing files, if any, and exit
    if missing_files != []:
        print 'The following input files could not be found:'
        for f in missing_files:
            print f
        print '\nExiting...\n'
        raise SystemExit

    # check that the input fastq directories exist and are not empty
    if not os.path.exists(panel_fastq_R1):
        print 'FASTQ input directory (%s) does not exist! Exiting...\n' % panel_fastq_R1
        raise SystemExit

    # create all other directories for this run
    # if it already exists, ignore and continue
    dirs = [all_vars[d] for d in all_vars if '_dir' in d]
    dirs.sort()
    for d in dirs:
        if not os.path.exists(d):
            os.mkdir(d)

    print '''

###################################################################################
Step 4: count read alignments to inserts
###################################################################################
'''

    # get r1 reads for all panel samples
    panel_r1_files = [panel_fastq_R1]

    # align r1 reads to inserts to obtain read counts across all barcodes
    all_tsv = barcode_dir + sample_basename + '.all.tsv'  # alignment counts for all barcodes
    resources.count_alignments(panel_r1_files, amplicon_file, human_fasta_file, all_tsv, temp_dir)
    #raise SystemExit
    print '''
###################################################################################
Step 5: call valid cells using selected method
###################################################################################
'''

    # call valid cells using cell_caller function
    valid_cells = os.listdir(by_cell_fastq_dir)

    # create SingleCell objects for each valid cell
    cells = [resources.SingleCell(barcode,
                                  by_cell_fastq_dir,
                                  by_cell_bam_dir,
                                  by_cell_gvcf_dir,
                                  by_cell_flt3_dir)
             for barcode in valid_cells]

    print '%s valid cells found!\n' % len(cells)


    print '''
####################################################################################
# Step 7: align, convert, sort, index panel reads (optional: call FLT3-ITDs)
####################################################################################
'''

    # limit number of cells to preprocess at a time (based on hardware limitations)
    n_preprocess = 300

    # create pool of workers and run through all samples
    preprocess_pool = ThreadPool(processes=n_preprocess)

    # align and index cells
    for c in cells:
        preprocess_pool.apply_async(resources.SingleCell.align_and_index, args=(c, bt2_ref,))

    preprocess_pool.close()
    preprocess_pool.join()

    # optionally, call FLT3-ITDs using ITDseek
    if not skip_flt3:

        preprocess_pool = ThreadPool(processes=n_preprocess)

        for c in cells:
            preprocess_pool.apply_async(resources.SingleCell.call_flt3, args=(c, human_fasta_file,))

        preprocess_pool.close()
        preprocess_pool.join()

    else:
        os.rmdir(by_cell_flt3_dir)

    print '''
####################################################################################
# Step 8: perform variant calling for all cells
####################################################################################
'''

    # limit number of cells to call variants at a time (based on hardware limitations)
    n_call_variants = 70

    # create pool of workers and run through all samples
    call_variants_pool = ThreadPool(processes=n_call_variants)

    for c in cells:
        call_variants_pool.apply_async(resources.SingleCell.call_variants, args=(c,
                                                                                 human_fasta_file,
                                                                                 interval_file,))

    call_variants_pool.close()
    call_variants_pool.join()

    ################################################################################

    print 'Pipeline complete!'
