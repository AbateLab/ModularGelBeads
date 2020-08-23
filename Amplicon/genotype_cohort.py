'''

dab-seq: single-cell dna genotyping and antibody sequencing
ben demaree 7.9.2019

script for performing joint genotyping of a cohort (multiple samples and/or longitudinal data)
the main processing script must be run beforehand
caution: this script uses A LOT of memory - run one or two samples at a time

'''

import os
import subprocess
import time
import sys
import argparse
# import functions from external files
import resources


def wait(processes):
    # waits for processes to finish
    return [process.communicate() for process in processes]

def db_import(db_path, sample_map_path, interval):
    # import single-cell gvcfs for one interval

    # make call to gatk for genomics db import
    db_import_cmd = 'gatk --java-options "-Xmx4g" GenomicsDBImport ' \
                    '--genomicsdb-workspace-path %s ' \
                    '--batch-size 50 ' \
                    '--reader-threads 2 ' \
                    '--validate-sample-name-map true ' \
                    '-L %s ' \
                    '--sample-name-map %s'% (db_path,
                                             interval,
                                             sample_map_path)

    process = subprocess.Popen(db_import_cmd, shell=True)

    return process

def joint_genotype(db_path, fasta, interval, output_vcf):
    # perform joint genotyping across a cohort using data from a genomicsdb store

    # make call to gatk for genotyping
    genotype_cmd = 'gatk --java-options "-Xmx4g" GenotypeGVCFs ' \
                   '-V %s ' \
                   '-R %s ' \
                   '-L %s ' \
                   '-O %s ' \
                   '--include-non-variant-sites' \
                   % (db_path,
                      fasta,
                      interval,
                      output_vcf)

    process = subprocess.Popen(genotype_cmd, shell=True)

    return process

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='usage same as first script', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('sample_name', type=str, help='sample basename')
    parser.add_argument('cfg_file', type=str, help='config filename')
    
    args = parser.parse_args()  # parse arguments

    sample_basename = args.sample_name
    cfg_f = args.cfg_file
    
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

    ####################################################################################
    # set variables and filepaths for sample cohorts
    ####################################################################################


    cohort_name = 'genotyping'
    output_dir = base_dir + cohort_name + '/'

    # create cohort folder if it does not exist
    # note that GATK will not overwrite existing genomicsdbs
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # option to analyze flt3-itd variants
    flt3_itd = True

    # minimum read depth and single-cell vaf for all flt3 variants
    # this helps filter false positives
    min_itd_dp = 15
    min_itd_vaf = 0.15

    # list of sample ids
    sample_ids = 0

    # sample map for gatk
    sample_map = output_dir + cohort_name + '.sample_map.tsv'

    # list of directories containing single-sample gvcfs
    gvcf_paths = base_dir + 'Cells/gvcf/'

    # list of directories containing flt3-itd vcf files
    itd_paths = base_dir + 'Cells/flt3/'

    # interval file path (EXCLUDING primer coordinates)
    interval_file = sys.path[0] + '/amplicon_regions/AML.bed'

    # human reference genome fasta file path
    human_fasta = '/drive2/igenomes/hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa'

    # snpeff config file path
    snpeff_config = '/usr/local/bin/snpEff/snpEff.config'

    # snpeff summary file path
    snpeff_summary = output_dir + cohort_name + '.snpeff_summary.html'

    # genotyped vcf file
    genotyped_vcf = output_dir + cohort_name + '.genotyped.vcf'

    # split vcf file
    split_vcf = output_dir + cohort_name + '.split.vcf'

    # snpeff annotated vcf
    snpeff_annot_vcf = output_dir + cohort_name + '.snpeff.annotated.vcf'

    # final annotated vcf
    annot_vcf = output_dir + cohort_name + '.annotated.vcf'

    # genotype hdf5 file
    geno_hdf5 = output_dir + cohort_name + '.genotypes.hdf5'

    # variant information table
    variants_tsv = output_dir + cohort_name + '.variants.tsv'

    # combined flt3 vcf
    flt3_vcf = output_dir + cohort_name + '.flt3.vcf'

    # clinvar vcf file path
    clinvar_vcf = '/drive2/snp_dbs/clinvar_20190805_hg19/clinvar_20190805.hg19.vcf.gz'

    # cosmic vcf file path
    cosmic_vcf = '/drive2/snp_dbs/cosmic_v89_hg19/CosmicAll_v89.hg19.vcf.gz'

    ####################################################################################
    # prepare gvcf and interval files for genotyping
    ####################################################################################

    # send slack notification
    start_time = time.time()

    # base path for genomics db
    db_dir = output_dir + 'dbs/'
    if not os.path.exists(db_dir):
        os.mkdir(db_dir)

    # base path for single-interval vcfs
    vcf_dir = output_dir + 'vcfs/'
    if not os.path.exists(vcf_dir):
        os.mkdir(vcf_dir)

    # list of lists of gvcf files for each sample
    gvcf_files = []
    gvcf_files.append([f for f in os.listdir(gvcf_paths) if f.endswith('.g.vcf')])

    # optional: list of lists of flt3 vcf files for each sample
    if flt3_itd:
        itd_files = []
        itd_files.append([itd_paths + f for f in os.listdir(itd_paths) if f.endswith('.vcf')])

    # create sample map file
    with open(sample_map, 'w') as f:
        for g in gvcf_files[0]:
            f.write(g.split('.')[0] + '-0' + '\t' + gvcf_paths + g + '\n')

    # extract intervals from bed file
    # option to exclude RUNX1_4 (used for antibodies)
    exclude_RUNX1_4 = True
    intervals = {}
    with open(interval_file, 'r') as f:
        for line in f:
            if exclude_RUNX1_4 and 'RUNX1_4' not in line:
                fields = line.strip().split('\t')
                intervals[fields[3]] = fields[0] + ':' + fields[1] + '-' + fields[2]

    # create a genomics db and output vcf for each interval
    db_paths = {}
    output_vcfs = {}
    for L in intervals:
        db_paths[L] = db_dir + L + '.genomics.db'
        output_vcfs[L] = vcf_dir + L + '.genotyped.vcf'

    ####################################################################################
    # import gvcfs into genomics DB
    ####################################################################################

    # parallelize by starting a process for each interval
    import_processes = []

    for L in intervals:
        import_processes.append(db_import(db_paths[L],
                                          sample_map,
                                          intervals[L]))

    # wait for processes to finish
    wait(import_processes)

    ####################################################################################
    # perform joint genotyping
    ####################################################################################

    # parallelize by starting a process for each interval
    genotype_processes = []

    for L in intervals:
        genotype_processes.append(joint_genotype('gendb://' + db_paths[L],
                                                 human_fasta,
                                                 intervals[L],
                                                 output_vcfs[L]))

    # wait for processes to finish
    wait(genotype_processes)

    ####################################################################################
    # merge single-interval vcfs
    ####################################################################################

    # call gatk to perform vcf merging
    vcfs_to_merge = ['-I ' + v for v in output_vcfs.values()]
    merge_cmd = 'gatk MergeVcfs %s -O %s' % (' '.join(vcfs_to_merge), genotyped_vcf)

    subprocess.call(merge_cmd, shell=True)

    ####################################################################################
    # split multiallelic sites and annotate vcf
    ####################################################################################

    # split multiallelics, left-align, and trim
    resources.left_align_trim(human_fasta, genotyped_vcf, split_vcf)

    # annotate vcf with snpeff (functional predictions)
    resources.snpeff_annotate(snpeff_summary, snpeff_config, split_vcf, snpeff_annot_vcf)

    # annotate with bcftools
    # use clinvar database
    resources.bcftools_annotate(clinvar_vcf, snpeff_annot_vcf, '-c INFO', annot_vcf)

    # convert vcf to variant matrix in hdf5 format
    # optional: include itd calls

    if flt3_itd:
        # combine flt3 itd vcf files from all samples
        # only include variants with sufficient read depth

        with open(flt3_vcf, 'w') as v:
            header = True
            for itd in itd_files[0]:
                with open(itd, 'r') as f:
                    for line in f:
                        if header and line[:2] == '##':
                            v.write(line)
                            continue

                        elif header and line[:2] == '#C':
                            v.write(line)

                        header = False
                        cell_barcode = itd.split('.')[0].split('/')[-1] + '-0'

                        if line[0] != '#':
                            depth = int(line.split('\t')[5])
                            vaf = float(line.strip().split('=')[-1])
                            if depth >= min_itd_dp and vaf >= min_itd_vaf:
                                # add cell barcode to id column
                                vcf_record = '\t'.join(line.split('\t')[:2] + [cell_barcode] + line.split('\t')[3:])
                                v.write(vcf_record)
                            else:
                                break

        resources.vcf_to_tables(annot_vcf, geno_hdf5, variants_tsv, flt3_vcf)

    else:
        resources.vcf_to_tables(annot_vcf, geno_hdf5, variants_tsv)

    # compress the final annotated vcf to reduce size
    subprocess.call('gzip %s' % annot_vcf, shell=True)

    ####################################################################################
    # clean up temporary files
    ####################################################################################

    try:
        os.remove(sample_map)
        os.remove(genotyped_vcf)
        os.remove(genotyped_vcf + '.idx')
        os.remove(split_vcf)
        os.remove(snpeff_annot_vcf + '.gz')
        os.remove(snpeff_annot_vcf + '.gz.tbi')

    except OSError:
        pass

    ####################################################################################

    print 'Pipeline complete!'

    # send slack notification
    elapsed_time = time.time() - start_time
    elapsed_time_fmt = str(time.strftime('%Hh %Mm %Ss', time.gmtime(elapsed_time)))
    print elapsed_time_fmt
