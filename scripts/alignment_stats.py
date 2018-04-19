import argparse
import os
import time
import logging as log
import pandas as pd
import yaml
import subprocess
#import matplotlib

from os import listdir
from os.path import isfile, join

#matplotlib.use('Agg')  # non-interactive backend
#import matplotlib.pyplot as plt

#import seaborn as sns
import numpy as np
from string import Template
from subprocess import Popen, PIPE

import pysam
import pybedtools
from scipy import stats

def getColumn(col_num, data):
    column = []
    for row in data:
        row_val = int(row[col_num])
        column.append(row_val)
    return column

def calculate_median_depth(c):
    column = getColumn(5 , c)
    return np.median(np.array(column))

def calculate_average_depth(c):
    count = 0;
    total = len(c);
    for row in c:
        count = count + int(row[5])

    average = float(count*1.0/total*1.0)
    return average

def calculate_zero_depth_intervals(c, depth):
    count = 0;
    total = len(c);
    for row in c:
        if(int(row[5]) >= 0 and int(row[5]) <= depth):
            count = count + 1

    percent = (count*1.0/total)*100.0
    #return the percentage
    return percent

def calculate_x_depth_intervals(c, depth):
    count = 0;
    total = len(c);
    for row in c:
        if(int(row[5]) >= depth):
            count = count + 1

    percent = (count*1.0/total)*100.0
    #return the percentage
    return percent

def calculate_x_depth_intervals_folds(c, start_depth, end_depth):
    count = 0;
    total = len(c);
    for row in c:
        if(int(row[5]) >= start_depth and int(row[5]) <= end_depth):
            count = count + 1

    percent = (count*1.0/total)*100.0
    #return the percentage
    return percent


def main():
    """The main function"""
    current_dir = os.getcwd()

    #mypath = current_dir + "/bams/"
    mypath = current_dir + "/alignments/"
    #bed_file = "/vlsci/UOM0040/shared/djp/rover_file/crc/CRC_10g_23May16.final.rover.bed"
    fastq_dir = current_dir + "/fastqs/"
    config_file = "pipeline.config"
    with open(config_file, 'r') as stream:
        try:
            bed_file = yaml.load(stream)['target_bed']
        except yaml.YAMLError as exc:
            print("Error with config file: " + exc)

    onlyfiles = []

    for root, dirs, files in os.walk(mypath):
        for file in files:
            if file.endswith(".primary.primerclipped.bam"):
                current_file = mypath + str(file)
                onlyfiles.append(os.path.join(root, file))

    #onlyfiles = [files for files in listdir(mypath) if (isfile(join(mypath, files)) and (files.endswith('.primary.primerclipped.bam')))]
    #file_paths = [join(mypath,files) for files in listdir(mypath) if (isfile(join(mypath, files)) and (files.endswith('.bam')))]
    #onlyfiles = [files for files in listdir(mypath) if (files.endswith(''))]

    #print onlyfiles
    #print len(onlyfiles)

    # stats list
    header = '\t'.join([ 'Sample_ID', 'Total_fastq_reads', 'Primary_reads', 'Reads_mapping_to_genome' , 'Reads_mapping_to_target', 'Percent_reads_mapping_to_genome', 'Percent_reads_mapping_to_target', 'Average_depth', \
        'Percent_target_not_covered', 'Percent_target_covered_at_<10X', 'Percent_target_covered_at_10X', 'Percent_target_covered_at_20X', 'Percent_target_covered_at_50X', 'Median_depth', \
        'Percent_target_covered_at_median', \
        'Percent_target_covered_at_median_10_fold', 'Percent_target_covered_at_median_20_fold', 'Percent_target_covered_at_median_30_fold', \
        'Percent_target_covered_at_median_40_fold', 'Percent_target_covered_at_median_50_fold'])
        #, 'Percent_target_covered_at_q50', \
        #'Percent_target_covered_at_q60', 'Percent_target_covered_at_q70', 'Percent_target_covered_at_q80'])

    #header = "Sample\tTotal_reads\tMapped_reads"

    print header

    for bam_file in onlyfiles:
        current_bam_file = join(mypath, bam_file)
        temp_bam_file = os.path.basename(current_bam_file)
        sample = temp_bam_file.replace(".primary.primerclipped.bam", "")

        fastq1 = fastq_dir + sample + "_L01_R1_001.fastq"
        fastq2 = fastq_dir + sample + "_L01_R2_001.fastq"


        fastq1_lc = int(subprocess.check_output(["wc", "-l", fastq1]).lstrip(' ').split(' ')[0])
        fastq2_lc = int(subprocess.check_output(["wc", "-l", fastq2]).lstrip(' ').split(' ')[0])

        total_fastq_lines = fastq1_lc + fastq2_lc
        total_fastq_reads = total_fastq_lines / 4

        flagstats = pysam.flagstat(current_bam_file)
        all_reads = int(flagstats.split('\n')[0].split('+')[0])
        reads_mapping_to_genome = int(flagstats.split('\n')[5].split('+')[0])

        x = pybedtools.example_bedtool(current_bam_file)
        b = pybedtools.example_bedtool(bed_file)
        y = x.intersect(b).moveto(join(mypath, 'temp.bam'))
        c = b.coverage(x)

        average_depth = calculate_average_depth(c)
        median_depth = calculate_median_depth(c)

        percent_target_not_covered = calculate_zero_depth_intervals(c, 0)
        percent_target_covered_at_L10X = calculate_zero_depth_intervals(c, 10)

        percent_target_covered_at_10X = calculate_x_depth_intervals(c, 10)
        percent_target_covered_at_20X = calculate_x_depth_intervals(c, 20)
        percent_target_covered_at_50X = calculate_x_depth_intervals(c, 50)


        percent_target_covered_at_median = calculate_x_depth_intervals_folds(c, median_depth, median_depth)
        # Using percentage from median
        #percent_target_covered_at_median_X10 = calculate_x_depth_intervals_folds(c, (median_depth - median_depth * (10.0/100)), (median_depth + median_depth * (10.0/100)))
        '''
        percent_target_covered_at_median_10_fold = calculate_x_depth_intervals_folds(c, (median_depth - (median_depth * 0.10)), (median_depth + (median_depth * 0.10)))
        percent_target_covered_at_median_20_fold = calculate_x_depth_intervals_folds(c, (median_depth - (median_depth * 0.20)), (median_depth + (median_depth * 0.20)))
        percent_target_covered_at_median_50_fold = calculate_x_depth_intervals_folds(c, (median_depth - (median_depth * 0.50)), (median_depth + (median_depth * 0.50)))
        percent_target_covered_at_median_60_fold = calculate_x_depth_intervals_folds(c, (median_depth - (median_depth * 0.60)), (median_depth + (median_depth * 0.60)))
        percent_target_covered_at_median_70_fold = calculate_x_depth_intervals_folds(c, (median_depth - (median_depth * 0.70)), (median_depth + (median_depth * 0.70)))
        percent_target_covered_at_median_80_fold = calculate_x_depth_intervals_folds(c, (median_depth - (median_depth * 0.80)), (median_depth + (median_depth * 0.80)))
        '''


        percent_target_covered_at_median_10_fold = calculate_x_depth_intervals(c, (median_depth / (10)))
        percent_target_covered_at_median_20_fold = calculate_x_depth_intervals(c, (median_depth / (20)))
        percent_target_covered_at_median_30_fold = calculate_x_depth_intervals(c, (median_depth / (30)))
        percent_target_covered_at_median_40_fold = calculate_x_depth_intervals(c, (median_depth / (40)))
        percent_target_covered_at_median_50_fold = calculate_x_depth_intervals(c, (median_depth / (50)))


        stats_temp = pysam.flagstat(join(mypath,'temp.bam'))
        on_target_reads = int(stats_temp.split('\n')[0].split('+')[0])
        reads_mapping_to_target = int(stats_temp.split('\n')[5].split('+')[0])

        #percent_reads_mapping_to_genome = ((reads_mapping_to_genome * 1.0)/all_reads)*100.0
        percent_reads_mapping_to_genome = ((reads_mapping_to_genome * 1.0)/total_fastq_reads)*100.0
        #percent_reads_mapping_to_target = ((reads_mapping_to_target * 1.0)/on_target_reads)*100.0
        percent_reads_mapping_to_target = ((reads_mapping_to_target * 1.0)/total_fastq_reads)*100.0

        os.remove(join(mypath,'temp.bam'))

        print("%s\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f" % (sample, total_fastq_reads, all_reads, reads_mapping_to_genome, reads_mapping_to_target, \
            percent_reads_mapping_to_genome, percent_reads_mapping_to_target, average_depth, \
            percent_target_not_covered, percent_target_covered_at_L10X, percent_target_covered_at_10X, percent_target_covered_at_20X, percent_target_covered_at_50X, \
            median_depth, percent_target_covered_at_median, \
            percent_target_covered_at_median_10_fold, percent_target_covered_at_median_20_fold, percent_target_covered_at_median_30_fold, \
            percent_target_covered_at_median_40_fold, percent_target_covered_at_median_50_fold))

        #print bam_file + "\t" + str(all_reads) + "\t" + str(reads_mapping_to_genome) + "\t" + str(reads_mapping_to_target) + "\t" + \
        #    str(percent_reads_mapping_to_genome) + "\t" + str(percent_reads_mapping_to_target) + "\t" + str(average_depth) + "\t" + \
        #    str(percent_target_not_covered) + "\t" + str(percent_target_covered_at_10X) + "\t" + str(percent_target_covered_at_20X) + "\t" + str(percent_target_covered_at_50X)

if __name__ == '__main__':
    main()
