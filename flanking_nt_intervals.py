#!/usr/bin/env python

import sys

if len(sys.argv) != 2:
    print("Must specify 1 filename")
    sys.exit(1)

def intervals(line_list):
    for file_line in line_list:
        first_num = int(file_line[(file_line.find(">") + 1):(file_line.find(":"))])
        second_num = int(file_line[(file_line.find(":") + 1):(file_line.find("<"))])
        corr_num_one = first_num - 1
        corr_num_two = second_num + 1
        gene = file_line[0:file_line.find(">")]
        ids = file_line[(file_line.find("<")+1):file_line.find("+")]
        pval = file_line[(file_line.find("+")+1):(len(file_line)-1)]
        if second_num > 1:
            print(gene + "\t" + str(corr_num_one) + "\t" + str(corr_num_two) + "\t" + ids + "\t" + pval)
        else:
            print(gene + "\t" + str(first_num) + "\t" + str(corr_num_two) + "\t" + ids + "\t" + pval)


bed_file = open(sys.argv[1], 'r')
lines = bed_file.readlines()
intervals(lines)

