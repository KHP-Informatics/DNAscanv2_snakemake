import os
import sys

non_human_stats = sys.argv[1]
non_human_list = sys.argv[2]
non_human_bam = sys.argv[3]
non_human_coverage_stats = sys.argv[4]
non_human_results = sys.argv[5]


os.system("awk \'{print $1}\' %s > %s ; for i in $(cat %s| grep ref); do printf \"$i \"; samtools depth %s -r $i | awk \'$3>0 {print $0}\' | wc -l; done > %s" % (non_human_stats, non_human_list, non_human_list, non_human_bam, non_human_coverage_stats))
coverage_stats = open('%s' % (non_human_coverage_stats), 'r')
coverage_stats_lines = coverage_stats.readlines()
coverage_stats_file = open('%s' % (non_human_stats), 'r')
coverage_stats_file_lines = coverage_stats_file.readlines()

i = 0

results = open('%s' % (non_human_results), 'w')
results.write("Id\tGenome_length\tNumber_of_reads\tCoverage\n")

while i < len(coverage_stats_lines):
    results.write("%s\t%s\t%s\t%s\n" % (coverage_stats_file_lines[i], coverage_stats_file_lines[i].split('\t')[1], coverage_stats_file_lines[i].split('\t')[2], coverage_stats_file_lines[i].split(' ')[1].strip()))
    i += 1

results.close()
