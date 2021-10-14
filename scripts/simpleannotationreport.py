import os
import sys
import re

gene_list_file = open(sys.argv[2])

gene_list_lines = gene_list_file.readlines()

gene_list = gene_list_lines

out_file_all = open(sys.argv[3], 'w')

counter = 0

for i in gene_list:

    with open(sys.argv[1]) as vcf:

        for j in vcf:

            check1 = re.search(
                r'(^chr)|(^[0-9,X,Y,M]+\t)', j, flags=0)

            check = re.search(
                "=%s;" % (i.strip().upper()), j, flags=0)

            if check and check1:

                infos = j.split('ANNOVAR_DATE')[1][12:].split('ALLELE_END')[0].replace(";", "\t")

                if counter == 0:

                    replaced_1 = re.sub(
                        '=[a-z,A-Z,0-9,\.,\_,\-,:,>,<]+', '',
                        infos)

                    out_file_all.write(
                        'CHR\tPosition\tRef\tAlt\tGenotype\t%s\n' %
                        (replaced_1))

                    counter = 1

                replaced = re.sub('[a-z,A-Z,0-9,\.,\_,\-,:,>,<]+=',
                                        '', infos)

                out_file_all.write(
                    '%s\t%s\t%s\t%s\t%s\t%s\n' %
                    (j.split('\t')[0], j.split('\t')[1],
                     j.split('\t')[3], j.split('\t')[4],
                     j.split('\t')[-1].split(':')[0], replaced))

out_file_all.close()
