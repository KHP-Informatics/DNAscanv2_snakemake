import os
import sys

annovar_expansion_variants = sys.argv[1]
temp_expansion_variants = sys.argv[2]
annovar_STR_variants = sys.argv[3]
temp_STR_variants = sys.argv[4]
concise_results_report = sys.argv[5]

os.system("cat %s | awk '{print $1 \"\t\" $2 \"\t\" $9 \"\t\" $5 \"\t\" $7 \"\t\" $6 \"\t\" $10 \"\t\" $91 \"\t\" $79 \"\t\" $77 \"\t\" $78 \"\t\" $83 \"\t\" $120 \"\t\" $121}' | awk -v OFS='\t' '{$3=\"RepeatExpansion\" ; print ;}' | awk -v OFS='\t' '{split($9,a,/:/);$9=a[5]}1' | awk -v OFS='\t' ' {NR==1?$8=\"Overlapping_Genes\t\"$8:$8=\"\t\"$8 } 1 ' | awk -v OFS='\t' ' {NR==1?$11=\"OMIM_Phenotype\t\"$11:$11=\"\t\"$11 } 1 ' > %s" % (annovar_expansion_variants, temp_expansion_variants))
os.system("cat %s | awk '{print $1 \"\t\" $2 \"\t\" $9 \"\t\" $5 \"\t\" $7 \"\t\" $6 \"\t\" $10 \"\t\" $91 \"\t\" $79 \"\t\" $77 \"\t\" $78 \"\t\" $83 \"\t\" $120 \"\t\" $121}' | awk -v OFS='\t' '{$3=\"STR\" ; print ;}' | awk -v OFS='\t' '{split($9,a,/:/);$9=a[5]}1' | awk -v OFS='\t' ' {NR==1?$8=\"Overlapping_Genes\t\"$8:$8=\"\t\"$8 } 1 ' | awk -v OFS='\t' ' {NR==1?$11=\"OMIM_Phenotype\t\"$11:$11=\"\t\"$11 } 1 ' > %s" % (annovar_STR_variants, temp_STR_variants))
os.system("tail -n +2 %s >> %s" % (temp_STR_variants, temp_expansion_variants))
os.system("tail -n +2 %s >> %s" % (temp_expansion_variants, concise_results_report))
