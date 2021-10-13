import os
import sys

annotsv_file = sys.argv[1]
temp_SV_MEI_variants = sys.argv[2] - out/reports/temp_sample_SV_MEI_variants.tsv
annovar_results_report = sys.argv[3]
temp_SNPindel_variants = sys.argv[4]
report_header = sys.argv[5]
concise_results_report = sys.argv[6]

os.system("cat %s | cut -f 2,3,6,15,18,28,29,34,73,87,88 | awk -v OFS='\t' '{split($4,a,/:/);$4=a[1]}1' | awk -v OFS='\t' ' {NR==1?$11=\"Clinvar_ID\t\"$11:$11=\"\t\"$11 } 1 ' | awk  -v OFS='\t' ' {NR==1?$13=\"Clinvar_Phenotype\t\"$13:$13=\"\t\"$13 } 1 ' | awk  -v OFS='\t' ' {NR==1?$14=\"Variant_Frequency_ExAC\t\"$14:$14=\"\t\"$14 } 1 ' | awk  -v OFS='\t' ' {NR==1?$14=\"Variant_Frequency_1000g\t\"$14:$14=\"\t\"$14 } 1 ' | awk  -v OFS='\t' ' {NR==1?$14=\"Variant_Frequency_gnomAD\t\"$14:$14=\"\t\"$14 } 1 ' | awk -F '\t' 'NR>1 {print \"chr\"$0}' > %s" % (annotsv_file, temp_SV_variants))

os.system("cat %s | awk '{print $1 \"\t\" $2 \"\t\" $9 \"\t\" $5 \"\t\" $7 \"\t\" $6 \"\t\" $10 \"\t\" $91 \"\t\" $79 \"\t\" $77 \"\t\" $78 \"\t\" $83 \"\t\" $120 \"\t\" $121}' | awk -v OFS='\t' '{$3=\"small_variant\" ; print ;}' | awk -v OFS='\t' '{split($9,a,/:/);$9=a[5]}1' | awk -v OFS='\t' ' {NR==1?$8=\"Overlapping_Genes\t\"$8:$8=\"\t\"$8 } 1 ' | awk -v OFS='\t' ' {NR==1?$11=\"OMIM_Phenotype\t\"$11:$11=\"\t\"$11 } 1 ' > %s" % (annovar_results_report, temp_SNPindel_variants))

os.system("cp %s %s/reports/%s_all_variants.tsv" % (report_header, concise_results_report))
os.system("tail -n +2 %s >> %s" % (temp_SNPindel_variants, concise_results_report))
os.system("tail -n +2 %s >> %s" % (temp_SV_variants, concise_results_report))
