import os
import sys

alsgenescanner_all = sys.argv[1]
alsgenescanner_alsod = sys.argv[2]
alsod_list = sys.argv[3]
alsgenescanner_clinvar = sys.argv[4]
clinvar_list = sys.argv[5]
alsgenescanner_manual_review = sys.argv[6]
manual_review_list = sys.argv[7]
alsgenescanner_ranked = sys.argv[8]

os.system("head -n 1 %s > %s ; cat %s | grep -iwf %s >> %s" % (alsgenescanner_all, alsgenescanner_alsod, alsgenescanner_all, alsod_list, alsgenescanner_alsod))
os.system("head -n 1 %s > %s ; cat %s | grep -iwf %s >> %s" % (alsgenescanner_all, alsgenescanner_clinvar, alsgenescanner_all, clinvar_list, alsgenescanner_clinvar))
os.system("head -n 1 %s > %s ; cat %s | grep -iwf %s >> %s" % (alsgenescanner_all, alsgenescanner_manual_review, alsgenescanner_all, manual_review_list, alsgenescanner_manual_review))
os.system("head -n 1 %s > %s ; cat %s | grep -i pathog | sed 's/ /_/g' | sort -k10nr >> %s ; cat %s | grep '^chr' | grep -iv pathog | sed 's/ /_/g' | sort -k10nr >> %s" % (alsgenescanner_all, alsgenescanner_ranked, alsgenescanner_all, alsgenescanner_ranked, alsgenescanner_all, alsgenescanner_ranked))
