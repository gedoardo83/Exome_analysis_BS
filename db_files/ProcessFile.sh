grep "Pathogenic\|pathogenic" variant_summary_110816.txt > variant_summary_110816_pathogenic.txt
head -n1 variant_summary_110816.txt | cat - variant_summary_110816_pathogenic.txt > variant_summary_110816_pathogenic.txt2
mv variant_summary_110816_pathogenic.txt2 variant_summary_110816_pathogenic.txt
cut -f5 variant_summary_110816_pathogenic.txt | sort -u | grep -v "^-" | grep -v "more than" > variant_summary_110816_pathogenic.genelist
grep -w --file variant_summary_110816_pathogenic.genelist gene_condition_id_110816 > gene_condition_id_110816_pathogenic
