#Split VCF file in INDELS and SNPS
vcftools --vcf input.vcf --recode --recode-INFO-all --remove-indels --out output_SNP.vcf
vcftools --vcf input.vcf --recode --recode-INFO-all --keep-only-indels --out output_INDELS.vcf

#Filter low quality INDELS variants from VCF file
perl /results/mywork/Resources/Scripts/FilterVCF.pl -i input.vcf -o output.vcf -g 8 -d 20 -q 30 -h 5 -f 4 -s 1

#Filter low quality SNPS variants from VCF file
perl /results/mywork/Resources/Scripts/FilterVCF.pl -i input.vcf -o output.vcf -g 5 -d 6 -q 20 -h 99 -f 2 -s 0.9

#Merge INDELS and SNPS vcf and sort
vcf-concat SNP.vcf INDELS.vcf | vcf-sort > output.vcf

#Decompose and Normalize VCF file with vt tools
/results/mywork/Tools/vt/vt decompose -s input.vcf | /results/mywork/Tools/vt/vt decompose_blocksub - | /results/mywork/Tools/vt/vt normalize -r /results/referenceLibrary/tmap-f3/hg19/hg19.fasta - > output.vcf

#Convert to Annovar input single sample VCF
/results/mywork/Tools/annovar/convert2annovar.pl -format vcf4 input.vcf > output.avinput

#Convert to Annovar multisample VCF
/results/mywork/Tools/annovar/convert2annovar.pl -format vcf4old input.vcf > output.avinput

#Annotate avinput file
perl /results/mywork/Tools/annovar/table_annovar.pl --nastring NA --remove --protocol refGene,avsnp147,genomicSuperDups,1000g2015aug_all,exac03,kaviar_20150923,dbnsfp30a,clinvar_20160302,gwasCatalog,dbnsfp31a_interpro,mcap,eigen,revel --buildver hg19 --operation g,f,r,f,f,f,f,f,r,f,f,f,f input.avinput /results/mywork/Tools/annovar/humandb/

#Add additional info with custom script
perl /results/mywork/Resources/Scripts/AddInfoFromVCF_v4.pl input.avinput.hg19_multianno.txt input.vcf /results/mywork/Resources/RVIS_ExAC_4KW.ok.txt /results/mywork/Resources/clinvar/gene_condition_source_id_110816_pathogenic /results/mywork/Resources/GDI_full_10282015.txt /results/mywork/Hi_FP_database/HiQ_recurrent.perfect.txt > output.INFO4.txt
