##On single subjects VCF
#Split VCF files in INDELS and SNPS
vcftools --vcf input.vcf --recode --recode-INFO-all --remove-indels --out output_SNP.vcf
vcftools --vcf input.vcf --recode --recode-INFO-all --keep-only-indels --out output_INDELS.vcf

#Filter low quality INDELS variants from VCF file
perl /results/mywork/Resources/Scripts/FilterVCF.pl -i input.vcf -o output.vcf -g 8 -d 20 -q 30 -h 5 -f 4 -s 1

#Filter low quality SNPS variants from VCF file
perl /results/mywork/Resources/Scripts/FilterVCF.pl -i input.vcf -o output.vcf -g 5 -d 6 -q 20 -h 99 -f 2 -s 0.9

#Merge INDELS and SNPS vcf and sort
vcf-concat SNP.vcf INDELS.vcf | vcf-sort > output.vcf

#Bgzip and index each files
bgzip sample1.vcf
tabix -p vcf sample1.vcf.gz

##Trio analysis
#Merge single samples VCF
vcf-merge sample1.vcf.gz sample2.vcf.gz sample3.vcf.gz > Trio_multisample.vcf

#Annotate with AF from 1000G and ExAc03
/results/mywork/Tools/annovar/table_annovar.pl Trio_multisample.vcf /results/mywork/Tools/annovar/humandb/ -buildver hg19 -remove -nastring . -protocol 1000g2015aug_all,exac03 -operation f,f -vcfinput 

#Segregation analysis
/results/mywork/Resources/Scripts/SegregationAnalysis_trio.pl -f 0.01 -s status_file.txt -m AR -i AnnotatedVCF.vcf

##For each candidate file
#Decompose and Normalize VCF file with vt tools
/results/mywork/Tools/vt/vt decompose -s input.vcf | /results/mywork/Tools/vt/vt decompose_blocksub - | /results/mywork/Tools/vt/vt normalize -r /results/referenceLibrary/tmap-f3/hg19/hg19.fasta - > output.vcf

#Convert to Annovar multisample VCF
/results/mywork/Tools/annovar/convert2annovar.pl -format vcf4old input.vcf > output.avinput

#Annotate avinput file
perl /results/mywork/Tools/annovar/table_annovar.pl --nastring NA --remove --protocol refGene,avsnp147,genomicSuperDups,1000g2015aug_all,exac03,kaviar_20150923,gnomad_exome,dbnsfp33a,clinvar_20170130,gwasCatalog,dbnsfp31a_interpro,revel --buildver hg19 --operation g,f,r,f,f,f,f,f,f,r,f,f input.avinput /results/mywork/Tools/annovar/humandb/

#Add additional info with custom script
perl /results/mywork/Resources/Scripts/AddInfoFromVCF_v4.pl input.avinput.hg19_multianno.txt input.vcf /results/mywork/Resources/RVIS_ExAC_4KW.ok.txt /results/mywork/Resources/clinvar/gene_condition_source_id_110816_pathogenic /results/mywork/Resources/GDI_full_10282015.txt /results/mywork/Hi_FP_database/HiQ_recurrent.perfect.txt > output.INFO4.txt
