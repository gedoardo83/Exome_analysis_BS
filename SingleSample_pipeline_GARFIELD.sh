#Decompose and Normalize VCF file with vt tools
/results/mywork/Tools/vt/vt decompose -s input.vcf | /results/mywork/Tools/vt/vt decompose_blocksub - | /results/mywork/Tools/vt/vt normalize -r /results/referenceLibrary/tmap-f3/hg19/hg19.fasta - > output.vcf

#Apply GARFIELD predictions

#Convert to Annovar input single sample VCF
/results/mywork/Tools/annovar/convert2annovar.pl -format vcf4 input.vcf > output.avinput

#Convert to Annovar multisample VCF
/results/mywork/Tools/annovar/convert2annovar.pl -format vcf4old input.vcf > output.avinput

#Annotate avinput file
perl /results/mywork/Tools/annovar/table_annovar.pl --nastring NA --remove --protocol refGene,avsnp147,genomicSuperDups,1000g2015aug_all,exac03,kaviar_20150923,gnomad_exome,dbnsfp33a,clinvar_20170130,gwasCatalog,dbnsfp31a_interpro,revel --buildver hg19 --operation g,f,r,f,f,f,f,f,f,r,f,f input.avinput /results/mywork/Tools/annovar/humandb/

#Add additional info with custom script
perl /results/mywork/Resources/Scripts/AddInfoFromVCF_GARFIELD.pl input.avinput.hg19_multianno.txt input_GARFIELD.vcf /results/mywork/Resources/RVIS_ExAC_4KW.ok.txt /results/mywork/Resources/clinvar/gene_condition_source_id_110816_pathogenic /results/mywork/Resources/GDI_full_10282015.txt /results/mywork/Hi_FP_database/HiQ_recurrent.perfect.txt > output.INFO4.txt
