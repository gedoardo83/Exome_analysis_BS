##On single subjects VCF
#Normalize vars using vt
/results/mywork/Tools/vt/vt decompose -s input.vcf | /results/mywork/Tools/vt/vt decompose_blocksub - | /results/mywork/Tools/vt/vt normalize -r /results/referenceLibrary/tmap-f3/hg19/hg19.fasta - > output.vcf

#Apply GARFIELD prediction

#Bgzip and index each normalized file (not the GARFIELD file)
bgzip sample1.vcf
tabix -p vcf sample1.vcf.gz

##Trio analysis
#Merge single samples VCF
vcf-merge sample1.vcf.gz sample2.vcf.gz sample3.vcf.gz > Trio_multisample.vcf

#Add GARFIELD values
perl /results/mywork/Resources/Scripts/Add_GARFIELD_to_Multisample.pl Trio_multisample.vcf file1.GAR.vcf,file2.GAR.vcf,file3.GAR.vcfi > Trio_multisample.GAR.vcf

#Annotate with AF from 1000G and ExAc03
perl /results/mywork/Tools/annovar/table_annovar.pl Trio_multisample.GAR.vcf /results/mywork/Tools/annovar/humandb/ -buildver hg19 -remove -nastring . -protocol 1000g2015aug_all,exac03 -operation f,f -vcfinput 

#Segregation analysis
perl /results/mywork/Resources/Scripts/SegregationAnalysis_trio.pl -f 0.01 -s status_file.txt -m AR -i AnnotatedVCF.vcf

##For each candidate file
#Normalize vars using vt
/results/mywork/Tools/vt/vt decompose -s input.vcf | /results/mywork/Tools/vt/vt decompose_blocksub - | /results/mywork/Tools/vt/vt normalize -r /results/referenceLibrary/tmap-f3/hg19/hg19.fasta - > output.vcf

#Convert to Annovar multisample VCF
perl /results/mywork/Tools/annovar/convert2annovar.pl -format vcf4old input.vcf > output.avinput

#Annotate avinput file
perl /results/mywork/Tools/annovar/table_annovar.pl --nastring NA --remove --protocol refGene,avsnp147,genomicSuperDups,1000g2015aug_all,exac03,kaviar_20150923,gnomad_exome,dbnsfp33a,clinvar_20170130,gwasCatalog,dbnsfp31a_interpro,revel --buildver hg19 --operation g,f,r,f,f,f,f,f,f,r,f,f input.avinput /results/mywork/Tools/annovar/humandb/

#Add additional info with custom script
perl /results/mywork/Resources/Scripts/AddInfoFromVCF_GARFIELD.pl input.avinput.hg19_multianno.txt Trio_multisample.GAR.vcf /results/mywork/Resources/RVIS_ExAC_4KW.ok.txt /results/mywork/Resources/clinvar/gene_condition_source_id_110816_pathogenic /results/mywork/Resources/GDI_full_10282015.txt /results/mywork/Hi_FP_database/HiQ_recurrent.perfect.txt > output.INFO4.txt
