# Exome_analysis_BS
Exome analysis pipeline used at UniBS

## ExomePipeline_TrioQuad.pl
This perform full analysis pipeline for trios or quad, starting directly from single samples VCFs.
Select segregating variants according to inheritance model and generate annotated candidate variants files.

#### Options as follow:
- input -- Comma separated vcf.gz files for the trio / quad
- status  --  Status file (tab-delimited file with samplesID as in VCF file and 1=affected 0=unaffected)
- settings -- settings file containing path to relevant tools (see example in test folder)
- reference  --  fasta file of reference genome
- remove -- activate this option to remove temporary files
- maf  --  MAF threshold for filtering (default 0.01)
- output  --  Output folder
- model  --  Model of inheritacne (Possible values AD / AR)

`ExomePipeline_TrioQuad.pl --input mother.vcf,fater.vcf,proband.vcf --status family_status.txt --settings settings_file.txt --reference hg19.fasta --maf 0.01 --output results_folder/ --model AR --remove`
 
## Other scripts
Desribe single steps to perform different kinds of analysis manually
