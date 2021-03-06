#!usr/bin/perl
#My full workflow for examine TRIOS or small pedigree data.

use Getopt::Long qw(GetOptions);
use File::Basename;

#####################
###  Get options  ###
#####################

#Set default values for options
$model= "AR" ;
$outpath= ".\\" ;
$maf=0.01;

#Parsing command line options
#Mandatory options
GetOptions(
	'input=s' => \$vcfs, #Comma separated vcf.gz files for the trio / quad
	'status=s' => \$pedfile, #Status file (tab-delimited file with samplesID as in VCF file and 1=affected 0=unaffected)
	'settings=s' => \$setfile,
	'reference=s' => \$reference,
	'remove!' => \$notemp,
	'maf:f' => \$maf, #MAF threshold for filtering
	'output:s' => \$outpath, #Output folder
	'model:s' => \$model, #Possible values AD / AR
) or die "Uncorrect or missing mandatory arguments\n";

$start_run=time();
print "Arguments as interpreted:\n";
print "\tInput files: $vcfs\n";
print "\tStatus file: $pedfile\n";
print "\tSettings file: $setfile\n";
print "\tReference genome: $reference\n";
print "\tMAF Threshold: $maf\n";
print "\tInheritance model: $model\n";
print "\tOutput folder: $outpath\n\n";

#Check input files existance and output files
if (!-f $setfile || !-f $pedfile) {die "FATAL! File $setfile does not exist!"}
if (!-f $pedfile) {die "FATAL! File $pedfile does not exist!"}
if (!-f $reference) {die "FATAL! File $reference does not exist!"}

if ($model eq "AR") {
  if (-f $outpath."/Recessive_dandidates.vcf" ||
      -f $outpath."/Compound_het_dandidates.vcf" ||
      -f $outpath."/DeNovo_variants.vcf") {die "FATAL! Output files already exist in $outpath folder!"}
} elsif ($model eq "AD") {
  if (-f $outpath."/Dominant_dandidates.vcf") {die "FATAL! Output files already exist in $outpath folder!"}
}
if (-f $outpath."/Multisample.vcf") {die "FATAL! File Multisample.vcf already exist in $outpath folder!"}

@vcffiles = split(",",$vcfs);
foreach (@vcffiles) {
  if (!-f $_) {die "FATAL! Input VCF file $_ does not exist!"}
}

system("which bgzip > /dev/null 2>&1");
if ($? != 0) {die "FATAL! Unable to launch bgzip. It is not installed or not in your path. bgzip needed to open .vcf.gz files\nError code: $!\n"}
system("which tabix > /dev/null 2>&1");
if ($? != 0) {die "FATAL! Unable to launch tabix. It is not installed or not in your path. tabix needed to index .vcf.gz files\nError code: $!\n"}

#######################################
###  Read Settings and Status file  ###
#######################################

#Reading resources locations from setting file
print "Reading settings file...\n";
open(IN, $setfile);
while($row=<IN>) {
  chomp($row);
  @line=split("=",$row);
  $line[0] =~ s/^\s+|\s+$//g;
  $line[1] =~ s/^\s+|\s+$//g;
  $settings{$line[0]}=$line[1];
}

#Check if all required values are present
@expected = ("annovar", "annovardb", "gdi", "rvis", "vcftools","vt", "clinvar", "platform"); #additional vardb
foreach (@expected) {
  if (exists $settings{$_}) {$count++} else {$missingset .= " $_"}
}
die "FATAL! Missing one of required settings:".$missingset." !\n" if ($count<8);

#Additional settings could be vardb and GARFIELD
#platform must be ion / illumina

#Read sample status file
print "INFO\tReading status file...\n";
open(PED, $pedfile);
while($row=<PED>) {
	chomp($row);
	@line = split("\t", $row);
	if ($line[1]==0) {$control{$line[0]}=0}
	if ($line[1]==1) {$affected{$line[0]}=0}
}
close(PED);

#####################################
###  Pre-process input VCF files  ###
#####################################

### STEP 1.Decompose and Normalize candidate VCF files with vt tools
print "INFO\tNormalizing variants with vt tools...\n";
foreach (@vcffiles) {
  my($vcfname, $vcfpath, $vcfsuffix) = fileparse($_);
  $command="$settings{vt} decompose -s $_ | $settings{vt} decompose_blocksub - | $settings{vt} normalize -r $reference - > $outpath/$vcfname.NORM.TEMP.vcf";
  print "INFO\tNormalizing vcf: $command\n";
  system($command);

  ##GARFIELD annotation if required
  if (exists $settings{GARFIELD}) {
    $command="perl $settings{GARFIELD} --input $outpath/$vcfname.NORM.TEMP.vcf --platform $settings{platform} --out $outpath/$vcfname.GAR.TEMP.vcf";
    print "INFO\tGARFIELD prediction: $command\n";
    system($command);
    push(@garfiles, "$outpath/$vcfname.GAR.TEMP.vcf");
  }
  
  system("bgzip $outpath/$vcfname.NORM.TEMP.vcf");
  system("tabix -p vcf $outpath/$vcfname.NORM.TEMP.vcf.gz");
  push(@normfiles, "$outpath/$vcfname.NORM.TEMP.vcf.gz");
}

### STEP2. Merge single sample VCF files into Multisample.vcf
$command= $settings{vcftools}."/vcf-merge ".join(" ", @normfiles)." > ".$outpath."/Multisample.vcf";
print "INFO\tRunning vcf-merge: $command\n";
system($command);


$command="perl $settings{annovar}/table_annovar.pl --nastring . --remove --protocol refGene,1000g2015aug_all,exac03 --vcfinput --buildver hg19 --operation g,f,f $outpath/Multisample.vcf $settings{annovardb}";
print "INFO\tRunning vcf annotation: $command\n";
system($command);
system("mv $outpath/Multisample.vcf.hg19_multianno.vcf $outpath/Multisample.vcf");

### STEP2. Add GARFIELD annotation if required
if (exists $settings{GARFIELD}) {
print "INFO\tAdding GARFIELD score to multisample\n";
  foreach (@garfiles) {
  	open(VCF, $_);
    while ($row=<VCF>) {
	    next if ($row =~ /^##/); #Skip header lines
	    if ($row =~ /^#/) {
		    chomp($row);
		    @line = split("\t", $row);
		    $sampleid = $line[9]."_true";
		    next;
	    }
		
	  chomp($row);
	  @line = split("\t", $row);
	  $varid = $line[0]."_".$line[1]."_".$line[3]."_".$line[4];
		
	  #Extract GARFIELD score
	  $line[7] =~ /($sampleid=[0-9.]+)/;
	  push(@{$GARscore{$varid}}, $1);
    }
  close(VCF);
  }

  open(IN, "$outpath/Multisample.vcf");
  open(OUT, ">>$outpath/Multisample.temp");
  while ($row=<IN>) {
	  if ($row =~ /^#/) {
		  print OUT $row; #Print out header lines
		  next;
	  }
		
	  chomp($row);
	  @line = split("\t", $row);
	  $varid = $line[0]."_".$line[1]."_".$line[3]."_".$line[4];
		
	  #Add GARFIELD scores to INFO column
	  $line[7] .= ";".join(";", @{$GARscore{$varid}});
	  print OUT join("\t", @line)."\n";
  }
  close(IN);
  close(OUT);
  system("mv $outpath/Multisample.temp $outpath/Multisample.vcf")
}

####################################
###  Select candidates variants  ###
####################################
#Read Multisample.vcf, select candidates vars and write them to files

#Define columns for affected and unaffected subjects
print "INFO\tReading multisample VCF file...\n";

$headline=`grep -m1 -n "^#CHROM" $outpath/Multisample.vcf`;
chomp($headline);
@line=split("\t",$headline);
@samplesid=@line[9..$#line];
for ($i=9; $i<=$#line; $i++) {
	if (exists $control{$line[$i]}) {push(@cntrindex, $i)}
	if (exists $affected{$line[$i]}) {push(@affectedindex, $i)}
}

#Store header of input VCF
$headline =~ /([0-9]+):#CHROM/;
$command= "sed '".$1."q' $outpath/Multisample.vcf";
$header=`$command`;

#Read variants and store informations
open (VCF, "$outpath/Multisample.vcf");
while ($row = <VCF>) {
  if ($row =~ /^##/) {next}
  if ($row =~ /^#CHROM/) {
    chomp($row);
    @line=split("\t", $row);
    @samplesid=@line[9..$#line];
    next;
  }
  
	$step++;
	print "\rReading variant $step...";
	chomp($row);
	@line = split("\t", $row);
	
	#Save entire var row
	$varid=$line[0]."_".$line[1]."_".$line[3]."_".$line[4];
	$var{$varid}{data}=$row;
  
  #Save AF
  if ($line[7] =~ /1000g2015aug_all=([01]\.[0-9]+)/ || $line[7] =~ /1000g2015aug_all=(1)/) {$var{$varid}{AF1KG}=$1} else {$var{$varid}{AF1KG}=0}
  if ($line[7] =~ /ExAC_ALL=([01]\.[0-9]+)/ || $line[7] =~ /ExAC_ALL=(1)/) {$var{$varid}{AFEXAC}=$1} else {$var{$varid}{AFEXAC}=0}
	
	#Save gene name
	$line[7] =~ /Gene.refGene=([^;]+);/;
	$var{$varid}{GENE} = $1;
	
	#Save genotypes
	foreach $myindex(@cntrindex) {
		@tags = split(":",$line[$myindex]);
		push (@{$var{$varid}{controls}}, $tags[0]);
	}
	foreach $myindex(@affectedindex) {
		@tags = split(":",$line[$myindex]);
		push (@{$var{$varid}{affected}}, $tags[0]);
	}
}
close(VCF);
$total=$step;
print "\nINFO\t$total variants loaded\n";

#Processing variants
print "INFO\tVariant analysis started\n";

$step=0;
if ($model eq "AR") {
  open(RECESSIVE, ">>$outpath/Recessive_candidates.vcf");
  open(COMPOUND, ">>$outpath/Compound_het_candidates.vcf");
  open(DENOVO, ">>$outpath/DeNovo_variants.vcf");
  print RECESSIVE "$header";
  print COMPOUND "$header";
  print DENOVO "$header";
  @outfiles = ("$outpath/Recessive_candidates.vcf", "$outpath/Compound_het_candidates.vcf", "$outpath/DeNovo_variants.vcf");
} elsif ($model eq "AD") {
  open(DOMINANT, ">>$outpath/Dominant_candidates.vcf");
  open(DENOVO, ">>$outpath/DeNovo_variants.vcf");
  print DOMINANT "$header";
  print DENOVO "$header";
  @outfiles = ("$outpath/Dominant_candidates.vcf","$outpath/DeNovo_variants.vcf");
}

foreach $myvar(keys %var) {
	$n=-1;
	$p=-1;
	$step++;
	print "\rAnalyzing segregation on variant $step of $total";

#AF filter
	if ($var{$myvar}{AF1KG} < $maf && $var{$myvar}{AFEXAC} < $maf) {
    
#Processing for Recessive model
    if ($model eq "AR") {

			foreach $mygeno(@{$var{$myvar}{affected}}) {
				if ($mygeno eq "1/1") {$n++}
			}
			foreach $mygeno(@{$var{$myvar}{controls}}) {
				if ($mygeno eq "0/1" || $mygeno eq "0/0" || $mygeno eq "\." || $mygeno eq "\./\.") {$p++}
			}
			if ($n==$#affectedindex && $p==$#cntrindex) {print RECESSIVE "$var{$myvar}{data}\n"}

			$n=-1;
			$p=-1;
	
			foreach $mygeno(@{$var{$myvar}{affected}}) {
				if ($mygeno eq "0/1") {$n++}
			}
			foreach $mygeno(@{$var{$myvar}{controls}}) {
				if ($mygeno eq "0/0" || $mygeno eq "\." || $mygeno eq "\./\.") {$p++}
			}
			if ($n==$#affectedindex && $p==$#cntrindex) {print DENOVO "$var{$myvar}{data}\n"}
      
      		$n=-1;
			$p=-1;
			
			foreach $mygeno(@{$var{$myvar}{affected}}) {
				if ($mygeno eq "0/1") {$n++}
			}
			foreach $mygeno(@{$var{$myvar}{controls}}) {
				if ($mygeno eq "0/1" || $mygeno eq "0/0" || $mygeno eq "\." || $mygeno eq "\./\.") {$p++}
			}
			if ($n==$#affectedindex && $p==$#cntrindex) {
			  if ($var{$myvar}{GENE} ne ".") {push (@{$compound{$var{$myvar}{GENE}}}, $var{$myvar}{data})}
			}
			 
    
#Processing for Dominant model
		} elsif ($model eq "AD") {

			foreach $mygeno(@{$var{$myvar}{affected}}) {
				if ($mygeno eq "0/1") {$n++}
			}
			foreach $mygeno(@{$var{$myvar}{controls}}) {
				if ($mygeno eq "0/0" || $mygeno eq "\." || $mygeno eq "\./\.") {$p++}
			}
			if ($n==$#affectedindex && $p==$#cntrindex) {print DOMINANT "$var{$myvar}{data}\n"}

			$n=-1;
			$p=-1;
	
			foreach $mygeno(@{$var{$myvar}{affected}}) {
				if ($mygeno eq "0/1") {$n++}
			}
			foreach $mygeno(@{$var{$myvar}{controls}}) {
				if ($mygeno eq "0/0" || $mygeno eq "\." || $mygeno eq "\./\.") {$p++}
			}
			if ($n==$#affectedindex && $p==$#cntrindex) {print DENOVO "$var{$myvar}{data}\n"}
		}
	}
}

#Print potential compound heterozygous for recessive model
if ($model eq "AR") {
  foreach (keys %compound) {
    if (@{$compound{$_}} > 1) {
    print COMPOUND join("\n", @{$compound{$_}})."\n";
  }
}
}

if ($model eq "AR") {
  close(RECESSIVE);
  close(COMPOUND);
  close(DENOVO);
} elsif ($model eq "AD") {
  close(DOMINANT);
  close(DENOVO);
}

#######################################################
###  Process and annotate candidate variants files  ###
#######################################################

### STEP 1.Annotating candidate variants with Annovar
#Convert to Annovar
print "INFO\tConverting to Annovar input...\n";

foreach (@outfiles) {
  $command="$settings{annovar}/convert2annovar.pl --includeinfo -format vcf4old $_ > $_.TEMP.avinput";
  system($command);
}

#Annotate avinput file
print "INFO\tAnnotating with Annovar...\n";
foreach (@outfiles) {
  $command="perl $settings{annovar}/table_annovar.pl --nastring NA --remove --otherinfo --protocol refGene,avsnp147,genomicSuperDups,1000g2015aug_all,exac03,kaviar_20150923,gnomad_exome,dbnsfp33a,clinvar_20170130,gwasCatalog,dbnsfp31a_interpro,revel --buildver hg19 --operation g,f,r,f,f,f,f,f,f,r,f,f $_.TEMP.avinput $settings{annovardb}";
  system($command);
}

### STEP 2. Add additional info with custom script
print "INFO\tAdditional annotations: gdi, rvis, clinvar...\n";

open(INTOL, $settings{rvis});
while ($row=<INTOL>) {
	chomp($row);
	@line=split("\t",$row);
	$rvis{$line[0]} = $line[1];
}
close(INTOL);

open(CLINVAR, $settings{clinvar});
while ($row=<CLINVAR>) {
	chomp($row);
	@line=split("\t",$row);
	$clinvar{$line[1]} = $line[3];
}
close(CLINVAR);

open(GDI, $settings{gdi});
while ($row=<GDI>) {
	chomp($row);
	@line=split("\t",$row);
	$gdi{$line[0]} = $line[2];
}
close(GDI);

if (exists $settings{vardb}) {
  open(VARDB, $settings{vardb});
  while ($row=<VARDB>) {
	  chomp($row);
	  @line=split("\t",$row);
	  $vardb{$line[0]} = $line[2];
  }
  close(VARDB);
} else {
  print "INFO\tNo variant db file provided\n";
}

foreach (@outfiles) {
  open(IN, "$_.TEMP.avinput.hg19_multianno.txt");
  open(OUT, ">>$_.Annotated.txt");
  
  $headline=<IN>;
  chomp($headline);
  @colnames = split("\t", $headline);
  $c=0;
  foreach (@colnames) {
  	$cols{$_}=$c;
  	$c++;
  }
  pop(@colnames);
  $headline = join("\t", @colnames);
  foreach (@samplesid) {
    $headline .= "\tGT_$_\tGQ_$_\tDP_$_\tAD_$_\tGAR_$_";
  }
  $headline .= "\tClinVar\tRVIS\tGDI\tVarDB";
  
  print OUT $headline."\n";
  
  while($row=<IN>) {
    chomp($row);
    @line=split("\t", $row);
    @vcfline=@line[$cols{Otherinfo}..$#line];

    print "########################\n";
    print join("\t", @vcfline)."\n##\n";
    print "$vcfline[0]\n$vcfline[1]\n";
    
    $newfields = "";
    $varid=$vcfline[0]."_".$vcfline[1]."_".$vcfline[3]."_".$vcfline[4];
    
    @format=split(":", $vcfline[8]);
		$i=9;
		foreach $mysample(@samplesid) {
		  if ($vcfline[7] =~ /${mysample}_true=([0-9.E\-]+)/) {$info{$mysample}{GARFIELD} = $1} else {$info{$mysample}{GARFIELD} = "-99"}
		  @column=split(":", $vcfline[$i]);
		  for ($n=0; $n<=$#format; $n++) {
		    $info{$mysample}{$format[$n]}=$column[$n];
		  }
		  $i++;
		}
    
    foreach (@samplesid) {
      $newfields .= "$info{$_}{GT}\t$info{$_}{GQ}\t$info{$_}{DP}\t$info{$_}{AD}\t$info{$_}{GARFIELD}\t";
    }
    
    $end = $cols{Otherinfo} - 1;
    
    if (exists $settings{vardb}) {
      print OUT join("\t", @line[0..$end])."\t".$newfields.$clinvar{$line[6]}."\t".$intol{$line[6]}."\t".$gdi{$line[6]}."\t".$vardb{$varid}."\n";
    } else {
      print OUT join("\t", @line[0..$end])."\t".$newfields."\t".$clinvar{$line[6]}."\t".$intol{$line[6]}."\t".$gdi{$line[6]}."\n";
    }
  }
  close(IN);
  close(OUT);
}

########################
###  Cleaning steps  ###
########################

#Remove temp files
if ($notemp) {
  print "INFO\tRemoving temp files...\n";
  system("rm $outpath/*.TEMP.*");
}

$end_run=time();
$run_time= $end_run-$start_run;
print "\nAll Done!\nAnalysis time: $run_time seconds\n";
