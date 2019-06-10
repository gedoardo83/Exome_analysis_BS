#Processing variants in a VCF file and find segregating variants according to AR or AD model
#Need multisample VCF and a status file
use Getopt::Long;

############################
###  Processing options  ###
############################

#Set default values for options
$model= "AR" ;
$outpath= ".\\" ;
$maf=0.01;

#Parsing command line options
GetOptions(
	'input=s' => \$inputvcf, #Multisample vcf file for the trio / quad
	'status=s' => \$pedfile, #Status file (tab-delimited file with samplesID as in VCF file and 1=affected 0=unaffected)
	'maf:f' => \$maf, #MAF threshold for filtering
	'output:s' => \$outpath, #Output folder
	'model:s' => \$model, #Possible values AD / AR
) or die "Uncorrect or missing mandatory arguments\n";

#################################
###  Read sample status file  ###
#################################

print "INFO\tReading status file...\n";
open(PED, $pedfile);
while($row=<PED>) {
	chomp($row);
	@line = split("\t", $row);
	if ($line[1]==0) {$control{$line[0]}=0}
	if ($line[1]==1) {$affected{$line[0]}=0}
}
close(PED);

############################################
###  Read variants from multisample VCF  ###
############################################

open (VCF, "$inputvcf");
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

######################################
###  Variant segregation analysis  ###
######################################

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
