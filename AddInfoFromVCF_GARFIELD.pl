#!/usr/bin/perl

#Argument order: 
#1.Annotated variant file from Annovar 
#2.VCF file containing samples 
#3.File with intolerance score 
#4.Additional gene list from ClinVar (gene_condition_source_id)
#5.File with GDI score 
#6.List of recurrent variants in HiQ chemistry
#Output order: FILTER,GENOTYPE,GQ,DP,Known genes R/D,Intolerance

open(VCF, $ARGV[1]);

while ($row = <VCF>) {
	$part1="";
	$part2="";
	if ($row =~ /^##/) {
	} elsif ($row =~ /^#CHROM/) {
		chomp($row);
		@line=split("\t",$row);
		for ($h=9; $h<=$#line; $h++) {
		$samplesid[$h]=$line[$h];
		#print "$h\t$line[$h]\t$samplesid[$h]\n";		
		}			
	} else {
	
	chomp($row);
	@line=split("\t",$row);
	for ($i=9; $i<=$#line; $i++) {
		
		#Create Hash for values in Sample colum
		@format=split(":", $line[8]);
		@column=split(":", $line[$i]);
		
		for ($n=0; $n<=$#format; $n++) {
			$info{$samplesid[$i]}{$format[$n]}=$column[$n];
			#print "$format[$i]\t$column[$i]\n";
		}

		#Save GARFIELD prediction value fro each sample
		if ($line[7] =~ /$samplesid[$i]_true=([0-9.E\-]+)/) {$info{$samplesid[$i]}{GARFIELD} = $1}
		
		#Create Hash for values in INFO colum
		#@infotags=split(";", $line[7]);
		#foreach $tag(@infotags) {
		#$tag =~ /([^=]+)=(.+)/;
		#$info{$1}=$2;
		#print "$1\t$2\n";
		#}
	
	}
	
	$var{$line[0]."_".$line[1]} = "$line[5]\t$line[3]\t$line[4]";

	foreach $mykey (keys %info) {
		$part1 .= "\t$info{$mykey}{GT}";
		$part2 .= "\t$info{$mykey}{GQ}\t$info{$mykey}{DP}\t$info{$mykey}{GARFIELD}";
	}
	$var{$line[0]."_".$line[1]} .= $part1.$part2;	
	}

}

close(VCF);

open(INTOL, $ARGV[2]);
while ($row=<INTOL>) {
	chomp($row);
	@line=split("\t",$row);
	$intol{$line[0]} = $line[1];	
}
close(INTOL);

open(CLINVAR, $ARGV[3]);
while ($row=<CLINVAR>) {
	chomp($row);
	@line=split("\t",$row);
	$clinvar{$line[1]} = $line[3];	
}
close(CLINVAR);

open(GDI, $ARGV[4]);
while ($row=<GDI>) {
	chomp($row);
	@line=split("\t",$row);
	$gene{$line[0]} = $line[2];	
}
close(GDI);

open(HIQ, $ARGV[5]);
while ($row=<HIQ>) {
	chomp($row);
	@line=split("\t",$row);
	$hiq{$line[0]} = $line[2];	
}
close(HIQ);


open(ANNO, $ARGV[0]);
$head = <ANNO>;
chomp($head);

foreach $mykey (keys %info){
	$infohead1 .= "$mykey\t";
	$infohead2 .= $mykey."_GQ\t".$mykey."_DP\t".$mykey."_GARFIELD\t";
	$nosample .= "NA\tNA\tNA\tNA\t";
}
print $head."\tQUAL\tVCF_REF\tVCF_ALT\t".$infohead1.$infohead2."RVIS_score\tClinVar_Disease\tGDI_Score\tFreq_HiQ\n";

while ($row = <ANNO>) {
	chomp($row);
	@line=split("\t",$row);
	if ($line[4] eq "-") {
		$pos = $line[1]-1;
		$perfect = $line[0]."_".$pos."_".$line[3]."_".$line[4];
	} else {
		$pos=$line[1];
		$perfect = $line[0]."_".$pos."_".$line[3]."_".$line[4];
	}
	if (exists $var{$line[0]."_".$pos}) {
		print $row."\t".$var{$line[0]."_".$pos}."\t".$intol{$line[6]}."\t".$clinvar{$line[6]}."\t".$gene{$line[6]}."\t".$hiq{$perfect}."\n";
	} else {
		print $row."\tNA\tNA\tNA\t".$nosample.$intol{$line[6]}."\t".$clinvar{$line[6]}."\t".$gene{$line[6]}."\t".$hiq{$perfect}."\n";
	}
}
close(ANNO);
