#!/usr/bin/perl
##Used to compare VCF files from Proton to NIST calls
# -i input .txt with list of files
# -o output file
# -f threshold of fraction in samples

use Getopt::Std;

#Parsing command line options
getopt('i:o:f:', \%opts);

#Store 1000F AF
#print "Reading 1000G AF...\n";
#open(KG, "1000G_indels_AF.txt");
#while ($row=<KG>) {
#	chomp($row);
#	@line=split("\t", $row);
#	$id=$line[0]."_".$line[1];
#	$freqKG{$id}=$line[4];
#}
#close(KG);

#Store ExAC AF
#print "Reading ExAc AF...\n";
#open(EXAC, "ExAC_0.2_indels_AF.txt");
#while ($row=<EXAC>) {
#	chomp($row);
#	@line=split("\t", $row);
#	$id=$line[0]."_".$line[1];
#	$freqEXAC{$id}=$line[4];
#}
#close(EXAC);

#open file with list of files
open(LIST, $opts{i});
while ($row=<LIST>) {
	chomp($row);
	push(@files, $row);
}
close(LIST);

#Analyze over files
foreach $myfile(@files) {

print "Reading $myfile...\n";
open(IN, $myfile);
while(<IN> =~ /^##/) {}
while($row=<IN>) {
	chomp($row);
	@line=split("\t", $row);
	$idperfect=$line[0]."_".$line[1]."_".$line[3]."_".$line[4];
	$idpos=$line[0]."_".$line[1];
	$perfect{$idperfect}{pos}=$idpos;
	$perfect{$idperfect}{count}++;
	$position{$idpos}++;
}
close(IN);

}

open(OUT, ">>$opts{o}.perfect.txt");
print OUT "Total files: ".($#files+1)."\n";
foreach $mykey(keys %perfect) {
	$fraction=$perfect{$mykey}{count}/($#files+1);
	if ($fraction > 1) {$fraction = 1}
	if ($fraction >= $opts{f}) {print OUT "$mykey\t$perfect{$mykey}{count}\t$fraction\n"}
}
close(OUT);

open(OUT, ">>$opts{o}.position.txt");
print OUT "Total files: ".($#files+1)."\n";
foreach $mykey(keys %position) {
	$fraction=$position{$mykey}{count}/($#files+1);
	if ($fraction > 1) {$fraction = 1}
	if ($fraction >= $opts{f}) {print OUT "$mykey\t$position{$mykey}\t$fraction\n"};
}
close(OUT);
