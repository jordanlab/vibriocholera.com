#!/usr/bin/perl -w
# Aroon Chande
# Extracted HMM profiled proteins from a proteins.faa file
# Input os output of  hmmscan --tblout
# ./protein_extractor.pl -hmm hmmdir -prot proteindir -o outdir
use strict;
use Getopt::Long;
use File::Basename;
use File::Temp;
my ($fasta,$gff,$outdir,$hmmFile,$protFile,$outfile,$prots,@prots,%prots,@hits,$predict);
$prots = '';
if (@ARGV < 1){print_usage();exit 1;}
GetOptions ('fasta=s' => \$fasta, 'gff=s' => \$gff, 'predict=s' => \$predict);
$outdir = dirname($fasta);
#system(`mkdir -p $outdir/tmp`) or die "Cannot create directory: $outdir/tmp\n";
#if ($predict eq 'yes'){
#	system(`prodigal -q -c -i $fasta -a $outdir/prots.faa -f gff -o $outdir/prots.gff 2> /dev/null` );
#	system(`sed -i 's/\\s\#.*.//g' $outdir/prots.faa`);
#}
#else{
#	system(`fastaFromBed -s -fi $fasta -bed $gff -fo $outdir/prots.fna`);
#	system(`sed -i 's/\\://g' $outdir/prots.fna`);
#	system(`transeq -sequence $outdir/prots.fna -outseq $outdir/prots.faa 2> /dev/null`);
#	system(`cp $gff $outdir/prots.gff`);
#}
my @hmms = glob ( "/home/blast/prediction_server/hmm_profiles/*.hmm" );
foreach my $hmm (@hmms){
	my $base = basename($hmm);
	$base =~ s/\.hmm$//g;
	my $tblout = join('.',$base,"out","txt");
	my $out = temp_filename();
	system(`hmmscan --cpu 4 -o $out --tblout $outdir/$tblout $hmm $outdir/prots.faa`) or die "Cannot parse an input option: $!\n";
		
}
@hits = ( );
$prots = '';
%prots = (	);
my $predictions=temp_filename();
opendir DIR, $outdir or die "cannot open dir $outdir: $!";
my @file= readdir DIR;
closedir DIR;
foreach (@file){
if ($_ =~ /out.txt$/){
	$hmmFile = join('/',$outdir,$_);
	$_ =~ s/\.txt//g;
	my ($strain) = $_ =~ m/([a-zA-Z]*[0-9]*[^.])/;
	#Load protein file
	open PROT, "$outdir/prots.faa" or die "Cannot open $outdir/prots.faa: $!\n";
	foreach (<PROT>){
		$prots .= $_;
	}
	close PROT;
	@prots = split(/\>/,$prots);
	shift @prots;
	foreach (@prots){
		my($desc, $seq) = split(/\r?\n/,$_,2);
	    $seq =~ s/\*[^ACDEFGHIJKLMNPQRSTVWY]//g;
	    $prots{$desc} = $seq;
	}
	open HMM, $hmmFile or die "Cannot open $hmmFile: $!\n";
	my $i = 0;
	foreach (<HMM>){
		if(($_ =~ /\#/)){next;}
	    else{
			my @columns = split(/\s+/,$_);
			$hits[$i] = $columns[2];
			$i++;
		}
	}
	close HMM;
	open OUT, ">>$predictions" or die "Cannot open $outdir/predictions.faa : $!";
	foreach my $hit (@hits){
		my $rpstemp=temp_filename();
		my $rpsresulttemp=temp_filename();
		open RPS, ">$rpstemp" or die "Cannot open $rpstemp: $!";
		print RPS ">$hit\n$prots{$hit}\n";
		print "$hit\t$rpsresulttemp\t$rpstemp\n";
		close RPS;
		`rpsblast -num_threads 4 -db /home/blast/prediction_server/cdd/Cdd -max_target_seqs 1 -outfmt "6" -out $rpsresulttemp -query $rpstemp`;
		open RPSR, "$rpsresulttemp" or die "Cannot open $rpsresulttemp: $!";
		foreach (<RPSR>){
			my @columns = split(/\s+/,$_);
            if (($columns[1] =~ /250847|254315|226200/) & ($strain eq "lipase")){
            	print OUT ">$hit|$strain\n$prots{$hit}\n";
            }
            elsif (($columns[1] =~ /255682|226199|178040|182849|177427/) & ($strain eq "hydrolase")){
            	print OUT ">$hit|$strain\n$prots{$hit}\n";
			}
            elsif (($columns[1] =~ /212030|197609/) & ($strain eq "lysm")){
            	print OUT ">$hit|$strain\n$prots{$hit}\n";
			}
            elsif (($columns[1] =~ /260130|227061|224093|223187|273628|224228|184001/) & (($strain eq "transferase") || ($strain eq "ntpase"))){
            	print OUT ">$hit|$strain\n$prots{$hit}\n";
			}
			elsif (($columns[1] =~ /274542|224482/) & ($strain eq "vgrg")){
            	print OUT ">$hit|$strain\n$prots{$hit}\n";
			}
            elsif ($strain eq "unknown"){
                print OUT ">$hit|$strain\n$prots{$hit}\n";
			}
		}
		close RPSR;
		
	}
}	else {next;}
}
close OUT;
$prots = '';
%prots = (	);
open PROT, $predictions or die "Cannot open $outdir/prots.faa: $!\n";
foreach (<PROT>){
	$prots .= $_;
}
close PROT;
@prots = split(/\>/,$prots);
shift @prots;
foreach (@prots){
	my($desc, $seq) = split(/\r?\n/,$_,2);
    $seq =~ s/\*[^ACDEFGHIJKLMNPQRSTVWY]//g;
    $prots{$desc} = $seq;
}
open OUT, ">$outdir/predictions.faa" or die "Cannot open $outdir/predictions.faa: $!";
foreach (sort keys %prots){
	print OUT ">$_\n$prots{$_}";
}
close OUT;
exit 0;
###############################################################################
sub temp_filename{
	my $file = File::Temp->new(
	TEMPLATE => 'tempXXXXX',
	DIR      => '/tmp/',
	);
} 

