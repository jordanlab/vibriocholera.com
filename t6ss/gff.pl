#!/usr/bin/perl -w
# Aroon Chande
use strict;
use Getopt::Long;
use File::Basename;
use File::Temp;
my ($fasta,$gff,$outdir,$hmmFile,$protFile,$outfile,$prots,@prots,%prots,@gff,%gff,@hits,$predict);
$prots = '';
GetOptions ('fasta=s' => \$fasta, 'gff=s' => \$gff, 'predict=s' => \$predict);
$outdir = dirname($fasta);
open PROT, "$outdir/predictions.faa" or die "Cannot open $outdir/predictions.faa $!\n";
foreach (<PROT>){
	$prots .= $_;
}
close PROT;
@prots = split(/\>/,$prots);
shift @prots;
foreach (@prots){
	my($desc, undef) = split(/\r?\n/,$_,2);
	my @desc = split(/\|/,$desc);
    $prots{$desc[0]} = $desc[1];
}
if (!(defined $gff)){$gff = "$outdir/prots.gff";}
open GFF, $gff or die "Cannot open $outdir/prots.gff $!\n";
foreach(<GFF>){
	my @cols = split(/\t/,$_);
	if ($cols[0] =~ /\#/){next;}
	else{
		#print "$cols[0]\t$cols[8]\n";
		my $base = $cols[0];
		my $id = $cols[8];
		my @idcols = split(/;/,$id);
		($id) = $idcols[0] =~ m/(_\d*)/;
		$id = join('',$cols[0],$id);
		$gff{$id}{start} = $cols[3];
		$gff{$id}{stop} = $cols[4];
		$gff{$id}{strand} = $cols[6];
		if(defined $prots{$id}){$gff{$id}{gene} = $prots{$id};}	
		else{$gff{$id}{gene} = "-";}
	}
}	
#foreach (keys %gff){	print "$_\t$gff{$_}{start}\t$gff{$_}{stop}\t$gff{$_}{strand}\t$gff{$_}{gene}\n"	}
no warnings 'uninitialized';
my %files;
my $j = 1;
my $scale;
my $max;
foreach my $key(keys %prots){
	if ($prots{$key} =~ /vgrg/){
		my ($num) = $key =~ m/(\d*$)/;
		$key =~ s/_\d*$//g;
		open OUT, ">$outdir/$j.ptt";
		print OUT "Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n";
		for(my $i = -5; $i <6; $i++){
			if($num + $i < 1){next;}
			else{
				my $gffkey = join("_",$key,$num+$i);
				if ((defined $gff{$gffkey}) & !(defined $scale)){$scale = $gff{$gffkey}{start}-1;}
				$max = $gff{$gffkey}{stop}-$scale if (defined $gff{$gffkey});
				print OUT $gff{$gffkey}{start}-$scale,"..",$gff{$gffkey}{stop}-$scale,"\t$gff{$gffkey}{strand}\t",abs($gff{$gffkey}{start}-$gff{$gffkey}{stop}),"\t-\t$gff{$gffkey}{gene}\t-\t-\t-\t-\n" if (defined $gff{$gffkey});
			}
		}
		close OUT;
		$files{$j}{max}=$max;
		undef $scale;
		$j++;
	}
}
foreach my $key (keys %files){
	print "$key\n";
	open IN, "$outdir/$key.ptt";
	my @temp = <IN>;
	close IN;
	open OUT, ">$outdir/$key.ptt";
	print OUT "Genome - 1..$files{$key}{max}\n",@temp-1," proteins\n";
	print OUT @temp;
	close OUT;	
}

