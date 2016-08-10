9a10
> use Parallel::ForkManager;
25a27
> my $manager = Parallel::ForkManager -> new ( 2 );
27a30
> 	$manager->start and next;
32,33c35,36
< 	system(`hmmscan --cpu 4 -o $out --tblout $outdir/$tblout $hmm $outdir/prots.faa`) or die "Cannot parse an input option: $!\n";
<
---
> 	system(`hmmscan --cpu 8 -o $out --tblout $outdir/$tblout $hmm $outdir/prots.faa`) or die "Cannot parse an input option: $!\n";
> 	$manager->finish;
34a38
> $manager->wait_all_children;
41a46
> $manager = Parallel::ForkManager -> new ( 4 );
72a78
> 		$manager->start and next;
79c85
< 		`rpsblast -num_threads 4 -db /home/blast/prediction_server/cdd/Cdd -max_target_seqs 1 -outfmt "6" -out $rpsresulttemp -query $rpstemp`;
---
> 		`rpsblast -num_threads 1 -db /home/blast/prediction_server/cdd/Cdd -max_target_seqs 1 -outfmt "6" -out $rpsresulttemp -query $rpstemp`;
103c109
<
---
> 		$manager->finish;
106a113
> $manager->wait_all_children;

