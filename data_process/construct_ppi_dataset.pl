use List::Util qw/sum/;
use POSIX;
use Statistics::R;
open PPI, "PPI/ppi_source/hippie_current.txt";
my %syn;
my %sy;
while(<PPI>){
	chomp;
	my @s=split/\t/;
	foreach my $i(split /,/,$s[0]){
		$syn{$i}=$s[1];
		$sy{$s[1]}{$i}=1;
	}
	foreach my $i(split /,/,$s[2]){
		$syn{$i}=$s[3];
		$sy{$s[3]}{$i}=1;
	}
	
}
close PPI;
open FASTA, "PPI/uniprot_protein/uniprot_sprot.fasta";
my $seq;
my $seqid;
my %seqE;
while(<FASTA>){
	chomp;
	if(/^>sp\|.+\|(.+?)\s/){
		my $len=length($seq);
		if($len>0){
			if(exists($syn{$seqid})){
				my @ids = keys %{$sy{$syn{$seqid}}};
				foreach my $y(@ids){
					$seqE{$y}{"len"}=$len;
					$seqE{$y}{"seq"}=$seq;
				}
			}else{
				$seqE{$seqid}{"len"}=$len;
				$seqE{$seqid}{"seq"}=$seq;
			}
		}
		$seqid=$1;
		$seq="";
	}
	else{
		$seq .= $_;
	}
}

open PPI, "PPI/ppi_source/hippie_current.txt";
my %protein;
my %ppi;
my $pnum=0;
my %pi;
open TR,">dataset/positive_set.csv";
open NR,">dataset/negative_set.csv";
my @pos;
while(<PPI>){
	chomp;
	my @s=split/\t/;
	$s[0] =~ s/^(.+?),.+$/$1/;
	$s[2] =~ s/^(.+?),.+$/$1/;
	unless(exists($seqE{$s[0]}{"seq"}) & exists($seqE{$s[2]}{"seq"})){
		next;
	}
	
	if($s[5]=~/in vivo/ | $s[5]=~/in vitro/){
		if (exists($seqE{$s[0]}{"seq"}) & exists($seqE{$s[2]}{"seq"})){
			my $feature="$s[0],$s[2]";
			printf TR $feature."\n";
			$ppi{"$s[0]:$s[2]"}=1;$ppi{"$s[2]:$s[0]"}=1;
			$protein{$s[0]}=1;$protein{$s[2]}=1;
			$pos[$pnum]=$feature;
			$pnum++;
		}
	}
	$pi{"$s[0]:$s[2]"}=1;$pi{"$s[2]:$s[0]"}=1;
}
close PPI;

my %test_ppi;
my @test_dataset = ("Venkatesan09","Stelzl05","Lehner04","Kaltenbach07","Bell09","Rual05","Goehler04","Colland04","Albers05","Nakayama02");
foreach my $ds(@test_dataset){
	open PPI, "PPI/ppi_source/hippie_current.txt";
	open TT, ">dataset/$ds"."_pos_test.csv";
	while(<PPI>){
		chomp;
		my @s=split/\t/;
		unless(exists($seqE{$s[0]}{"seq"}) & exists($seqE{$s[2]}{"seq"})){
			next;
		}
		if($s[5]=~/$ds/ & !exists($ppi{"$s[0]:$s[2]"})){
			if(exists($seqE{$s[0]}{"seq"}) & exists($seqE{$s[2]}{"seq"})){
				my $feature="$s[0],$s[2]";
				print TT $feature."\n";
				$test_ppi{$ds} += 1;
			}
		}
	}
	close TT;
	close PPI;
}



my @proteins=keys %protein;
my $np=@proteins;
my %nppi;
my @neg;
my $nnum=0;
while($nnum<200*$pnum){
	my $pA=$proteins[int(rand($np))];
	my $pB=$proteins[int(rand($np))];
	unless(exists($seqE{$pA}{"seq"}) & exists($seqE{$pB}{"seq"})){
		next;
	}
	if(!exists($pi{"$pA:$pB"}) & !exists($nppi{"$pA:$pB"})){
		if(exists($seqE{$pA}{"seq"}) & exists($seqE{$pB}{"seq"})){
			my $feature="$pA,$pB";
			$nppi{"$pA:$pB"}=1;$nppi{"$pB:$pA"}=1;
			print NR $feature."\n";
			$neg[$nnum]=$feature;
			$nnum++;
		}
	}
}
close TR;close NR;

sub ran{
	my $max = shift;
	my $n = shift;
	my @ran;
	for(my $i=0;$i<$n;$i++){
		my $num=int(rand($max));
		$ran[$i]=$num;
	}
	return \@ran;
}

sub urran{
	my $max = shift;
	my $n = shift;
	my %sns;
	for(my $i=0;$i<$n;$i++){
	    do {
	       $num=int(rand($max));
	    } while (exists $sns{$num});
	    $sns{$num} = 1;
	}
	my @ran=keys %sns;
	return \@ran;
}
foreach my $ds(@test_dataset){
	open TT, ">dataset/$ds"."_neg_test.csv";
	my @ni = @{urran($nnum,$test_ppi{$ds})};
	my @rn = @neg[@ni];
	foreach my $j(@rn){
		print TT "$j\n";
	}
	close TT;
}
for(my $i=0; $i<100; $i++){
	open F, ">dataset/train_pos_set_$i.csv";
	open F1, ">dataset/train_neg_set_$i.csv";
	my @p = @{ran($pnum,$pnum)};
	
	my @n = @{ran($nnum,$nnum)};
	my @pi = @{urran($pnum,800)};
	my @ni = @{urran($nnum,800)};
	
	my @rp = @pos[@p[@pi]];
	#print "@rp\n";
	my @rn = @neg[@n[@ni]];
	#print "@rn\n";
	foreach my $i(@rp){
		print F "$i\n";
	}
	foreach my $j(@rn){
		print F1 "$j\n";
	}
	close F;close F1;
}
