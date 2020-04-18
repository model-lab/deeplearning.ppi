use List::Util qw/sum/;
use POSIX;
use Statistics::R;
open AI, "PPI/AAindex/aaindex1";

my %index;my $id; my $value;
my $class =7;

while(<AI>){
	chomp;
	if(/^H\s(.+)/){
		
		if(length($value)>0){
			$value =~ s/^\s+|\s+$//;
			my @d=split /\s+/,$value;
			my @cent=random(\@d,$class);
			unless(@cent<$class){
				#print "$id\n@d\n@cent\n";
				my $R = Statistics::R->new();
				$R->set("cent",\@cent);
				$R->set("x",\@d);
				$R->set('i',$class);
				$R->run(q`y<-kmeans(x,centers=cent,iter.max = 10)`);
				my @d=@{$R->get('y$cluster')};
				#print "@d\n";
				$R->stop;
				$index{$id}{"A"}=$d[0];$index{$id}{"R"}=$d[1];$index{$id}{"N"}=$d[2];$index{$id}{"D"}=$d[3];
				$index{$id}{"C"}=$d[4];$index{$id}{"Q"}=$d[5];$index{$id}{"E"}=$d[6];$index{$id}{"V"}=$d[19];
				$index{$id}{"G"}=$d[7];$index{$id}{"H"}=$d[8];$index{$id}{"I"}=$d[9];$index{$id}{"L"}=$d[10];
				$index{$id}{"K"}=$d[11];$index{$id}{"M"}=$d[12];$index{$id}{"F"}=$d[13];$index{$id}{"P"}=$d[14];
				$index{$id}{"S"}=$d[15];$index{$id}{"T"}=$d[16];$index{$id}{"W"}=$d[17];$index{$id}{"Y"}=$d[18];
			}
		}
		$id=$1;
		$value="";
	}
	else{
		if(/^I\s/){
			my $s1=<AI>;my $s2=<AI>;chomp($s1);chomp($s2);
			$value =$s1.$s2;
		}
	}
}

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

my $inter=3;# 3:ATC or 4:ATCG
my %dict_motif;
my @element=("1","2","3","4","5","6","7");
for(my $i=0;$i<$inter;$i++){
    if($i==0){
        foreach my $base(@element[0..$class-1]){
            $dict_motif{$base}=0;
        }
    }
    else {
        foreach my $old(keys %dict_motif){
            foreach my $base(@element[0..$class-1]){
        		my $new=$old.$base;
        		$dict_motif{$new}=0;
            }
            delete $dict_motif{$old};
        }
    }
}
my @aaindex=("JANJ780101","MIYS990102","NAGK730102","NAKH920106","RICJ880109","OOBM770104");
my @motif = keys %dict_motif;
#my @aaindexs=keys %index;
#my @aaindex = @aaindexs[0..5];
for my $ii(0 .. 99){
    open PS, "dataset/train_pos_set_$ii.csv";
    open NS, "dataset/train_neg_set_$ii.csv";
    open TR,">encoding_dataset/train_dataset_$ii.csv";
    while(<PS>){
        chomp;
		my @s=split/,/;
        my $feature="1";
		my @sq1=split //,$seqE{$s[0]}{'seq'};
		my @sq2=split //,$seqE{$s[1]}{'seq'};
		
		foreach my $i(@aaindex){
			my $sequence="";
			foreach my $a(@sq1){
				$sequence .= $index{$i}{$a};
				#print "$a\n$index{$i}{$a}\n";
			}
			
			my %first_motif;
			my $max=0;my $min=100000;
			foreach my $m(@motif){
				$first_motif{$m} =()=($sequence =~ /$m/g);
				if($first_motif{$m}>$max){
					$max=$first_motif{$m};
				}
				if($first_motif{$m}<$min){
					$min=$first_motif{$m};
				}
				#print "$m\t$first_motif{$m}\n";
			}
			
			foreach my $m(@motif){
				$first_motif{$m} =($first_motif{$m}-$min)/($max-$min);
			}
			
			$sequence="";
			foreach my $a(@sq2){
				$sequence .= $index{$i}{$a};
			}
			
			my %second_motif;
			my $max=0;my $min=100000;
			foreach my $m(@motif){
				$second_motif{$m} =()=($sequence =~ /$m/g);
				if($second_motif{$m}>$max){
					$max=$second_motif{$m};
				}
				if($second_motif{$m}<$min){
					$min=$second_motif{$m};
				}
			}
			
			foreach my $m(@motif){
				$second_motif{$m}=($second_motif{$m}-$min)/($max-$min);
				my $freq = sprintf("%.4f",($first_motif{$m}+$second_motif{$m})**2);
				$feature .= ",$freq";
			}
		}
		printf TR $feature."\n";
    }
    close PS;
    while(<NS>){
        chomp;
		my @s=split/,/;
        my $feature="0";
		my @sq1=split //,$seqE{$s[0]}{'seq'};
		my @sq2=split //,$seqE{$s[1]}{'seq'};
		foreach my $i(@aaindex){
			my $sequence;
			foreach my $a(@sq1){
				$sequence .= $index{$i}{$a};
			}
			my %first_motif;
			my $max=0;my $min=100000;
			foreach my $m(@motif){
				$first_motif{$m} =()=($sequence =~ /$m/g);
				if($first_motif{$m}>$max){
					$max=$first_motif{$m};
				}
				if($first_motif{$m}<$min){
					$min=$first_motif{$m};
				}
			}
			
			foreach my $m(@motif){
				$first_motif{$m} =($first_motif{$m}-$min)/($max-$min);
			}
			
			$sequence="";
			foreach my $a(@sq2){
				$sequence .= $index{$i}{$a};
			}
			
			my %second_motif;
			my $max=0;my $min=100000;
			foreach my $m(@motif){
				$second_motif{$m} =()=($sequence =~ /$m/g);
				if($second_motif{$m}>$max){
					$max=$second_motif{$m};
				}
				if($second_motif{$m}<$min){
					$min=$second_motif{$m};
				}
			}
			
			foreach my $m(@motif){
				$second_motif{$m}=($second_motif{$m}-$min)/($max-$min);
				my $freq = sprintf("%.4f",($first_motif{$m}+$second_motif{$m})**2);
				$feature .= ",$freq";
			}
		}
		print TR $feature."\n";
    }
    close TR;
    close NS;
}
my @dataset = ("Venkatesan09","Stelzl05","Lehner04","Kaltenbach07","Bell09","Rual05","Goehler04","Colland04","Albers05","Nakayama02","Lim06");
foreach my $ii(@dataset){
    open PS, "dataset/$ii"."_pos_test.csv";
    open NS, "dataset/$ii"."_neg_test.csv";
    open TR,">encoding_dataset/$ii"."_test_dataset.csv";
    while(<PS>){
        chomp;
		my @s=split/,/;
        my $feature="1";
		my @sq1=split //,$seqE{$s[0]}{'seq'};
		my @sq2=split //,$seqE{$s[1]}{'seq'};
		foreach my $i(@aaindex){
			my $sequence="";
			foreach my $a(@sq1){
				$sequence .= $index{$i}{$a};
				#print "$a\n$index{$i}{$a}\n";
			}
			
			my %first_motif;
			my $max=0;my $min=100000;
			foreach my $m(@motif){
				$first_motif{$m} =()=($sequence =~ /$m/g);
				if($first_motif{$m}>$max){
					$max=$first_motif{$m};
				}
				if($first_motif{$m}<$min){
					$min=$first_motif{$m};
				}
				#print "$m\t$first_motif{$m}\n";
			}
			#print "$seqE{$s[0]}{'seq'}\n@sq1\n$sequence\n$min\t$max\n";
			foreach my $m(@motif){
				$first_motif{$m} =($first_motif{$m}-$min)/($max-$min);
			}
			
			$sequence="";
			foreach my $a(@sq2){
				$sequence .= $index{$i}{$a};
			}
			
			my %second_motif;
			my $max=0;my $min=100000;
			foreach my $m(@motif){
				$second_motif{$m} =()=($sequence =~ /$m/g);
				if($second_motif{$m}>$max){
					$max=$second_motif{$m};
				}
				if($second_motif{$m}<$min){
					$min=$second_motif{$m};
				}
			}
			
			foreach my $m(@motif){
				$second_motif{$m}=($second_motif{$m}-$min)/($max-$min);
				my $freq = sprintf("%.4f",($first_motif{$m}+$second_motif{$m})**2);
				$feature .= ",$freq";
			}
		}
		printf TR $feature."\n";
    }
    close PS;
    while(<NS>){
        chomp;
		my @s=split/,/;
        my $feature="0";
		my @sq1=split //,$seqE{$s[0]}{'seq'};
		my @sq2=split //,$seqE{$s[1]}{'seq'};
		foreach my $i(@aaindex){
			my $sequence;
			foreach my $a(@sq1){
				$sequence .= $index{$i}{$a};
			}
			my %first_motif;
			my $max=0;my $min=100000;
			foreach my $m(@motif){
				$first_motif{$m} =()=($sequence =~ /$m/g);
				if($first_motif{$m}>$max){
					$max=$first_motif{$m};
				}
				if($first_motif{$m}<$min){
					$min=$first_motif{$m};
				}
			}
			
			foreach my $m(@motif){
				$first_motif{$m} =($first_motif{$m}-$min)/($max-$min);
			}
			
			$sequence="";
			foreach my $a(@sq2){
				$sequence .= $index{$i}{$a};
			}
			
			my %second_motif;
			my $max=0;my $min=100000;
			foreach my $m(@motif){
				$second_motif{$m} =()=($sequence =~ /$m/g);
				if($second_motif{$m}>$max){
					$max=$second_motif{$m};
				}
				if($second_motif{$m}<$min){
					$min=$second_motif{$m};
				}
			}
			
			foreach my $m(@motif){
				$second_motif{$m}=($second_motif{$m}-$min)/($max-$min);
				my $freq = sprintf("%.4f",($first_motif{$m}+$second_motif{$m})**2);
				$feature .= ",$freq";
			}
		}
		print TR $feature."\n";
    }
    close TR;
    close NS;
}

sub random{
	my $g = shift;
	my $cls = shift;
	my %lis;
	my @centers;
	foreach my $l(@{$g}){
		$lis{$l}=1;
	}
	unless(exists($lis{'NA'})){
	    my @list = sort{$a<=>$b} @{$g};
	    my $count = @list;
		my %cent;
		for (my $i=1;$i<=$cls;$i++){
			my $a = ($count-1)*$i/$cls;
			$cent{$list[int($a)]}=1;
		}
		my @center=keys %cent;
		@centers= sort{$a<=>$b} @center;
	}
	return @centers;
}