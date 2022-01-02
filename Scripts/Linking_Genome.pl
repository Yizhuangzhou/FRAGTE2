#!usr/bin/perl -w
use strict;
use POSIX qw(strftime);
use Time::Local;
use Cwd qw(abs_path getcwd);
use Statistics::R;
die "perl $0 [GenomeInfo][outdir][output][len]" unless @ARGV==4;
my $t1=timelocal(localtime());

############################## Intializing for Zvalue ##############################################
my $km = 4;
my @mono = ("A","T","G","C");
my @nuc = @mono;
my @oligo_kmer = (); 
my $index = 1;
while(){
	my @tmpnuc = ();
	$index++;
	foreach(@mono){
		my $word = $_;
		foreach (@nuc){
			my $nuc = $word."$_";
			push @tmpnuc,$nuc;
		}
	}
	@nuc = @tmpnuc;
	if($index == $km){
		@oligo_kmer = @nuc;
		last;
	}
}
my @UsedKmer=();
my %tag=();
foreach my $k1 (@oligo_kmer) {
	my $k2=reverse($k1);
	$k2=~tr/ATGC/TACG/;
	if($tag{$k2}){
		next;
	}
	else{
		push @UsedKmer,$k1;
		$tag{$k1}=1;
	}
}

############################## Main ############################################################
open(IF,$ARGV[0])||die;
my $dir=$ARGV[1];
$dir = abs_path($dir);
mkdir $dir,0755 unless (-e $dir);
open(OUT,">$ARGV[2]")||die;
my $FragLen = $ARGV[3];
while(<IF>){
	chomp;
	my @a=split /\t/;
	next if($a[5]<10000);
	my $file=$a[6];
	open(FA,$file)||die;
	my $MAS="";
	my $SeqID="";
	my $Size=0;
	my $seq="";
	my %MAS=();
	my %MASlen=();
	my @MASlen=();
	while(<FA>){
		chomp;
		if(/^>(\S+)/){
			if($SeqID ne ""){
				my $MASlen=length $MAS;
				if($MASlen>=1000){
					$Size+=$MASlen;
					push @MASlen,$MASlen;
					$seq.=$MAS;
					$MAS{$SeqID}=$MAS;
					$MASlen{$SeqID}=$MASlen;
				}
			}
			$SeqID=$1;
			$MAS="";
		}
		else{
			s/[^ATGC]//gi;
			$MAS.=uc($_);
		}
	}
	close FA;
	my $MASlen=length $MAS;
	if($MASlen >=1000){
		$Size+=$MASlen;
		push @MASlen,$MASlen;
		$seq.=$MAS;
		$MAS{$SeqID}=$MAS;
		$MASlen{$SeqID}=$MASlen;
	}
	next if($Size<10000);
	if(scalar keys %MAS >1){
		my $LinkedSeq=&Bayesian(\%MAS,$seq,$Size,\%MASlen);
		$seq=$LinkedSeq;
	}
	my $MeanLen=&mean(\@MASlen);
	open(OF,">$dir/$a[0]\.fna")||die;
	print OF ">$a[0]\n$seq\n";
	close OF;
	print OUT "$a[0]\t$a[1]\t$a[2]\t$MeanLen\t$a[4]\t$Size\t$dir/$a[0]\.fna\n";
	my $current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
	print STDERR "$current_time: $a[0]\n";
}
close IF;
close OUT;
my $t2=timelocal(localtime());
my $s=$t2-$t1;
my $hour=int($s / 3600);
my $ms=$s % 3600;
my $minutes=int($ms / 60);
my $second=$ms % 60;
print STDERR "Time elaspe in seconds for determining phase: $s\n";
print STDERR "Time elaspe for determining phase: $hour\:$minutes\:$second\n";

################################################################# Bayesian #####################################################
sub Bayesian{
	my ($SeqRef,$seq,$len,$LenRef)=@_;
	my %seq=%{$SeqRef};
	my %len=%{$LenRef};
	my $addlen=$FragLen;
	###TNF
	my $r=Statistics::R->new();
	$r->run(q'TNF<-c()');
	for(my $i=0;$i<=$len-$addlen;$i+=$addlen){
		my $fragment=substr($seq,$i,$addlen);
		my %TNF=();
		for(my $j=0;$j<=$addlen-$km;$j++){
			my $sub=substr($fragment,$j,$km);
			$TNF{$sub}++;
		}
		my @TNF=();
		foreach my $k1 (@UsedKmer){
			my $k2=reverse($k1);
			$k2=~tr/ATGC/TACG/;
			my $kmer=0;
			if($TNF{$k1}){
				$kmer+=$TNF{$k1};
			}
			if($k1 ne $k2){
				if($TNF{$k2}){
					$kmer+=$TNF{$k2};
				}
				if(! $kmer){
					$kmer=1;
				}
				else{
					$kmer=$kmer/2;
				}
			}
			else{
				if(! $kmer){
					$kmer=1;
				}
			}
			push @TNF,$kmer;
		}
		$r->set("x",\@TNF);
		$r->run(q'TNF<-append(TNF,x)');
	}

	$r->run(qq'mTNF<-matrix(TNF,ncol=136,byrow=T)
		Median_TNF<-apply(mTNF,2,median)');
	my $MedianTNFref=$r->get('Median_TNF');
	my @MedianTNF=@{$MedianTNFref};
	my $sum=0;
	foreach(@MedianTNF){
		$sum+=$_;
	}
	my %prob=();
	for(my $i=0;$i<=$#UsedKmer;$i++){
		my $tetra=$UsedKmer[$i];
		my $kmer=$MedianTNF[$i];
		my $val=$kmer/$sum;
		$prob{$tetra}=log($val)/log(10);
	}

	my %PP=();
	my $LinkedSeq="";
	foreach my $id (keys %seq) {
		my $seq=$seq{$id};
		my $len=$len{$id};
		my %TNF=();
		for(my $i=0;$i<=$len-$km;$i++){
			my $tetra=substr($seq,$i,$km);
			$TNF{$tetra}++;
		}
		my $PP=0;
		foreach my $tetra (@UsedKmer) {
			if($TNF{$tetra}){
				$PP +=$TNF{$tetra}*$prob{$tetra};
			}
		}
		$PP{$id}=$PP*500/$len;
	}
	my @SeqID=sort {$PP{$b} <=> $PP{$a}} keys %seq;
	foreach my $id (@SeqID) {
		$LinkedSeq.=$seq{$id};
	}
	return $LinkedSeq;
}

 ###################################################### Subroutine for computing mean##############################################
sub mean{
	my ($ref)=@_;
	my @a=@{$ref};
	my $sum=0;
	foreach(@a){
		$sum+=$_;
	}
	my $mean=int($sum/@a);
	return $mean;
}
