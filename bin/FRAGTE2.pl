#!/usr/bin/perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Basename;
use POSIX qw(strftime);
use Time::Local;
use Cwd qw(abs_path getcwd);
use Getopt::Long;
use Statistics::R;
=head1 Name
 
 -------------------------------------------------------------------------------------------------------------
 FRAGTE2.pl  --  fragment tetranucleotide frequency correlation coefficient (FRAGTE) v2
 -------------------------------------------------------------------------------------------------------------

=head1 Description
 
 -------------------------------------------------------------------------------------------------------------
 The FRAGTE2 approach is an enhanced algorithm to sieve genome pairs for species delineation in prokaryotes
 --------------------------------------------------------------------------------------------------------------

=head1 Contact & Version
 
 ------------------------------------------------------------------------------------------------------------
  Author: Zhou Yizhuang, zhouyizhuang3@163.com
  Version: 2.0,  Date: 2021-12-05
 ------------------------------------------------------------------------------------------------------------

=head1 Command-line Option
 
 ------------------------------------------------------------------------------------------------------------
  perl $0 [Options]
  --Qfile    String    Input file for query tested genomes
  --Rfile    String    Input file for reference tested genomes
  --Afile    String    Input file for all testded genomes including query and reference genomes, 
                       intead for the combination of "--Qfile" and "--Rfile"
  --QPrefix  String    The string for prefix of the output files for queries,default is "Query"
  --RPrefix  String    The string for prefix of the output files for references,default is "Ref"
  --APrefix  String    The string for prefix of the output files for references,default is "All"
  --Outdir   String    The directory for output files, the default is the currect working directory
  --verbose            output verbose information to screen  
  --h or help          Display this message
 ------------------------------------------------------------------------------------------------------------

=head1 Usage Exmples
  
  -----------------------------------------------------------------------------------------------------------
  perl ./FRAGTE2.pl --Qfile  <Qfile> --Rfile <Rfile>
  perl ./FRAGTE2.pl --Qfile  <Qfile> --Rfile <Rfile> --QPrefix <Prefix> --RPrefix <Prefix>
  perl ./FRAGTE2.pl --Qfile  <Qfile> --Rfile <Rfile> --QPrefix <Prefix> --RPrefix <Prefix> --Outdir <outdir> 
 ------------------------------------------------------------------------------------------------------------

=head1 Note
 
 ------------------------------------------------------------------------------------------------------------
  (1)The input file for "--Qfile" or "--Rfile" has 7 fields separatted by Tab, including: 
  Assembly_accession,species_taxid,organism_name,Average Size,assembly_level,total genome size (in bp)
  and file for genome
  -----------------------------------------------------------------------------------------------------------
  (2) file for genome,the file with absolute path containing sequences in fasta format.
  ------------------------------------------------------------------------------------------------------------
  (3) "--Qfile" and "--Rfile" should be used together. Otherwise, "--Afile" is used instead which includes both 
  queries and references.
  -----------------------------------------------------------------------------------------------------------
  (4) the output file named "*_Pairs_byFRAGTE2.xls" for pairs sieved by FRAGTE2 includes 5 fields separatted by
  Tab as follows: Query ID, Reference ID, PCCD, GSCq and GSCr.
  ------------------------------------------------------------------------------------------------------------
  (5) the output file named "*_Pairs_byFRAGTE2.log" for pairs unsieve by FRAGTE2 also includes 5 fields separatted
  includes 5 fields separatted by Table as follows: Query ID, Reference ID, PCCD, GSCq and GSCr.
 ------------------------------------------------------------------------------------------------------------

=head1 Please cite

 ------------------------------------------------------------------------------------------------------------
 If you use FRAGTE2 in your publication, please cite:
 Zhou, Y., et al., A completeness-independent method for pre-selection of closely related genomes for species
 delineation in prokaryotes. BMC Genomics, 2020. 21(1): p. 183.
 ------------------------------------------------------------------------------------------------------------
=cut

die `pod2text $0` if (@ARGV == 0);
## define variable 
my $Qfile;
my $Rfile;
my $Afile;
my $QPrefix;
my $RPrefix;
my $APrefix;
my $Outdir;
my $verbose;
my $Help;
## parse options from @ARGV
GetOptions(
	"Qfile:s"=>\$Qfile,         ## string
	"Rfile:s"=>\$Rfile,         ## string
	"Afile:s"=>\$Afile,         ## string
	"QPrefix:s"=>\$QPrefix,     ## string
	"RPrefix:s"=>\$RPrefix,     ## string
	"APrefix:s"=>\$RPrefix,     ## string
	"Outdir:s"=>\$Outdir,       ## string
	"verbose"=>\$verbose,       ## flag
	"help|h"=>\$Help            ## flag
) || die "Please use --help option to get help\n";

## option variables with default value
$QPrefix ||= "Query";
$RPrefix ||= "Ref";
$APrefix ||= "All";
$Outdir ||= getcwd();

$Outdir = abs_path($Outdir);
mkdir $Outdir,0755 unless (-e $Outdir);

die `pod2text $0` if ($Help);
if($Qfile && !$Rfile){
	print STDERR "\n";
	print STDERR "Error: --Rfile option is also required\n";
	print STDERR "\n";
	print STDERR "The following is the manual:\n";
	print STDERR "\n";
	die `pod2text $0`;
}
if($Rfile && !$Qfile){
	print STDERR "\n";
	print STDERR "Error: --Qfile option is also required\n";
	print STDERR "\n";
	print STDERR "The following is the manual:\n";
	print STDERR "\n";
	die `pod2text $0`;
}
if(($Rfile || $Qfile) && $Afile){
	print STDERR "\n";
	print STDERR "Error: --Rfile and --Qfile options should be excluded\n";
	print STDERR "\n";
	print STDERR "The following is the manual:\n";
	print STDERR "\n";
	die `pod2text $0`;
}

my $start_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
my $t1=timelocal(localtime());
if(!$verbose){
	open(LOG,">$Outdir/Pairs_byFRAGTE2.log")||die;
}
if(!$verbose){
	print LOG "$start_time: Start FRAGTE2.pl ......\n";
}
else{
	print STDERR "$start_time: Start FRAGTE2.pl ......\n";
}

######################################################## Intializing  ##############################################
my $current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
if(!$verbose){
	print LOG "$current_time: Start Intializing\n";
}
else{
	print STDERR "$current_time: Start Intializing\n";
}
my $km = 4;
my @mono = ("A","T","G","C");
my @nuc = @mono;
#for dinucleotide
my @oligo_k_2mer = (); 
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
	@nuc=@tmpnuc;
	if($index == $km-2){
		@oligo_k_2mer = @nuc;
		last;
	}
}
my @UsedKmer2=();
my %tag=();
foreach my $k1 (@oligo_k_2mer) {
	my $k2=reverse($k1);
	$k2=~tr/ATGC/TACG/;
	if($tag{$k2}){
		next;
	}
	else{
		push @UsedKmer2,$k1;
		$tag{$k1}=1;
	}
}

#for trinucleotide
my @oligo_k_1mer = (); 
$index = 1;
@nuc = @mono;
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
	if($index == $km-1){
		@oligo_k_1mer = @nuc;
		last;
	}
}
my @UsedKmer1=();
%tag=();
foreach my $k1 (@oligo_k_1mer) {
	my $k2=reverse($k1);
	$k2=~tr/ATGC/TACG/;
	if($tag{$k2}){
		next;
	}
	else{
		push @UsedKmer1,$k1;
		$tag{$k1}=1;
	}
}

#for tetranucleotide
my @oligo_kmer = (); 
$index = 1;
@nuc = @mono;
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
%tag=();
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

$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
if(!$verbose){
	print LOG "$current_time: Finish Intializing\n";
}
else{
	print STDERR "$current_time: Finish Intializing\n";
}

####################################################### Reading Cutoff #############################################
open(TMP,"$Bin/Cutoff.xls")||die;
<TMP>;
my %LSC;
my %LSC2;
while(<TMP>){
	chomp;
	my @a=split /\t/;
	my @b=split /\s*vs\s*/,$a[0];
	my $LSC2=$a[5];
	if($LSC2=~/(0\.\d{2})/){
		$LSC2{$b[0]}{$b[1]}=$1;
		$LSC2{$b[1]}{$b[0]}=$1;
	}
	if($b[0] eq $b[1]){
		my $cutoff=$a[1]-$a[2]-0.01;
		if($cutoff=~/(0\.\d{2})/){
			$LSC{$b[0]}=$1;
		}
	}
}
close TMP;
$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
if(!$verbose){
	print LOG "$current_time: Finish reading $Bin/Cutoff.xls\n";
}
else{
	print STDERR "$current_time: Finish reading $Bin/Cutoff.xls\n";
}

########################################### Fragmenting phase for references #######################################
if($Qfile && $Rfile){
	$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
	my $rt1=timelocal(localtime());
	if(!$verbose){
		print LOG "$current_time: Start fragmenting phase for references\n";
	}
	else{
		print STDERR "$current_time: Start fragmenting phase for references\n";
	}
	open(RF,$Rfile)||die;
	open(OR,">$Outdir/$RPrefix\_RepZvalue.xls")||die;
	while(<RF>){
		chomp;
		my @a=split /\t/;
		next if($a[5]<10000);
		my $id=$a[0];
		my $spid=$a[1];
		my $AgeLen=$a[3];
		my $file=$a[6];
		open(FA,$file)||die;
		my $MAS="";
		my $SeqID="";
		my $seq="";
		my $Size=0;
		while(<FA>){
			chomp;
			if(/^>(\S+)/){
				if($SeqID ne ""){
					my $MASlen=length $MAS;
					if($MASlen>=1000){
						$seq.=$MAS;
						$Size+=$MASlen;
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
			$seq.=$MAS;
			$Size+=$MASlen;
		}
		next if($Size<10000);

		if($Size <40000){
			my $kb=int($Size/10000)*10;
			my $threshold=$LSC{$kb};
			my $ref=&Fragment1($seq,$Size);
			print OR "$id\t$Size\t$kb\t$threshold\t$ref\n";
		}
		else{
			my ($threshold,$ref,$kb)=&Fragment2($id,$seq,$Size,$AgeLen);
			print OR "$id\t$Size\t$kb\t$threshold\t$ref\n";
		}
	}
	close RF;
	close OR;
	$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
	my $rt2=timelocal(localtime());
	if(!$verbose){
		print LOG "$current_time: Finish fragmenting phase for references\n";
		print LOG "$current_time: Produce $Outdir/$RPrefix\_RepZvalue.xls\n";
	}
	else{
		print STDERR "$current_time: Finish fragmenting phase for references\n";
		print STDERR "$current_time: Produce $Outdir/$RPrefix\_RepZvalue.xls\n";
	}
	my $rs=$rt2-$rt1;
	my $rhour=int($rs / 3600);
	my $rms=$rs % 3600;
	my $rminutes=int($rms / 60);
	my $rsecond=$rms % 60;
	if(!$verbose){
		print LOG "Time elaspe in seconds for fragmenting phase for references: $rs\n";
		print LOG "Time elaspe for fragmenting phase for references: $rhour\:$rminutes\:$rsecond\n";
	}
	else{
		print STDERR "Time elaspe in seconds for fragmenting phase for references: $rs\n";
		print STDERR "Time elaspe for fragmenting phase for references: $rhour\:$rminutes\:$rsecond\n";
	}

	########################################### Fragmenting phase for queries #####################################
	$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
	my $qt1=timelocal(localtime());
	if(!$verbose){
		print LOG "$current_time: Start fragmenting phase for queries\n";
	}
	else{
		print STDERR "$current_time: Start fragmenting phase for queries\n";
	}
	open(QF,$Qfile)||die;
	open(OQ,">$Outdir/$QPrefix\_RepZvalue.xls")||die;
	while(<QF>){
		chomp;
		my @a=split /\t/;
		next if($a[5]<10000);
		my $id=$a[0];
		my $spid=$a[1];
		my $AgeLen=$a[3];
		my $file=$a[6];
		open(FA,$file)||die;
		my $MAS="";
		my $SeqID="";
		my $seq="";
		my $Size=0;
		while(<FA>){
			chomp;
			if(/^>(\S+)/){
				if($SeqID ne ""){
					my $MASlen=length $MAS;
					if($MASlen>=1000){
						$seq.=$MAS;
						$Size+=$MASlen;
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
			$seq.=$MAS;
			$Size+=$MASlen;
		}
		next if($Size<10000);

		if($Size <40000){
			my $kb=int($Size/10000)*10;
			my $threshold=$LSC{$kb};
			my $ref=&Fragment1($seq,$Size);
			print OQ "$id\t$Size\t$kb\t$threshold\t$ref\n";
		}
		else{
			my ($threshold,$ref,$kb)=&Fragment2($id,$seq,$Size,$AgeLen);
			print OQ "$id\t$Size\t$kb\t$threshold\t$ref\n";
		}
	}
	close QF;
	close OQ;
	$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
	my $qt2=timelocal(localtime());
	if(!$verbose){
		print LOG "$current_time: Finish fragmenting phase for queries\n";
		print LOG "$current_time: Produce $Outdir/$QPrefix\_RepZvalue.xls\n";
	}
	else{
		print STDERR "$current_time: Finish fragmenting phase for queries\n";
		print STDERR "$current_time: Produce $Outdir/$QPrefix\_RepZvalue.xls\n";
	}
	my $qs=$qt2-$qt1;
	my $qhour=int($qs / 3600);
	my $qms=$qs % 3600;
	my $qminutes=int($qms / 60);
	my $qsecond=$qms % 60;
	if(!$verbose){
		print LOG "Time elaspe in seconds for fragmenting phase for queries: $qs\n";
		print LOG "Time elaspe for fragmenting phase for queries: $qhour\:$qminutes\:$qsecond\n";
	}
	else{
		print STDERR "Time elaspe in seconds for fragmenting phase for queries: $qs\n";
		print STDERR "Time elaspe for fragmenting phase for queries: $qhour\:$qminutes\:$qsecond\n";
	}
	############################################### Determining phase  ###############################################
	$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
	my $dt1=timelocal(localtime());
	if(!$verbose){
		print LOG "$current_time: Start determining phase\n";
	}
	else{
		print STDERR "$current_time: Start determining phase\n";
	}
	my @ID=();
	my %Size=();
	my %kb=();
	my %GSC=();
	my $file="$Outdir/$RPrefix\_RepZvalue.xls";
	my $r=Statistics::R->new();
	$r->set("file",$file);
	$r->run(qq'd<-read.delim(file,header=F)
	    zvalue<-d[,5:140]');
	open(IN,$file)||die;
	while(<IN>){
		chomp;
		my @a=split /\t/;
		my $id=shift @a;
		push @ID,$id;
		my $Size=shift @a;
		$Size{$id}=$Size;
		my $kb=shift @a;
		$kb{$id}=$kb;
		my $GSC=shift @a;
		$GSC{$id}=$GSC;
	}
	close IN;

	open(IF,"$Outdir/$QPrefix\_RepZvalue.xls")||die;
	open(OUT,">$Outdir/$QPrefix\_$RPrefix\_Pairs_byFRAGTE2.xls")||die;
	open(OF,">$Outdir/$QPrefix\_$RPrefix\_Pairs_byFRAGTE2.log")||die;
	while(<IF>){
		chomp;
		my @a=split /\t/;
		my $id1=shift @a;
		my $Size1=shift @a;
		my $kb1=shift @a;
		my $GSC1=shift @a;
		$r->set("query",\@a);
		$r->run(qq'PCCD<-apply(zvalue,1,function(i){cor(x=query,y=i)})
            PCCD<-round(PCCD,digits=2)');
		my @PCCD=();
		if(@ID <=1){
			my $PCCD=$r->get("PCCD");
			push @PCCD,$PCCD;
		}
		else{
			my $PCCDref=$r->get("PCCD");
			@PCCD=@{$PCCDref};
		}
		my $flag =0;
		my %PCCD=();
		my @tmpID=();
		my %tag=();
		for(my $i=0;$i<=$#ID;$i++){
			my $id2=$ID[$i];
			my $Size2=$Size{$id2};
			my $kb2=$kb{$id2};
			my $GSC2=$GSC{$id2};
			my $PCCD=$PCCD[$i];
			if($PCCD >=$LSC2{$kb1}{$kb2}){
				push @tmpID,$id2;
				$PCCD{$id2}=$PCCD;
			}
			if($PCCD>=$GSC1 || $PCCD>=$GSC2){
				$flag++;
				$tag{$id2} = 1;
				print OUT "$id1\t$id2\t$PCCD\t$GSC1\t$GSC2\n";
			}
			else{
				next;
			}
		}
		if($flag == 0){
			@tmpID=sort {$PCCD{$b} <=> $PCCD{$a}} @tmpID;
			my $n=0;
			my $Fid=$tmpID[99];
			foreach my $id2 (@tmpID) {
				$n++;
				my $GSC2=$GSC{$id2};
				if($n<=100){
					print OUT "$id1\t$id2\t$PCCD{$id2}\t$GSC1\t$GSC2\n";
				}
				elsif($PCCD{$id2} == $PCCD{$Fid}){
					print OUT "$id1\t$id2\t$PCCD{$id2}\t$GSC1\t$GSC2\n";
				}
				else{
					print OF "$id1\t$id2\t$PCCD{$id2}\t$GSC1\t$GSC2\n";
				}
			}
		}
		else{
			foreach my $id2 (@tmpID) {
				if($tag{$id2}){
					next;
				}
				else{
					my $GSC2=$GSC{$id2};
					print OF "$id1\t$id2\t$PCCD{$id2}\t$GSC1\t$GSC2\n";
				}
			}
		}
	}
	close IN;
	close OUT;
	close OF;
	$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
	my $t2=timelocal(localtime());
	if(!$verbose){
		print LOG "$current_time: Finish determining phase\n";
		print LOG "$current_time: Produce $Outdir/$QPrefix\_$RPrefix\_Pairs_byFRAGTE2.xls\n";
		print LOG "$current_time: Produce $Outdir/$QPrefix\_$RPrefix\_Pairs_byFRAGTE2.log\n";
	}
	else{
		print STDERR "$current_time: Finish determining phase\n";
		print STDERR "$current_time: Produce $Outdir/$QPrefix\_$RPrefix\_Pairs_byFRAGTE2.xls\n";
		print STDERR "$current_time: Produce $Outdir/$QPrefix\_$RPrefix\_Pairs_byFRAGTE2.log\n";
	}
	my $s=$t2-$dt1;
	my $hour=int($s / 3600);
	my $ms=$s % 3600;
	my $minutes=int($ms / 60);
	my $second=$ms % 60;
	if(!$verbose){
		print LOG "Time elaspe in seconds for determining phase: $s\n";
		print LOG "Time elaspe for determining phase: $hour\:$minutes\:$second\n";
	}
	else{
		print STDERR "Time elaspe in seconds for determining phase: $s\n";
		print STDERR "Time elaspe for determining phase: $hour\:$minutes\:$second\n";
	}
	$s=$t2-$t1;
	$hour=int($s / 3600);
	$ms=$s % 3600;
	$minutes=int($ms / 60);
	$second=$ms % 60;
	if(!$verbose){
		print LOG "Time elaspe in seconds for all: $s\n";
		print LOG "Time elaspe for all: $hour\:$minutes\:$second\n";
	}
	else{
		print STDERR "Time elaspe in seconds for all: $s\n";
		print STDERR "Time elaspe for all: $hour\:$minutes\:$second\n";
	}
	if(!$verbose){
		close LOG;
	}
}
elsif($Afile){
	########################################### Fragmenting phase for All seqs #####################################
	$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
	my $qt1=timelocal(localtime());
	if(!$verbose){
		print LOG "$current_time: Start fragmenting phase for all seqs\n";
	}
	else{
		print STDERR "$current_time: Start fragmenting phase for all seqs\n";
	}
	open(AF,$Afile)||die;
	open(OA,">$Outdir/$APrefix\_RepZvalue.xls")||die;
	while(<AF>){
		chomp;
		my @a=split /\t/;
		next if($a[5]<10000);
		my $id=$a[0];
		my $spid=$a[1];
		my $AgeLen=$a[3];
		my $file=$a[6];
		open(FA,$file)||die;
		my $MAS="";
		my $SeqID="";
		my $seq="";
		my $Size=0;
		while(<FA>){
			chomp;
			if(/^>(\S+)/){
				if($SeqID ne ""){
					my $MASlen=length $MAS;
					if($MASlen>=1000){
						$seq.=$MAS;
						$Size+=$MASlen;
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
			$seq.=$MAS;
			$Size+=$MASlen;
		}
		next if($Size<10000);

		if($Size <40000){
			my $kb=int($Size/10000)*10;
			my $threshold=$LSC{$kb};
			my $ref=&Fragment1($seq,$Size);
			print OA "$id\t$Size\t$kb\t$threshold\t$ref\n";
		}
		else{
			my ($threshold,$ref,$kb)=&Fragment2($id,$seq,$Size,$AgeLen);
			print OA "$id\t$Size\t$kb\t$threshold\t$ref\n";
		}
	}
	close AF;
	close OA;
	$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
	my $qt2=timelocal(localtime());
	if(!$verbose){
		print LOG "$current_time: Finish fragmenting phase for queries\n";
		print LOG "$current_time: Produce $Outdir/$APrefix\_RepZvalue.xls\n";
	}
	else{
		print STDERR "$current_time: Finish fragmenting phase for queries\n";
		print STDERR "$current_time: Produce $Outdir/$APrefix\_RepZvalue.xls\n";
	}
	my $qs=$qt2-$qt1;
	my $qhour=int($qs / 3600);
	my $qms=$qs % 3600;
	my $qminutes=int($qms / 60);
	my $qsecond=$qms % 60;
	if(!$verbose){
		print LOG "Time elaspe in seconds for fragmenting phase for all seqs: $qs\n";
		print LOG "Time elaspe for fragmenting phase for all seqs: $qhour\:$qminutes\:$qsecond\n";
	}
	else{
		print STDERR "Time elaspe in seconds for fragmenting phase for all seqs: $qs\n";
		print STDERR "Time elaspe for fragmenting phase for all seqs: $qhour\:$qminutes\:$qsecond\n";
	}

	############################################### Determining phase  ###############################################
	$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
	my $dt1=timelocal(localtime());
	if(!$verbose){
		print LOG "$current_time: Start determining phase\n";
	}
	else{
		print STDERR "$current_time: Start determining phase\n";
	}
	open(OUT,">$Outdir/$APrefix\_Pairs_byFRAGTE2.xls")||die;
	open(OF,">$Outdir/$APrefix\_Pairs_byFRAGTE2.log")||die;
	my $file="$Outdir/$APrefix\_RepZvalue.xls";
	my $r=Statistics::R->new();
	$r->set("file",$file);
	$r->run(qq'd<-read.delim(file,header=F)
			zvalue<-d[,5:140]');
	open(IN,$file)||die;
	my @ID=();
	my %Size=();
	my %kb=();
	my %GSC=();
	while(<IN>){
		chomp;
		my @a=split /\t/;
		my $id=shift @a;
		push @ID,$id;
		my $Size=shift @a;
		$Size{$id}=$Size;
		my $kb=shift @a;
		$kb{$id}=$kb;
		my $GSC=shift @a;
		$GSC{$id}=$GSC;
	}
	close IN;
	open(IN,$file)||die;
	while(<IN>){
		chomp;
		my @a=split /\t/;
		my $id1=shift @a;
		my $Size1=shift @a;
		my $kb1=shift @a;
		my $GSC1=shift @a;
		$r->set("query",\@a);
		$r->run(qq'PCCD<-apply(zvalue,1,function(i){cor(x=query,y=i)})
				PCCD<-round(PCCD,digits=2)');
		my @PCCD=();
		if(@ID <=1){
			my $PCCD=$r->get("PCCD");
			push @PCCD,$PCCD;
		}
		else{
			my $PCCDref=$r->get("PCCD");
			@PCCD=@{$PCCDref};
		}
		my $flag =0;
		my %PCCD=();
		my @tmpID=();
		my %tag=();
		for(my $i=0;$i<=$#ID;$i++){
			my $id2=$ID[$i];
			my $Size2=$Size{$id2};
			my $kb2=$kb{$id2};
			my $GSC2=$GSC{$id2};
			my $PCCD=$PCCD[$i];
			next if($id1 eq $id2);
			if($PCCD >=$LSC2{$kb1}{$kb2}){
				push @tmpID,$id2;
				$PCCD{$id2}=$PCCD;
			}
			if($PCCD>=$GSC1 || $PCCD>=$GSC2){
				$flag++;
				print OUT "$id1\t$id2\t$PCCD\t$GSC1\t$GSC2\n";
				$tag{$id2}=1;
			}
			else{
				next;
			}
		}
		if($flag == 0){
			@tmpID=sort {$PCCD{$b} <=> $PCCD{$a}} @tmpID;
			my $n=0;
			my $Fid=$tmpID[99];
			foreach my $id2 (@tmpID) {
				$n++;
				my $GSC2=$GSC{$id2};
				if($n<=100){
					print OUT "$id1\t$id2\t$PCCD{$id2}\t$GSC1\t$GSC2\n";
				}
				elsif($PCCD{$id2} == $PCCD{$Fid}){
					print OUT "$id1\t$id2\t$PCCD{$id2}\t$GSC1\t$GSC2\n";
				}
				else{
					print OF "$id1\t$id2\t$PCCD{$id2}\t$GSC1\t$GSC2\n";
				}
			}
		}
		else{
			foreach my $id2 (@tmpID) {
				if($tag{$id2}){
					next;
				}
				else{
					my $GSC2=$GSC{$id2};
					print OF "$id1\t$id2\t$PCCD{$id2}\t$GSC1\t$GSC2\n";
				}
			}
		}
	}
	close IN;
	close OUT;
	close OF;
	$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
	my $t2=timelocal(localtime());
	if(!$verbose){
		print LOG "$current_time: Finish determining phase\n";
		print LOG "$current_time: Produce $Outdir/$APrefix\_Pairs_byFRAGTE2.xls\n";
		print LOG "$current_time: Produce $Outdir/$APrefix\_Pairs_byFRAGTE2.log\n";
	}
	else{
		print STDERR "$current_time: Finish determining phase\n";
		print STDERR "$current_time: Produce $Outdir/$APrefix\_Pairs_byFRAGTE2.xls\n";
		print STDERR "$current_time: Produce $Outdir/$APrefix\_Pairs_byFRAGTE2.log\n";
	}
	my $s=$t2-$dt1;
	my $hour=int($s / 3600);
	my $ms=$s % 3600;
	my $minutes=int($ms / 60);
	my $second=$ms % 60;
	if(!$verbose){
		print LOG "Time elaspe in seconds for determining phase: $s\n";
		print LOG "Time elaspe for determining phase: $hour\:$minutes\:$second\n";
	}
	else{
		print STDERR "Time elaspe in seconds for determining phase: $s\n";
		print STDERR "Time elaspe for determining phase: $hour\:$minutes\:$second\n";
	}
	$s=$t2-$t1;
	$hour=int($s / 3600);
	$ms=$s % 3600;
	$minutes=int($ms / 60);
	$second=$ms % 60;
	if(!$verbose){
		print LOG "Time elaspe in seconds for all: $s\n";
		print LOG "Time elaspe for all: $hour\:$minutes\:$second\n";
	}
	else{
		print STDERR "Time elaspe in seconds for all: $s\n";
		print STDERR "Time elaspe for all: $hour\:$minutes\:$second\n";
	}
	if(!$verbose){
		close LOG;
	}
}
else{
	print STDERR "\n";
	print STDERR "Error: --Afile (or --Rfile and --Qfile) should be given\n";
	print STDERR "\n";
	print STDERR "The following is the manual:\n";
	print STDERR "\n";
	die `pod2text $0`;
}

############################################################## Computing zvalue ##########################################################
sub zvalue{
	my ($ref1,$ref2,$ref3)=@_;
	my %kmer = ();
	my %k_1mer = ();
	my %k_2mer = ();
	%kmer=%{$ref1};
	%k_1mer=%{$ref2};
	%k_2mer=%{$ref3};
	my @Zvalue=();
	foreach my $k1 (@UsedKmer){
		my $k2=reverse($k1);
		$k2=~tr/ATGC/TACG/;
		my $N_koligo=$kmer{$k1};
		my $YZ=substr($k1,1,$km-2);
		my $XYZ=substr($k1,0,$km-1);
		my $YZW=substr($k1,1,$km-1);
		my $N_former=$k_1mer{$XYZ};
		my $N_latter=$k_1mer{$YZW};
		my $N_midder=$k_2mer{$YZ};
		my $denominator=$N_midder-$N_former;
		$denominator *= ($N_midder-$N_latter);
		$denominator *= $N_former;
		$denominator *= $N_latter;
		my $sqrt=$denominator**0.5;
		if($sqrt !=0 and $N_midder !=0){
			#my $zvalue=$N_midder**0.5*($N_koligo*$N_midder-$N_former*$N_latter)/$sqrt;
			my $zvalue=$N_midder**0.5;
			$zvalue *= ($N_koligo*$N_midder-$N_former*$N_latter);
			$zvalue /= $sqrt;
			push @Zvalue,$zvalue;
		}
		else{
			push @Zvalue,0;
		}
	}
	return \@Zvalue;
}

###################################################### Subroutine for computing mean##############################################
sub mean{
	my ($ref)=@_;
	my @a=@{$ref};
	my $sum=0;
	foreach(@a){
		$sum+=$_;
	}
	my $mean=$sum/@a;
	return $mean;
}

####################################################################### Fragment #####################################################
sub Fragment1{
	my ($seq,$len)=@_;
	my %TNF=();
	for(my $j=0;$j<=$len-$km;$j++){
		my $sub=substr($seq,$j,$km);
		$TNF{$sub}++;
	}

	my %SUM = ();
	my %SUM1 = ();
	my %SUM2 = ();
	foreach my $k1 (@UsedKmer){
		my $k2=reverse($k1);
		$k2=~tr/ATGC/TACG/;
		if($TNF{$k1}){
			$SUM{$k1}+=$TNF{$k1};
		}
		if($TNF{$k2}){
			$SUM{$k1}+=$TNF{$k2};
		}
		if(!$SUM{$k1}){
			$SUM{$k1}=1;
		}
		$SUM{$k2}=$SUM{$k1};
	}
	foreach my $kmer(@UsedKmer) {
		for(my $j=0;$j<=$km-($km-1);$j++){
			my $sub=substr($kmer,$j,$km-1);
			if($SUM{$kmer} >1){
				$SUM1{$sub}+=$SUM{$kmer};
			}
		}
		for(my $j=0;$j<=$km-($km-2);$j++){
			my $sub=substr($kmer,$j,$km-2);
			if($SUM{$kmer}>1){
				$SUM2{$sub}+=$SUM{$kmer};
			}
		}
	}

	foreach my $kmer(@oligo_kmer) {
		for(my $j=0;$j<=$km-($km-1);$j++){
			my $sub=substr($kmer,$j,$km-1);
			if($SUM{$kmer} >1){
				$SUM1{$sub}+=$SUM{$kmer};
			}
		}
		for(my $j=0;$j<=$km-($km-2);$j++){
			my $sub=substr($kmer,$j,$km-2);
			if($SUM{$kmer}>1){
				$SUM2{$sub}+=$SUM{$kmer};
			}
		}
	}
	foreach my $k (@oligo_k_1mer){
		if(!$SUM1{$k}){
			$SUM1{$k}=1;
		}
		else{
			$SUM1{$k}=int(($SUM1{$k}+1)/2);
		}
	}
	foreach my $k (@oligo_k_2mer){
		if(!$SUM2{$k}){
			$SUM2{$k}=1;
		}
		else{
			$SUM2{$k}=int(($SUM2{$k}+2)/3);
		}
	}
	my $zvalueref=&zvalue(\%SUM,\%SUM1,\%SUM2);
	my $mean=&mean($zvalueref);
	my @f=();
	foreach(@{$zvalueref}){
		my $tmp=$_-$mean;
		push @f,$tmp;
	}
	my $refxx=\@f;
	my $ref=join("\t",@f);
	return ($ref);
}

####################################################################### Fragment #####################################################
sub Fragment2{
	my ($id,$seq,$len,$AgeLen)=@_;
	my $ref;
	my $threshold;
	my $sublen=int($len/4);
	if($AgeLen < 200000){
		if($sublen >200000){
			$sublen=200000;
		}
	}
	my $addlen=int($sublen/2);

	###TNF
	my %TNF=();
	my $index=0;
	my @index=();

	for(my $i=0;$i<=$len-$addlen;$i+=$addlen){
		$index++;
		push @index,$index;
		my $fragment=substr($seq,$i,$addlen);
		for(my $j=0;$j<=$addlen-$km;$j++){
			my $sub=substr($fragment,$j,$km);
			$TNF{$index}{$sub}++;
		}
	}

	## To obtain the representative fragment
	my $index2="";
	my $r=Statistics::R->new();
	$r->run(q'zvalue<-c()');
	for (my $ind=0;$ind<$#index;$ind++){
		my $index=$index[$ind];
		my %kmer=();
		my %k_1mer=();
		my %k_2mer=();
		### To obtain km frequency
		foreach my $k1 (@UsedKmer){
			my $k2=reverse($k1);
			$k2=~tr/ATGC/TACG/;
			foreach (0 .. 1) {
				$index2=$index+$_;
				if($index2 == @index +1){#To link the start and end sequences
					$index2 =1;
				}
				if($TNF{$index2}{$k1}){
					$kmer{$k1}+=$TNF{$index2}{$k1};
				}
				if($TNF{$index2}{$k2}){
					$kmer{$k1}+=$TNF{$index2}{$k2};
				}
			}
			if(! $kmer{$k1}){
				$kmer{$k1}=1;
			}
			$kmer{$k2}=$kmer{$k1};
		}
			
		### To obtain km-1 and km-2 frequencies
		foreach my $kmer(@oligo_kmer) {
			for(my $j=0;$j<=$km-($km-1);$j++){
				my $sub=substr($kmer,$j,$km-1);
				if($kmer{$kmer} >1){
					$k_1mer{$sub}+=$kmer{$kmer};
				}
			}
			for(my $j=0;$j<=$km-($km-2);$j++){
				my $sub=substr($kmer,$j,$km-2);
				if($kmer{$kmer}>1){
					$k_2mer{$sub}+=$kmer{$kmer};
				}
			}
		}
		foreach my $k (@oligo_k_1mer){
			if(!$k_1mer{$k}){
				$k_1mer{$k}=1;
			}
			else{
				$k_1mer{$k}=int(($k_1mer{$k}+1)/2);
			}
		}
		foreach my $k (@oligo_k_2mer){
			if(!$k_2mer{$k}){
				$k_2mer{$k}=1;
			}
			else{
				$k_2mer{$k}=int(($k_2mer{$k}+2)/3);
			}
		}
		my $zvalueref=&zvalue(\%kmer,\%k_1mer,\%k_2mer);
		$r->set("x",$zvalueref);
		$r->run(q'zvalue<-append(zvalue,x)');
	}
	$r->run(qq'mz<-matrix(zvalue,nrow=136,byrow=F)
			ZRF<-apply(mz,1,function(i){
				if(ncol(mz) %% 2 == 0){
					TNFmean<-mean(i,trim=(ncol(mz)-4)/ncol(mz))
				}else{
					TNFmean<-mean(i,trim=(ncol(mz)-3)/ncol(mz))
				}
				return(TNFmean)
			})
			PCCD<-apply(mz,2,function(i){cor(x=ZRF,y=i)})
			GSC<-floor((mean(PCCD)-sd(PCCD))*100)/100-0.01');
	my $ZRFref=$r->get('ZRF');
	my @ZRF=@{$ZRFref};
	$threshold=$r->get("GSC");
	my $kb=int($len/10000)*10;
	if($kb > 200){
		$kb = 200;
	}
	my $current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
	if(!$verbose){
		print LOG "$current_time: $id\t$threshold\t$LSC{$kb}\t$LSC2{$kb}{$kb}\n";
	}
	else{
		print STDERR "$current_time: $id\t$threshold\t$LSC{$kb}\t$LSC2{$kb}{$kb}\n";
	}
	$ref=join("\t",@ZRF);
	if($threshold <=$LSC{$kb}){
		if($threshold<=$LSC2{$kb}{$kb}){
			$addlen=10000;
			if($sublen/$addlen >=2){
				###TNF
				my %TNF=();
				my $index=0;
				my @index=();
				for(my $i=0;$i<=$len-$addlen;$i+=$addlen){
					$index++;
					push @index,$index;
					my $fragment=substr($seq,$i,$addlen);
					for(my $j=0;$j<=$addlen-$km;$j++){
						my $sub=substr($fragment,$j,$km);
						$TNF{$index}{$sub}++;
					}
				}

				## To obtain the representative fragment
				my $index2="";
				my $r=Statistics::R->new();
				$r->run(q'zvalue<-c()');
				for (my $ind=0;$ind<$#index;$ind++){
					my $index=$index[$ind];
					my %kmer=();
					my %k_1mer=();
					my %k_2mer=();
					### To obtain km frequency
					foreach my $k1 (@UsedKmer){
						my $k2=reverse($k1);
						$k2=~tr/ATGC/TACG/;
						foreach (0 .. 1) {
							$index2=$index+$_;
							if($index2 == @index +1){#To link the start and end sequences
								$index2 =1;
							}
							if($TNF{$index2}{$k1}){
								$kmer{$k1}+=$TNF{$index2}{$k1};
							}
							if($TNF{$index2}{$k2}){
								$kmer{$k1}+=$TNF{$index2}{$k2};
							}
						}
						if(! $kmer{$k1}){
							$kmer{$k1}=1;
						}
						$kmer{$k2}=$kmer{$k1};
					}
						
					### To obtain km-1 and km-2 frequencies
					foreach my $kmer(@oligo_kmer) {
						for(my $j=0;$j<=$km-($km-1);$j++){
							my $sub=substr($kmer,$j,$km-1);
							if($kmer{$kmer} >1){
								$k_1mer{$sub}+=$kmer{$kmer};
							}
						}
						for(my $j=0;$j<=$km-($km-2);$j++){
							my $sub=substr($kmer,$j,$km-2);
							if($kmer{$kmer}>1){
								$k_2mer{$sub}+=$kmer{$kmer};
							}
						}
					}
					foreach my $k (@oligo_k_1mer){
						if(!$k_1mer{$k}){
							$k_1mer{$k}=1;
						}
						else{
							$k_1mer{$k}=int(($k_1mer{$k}+1)/2);
						}
					}
					foreach my $k (@oligo_k_2mer){
						if(!$k_2mer{$k}){
							$k_2mer{$k}=1;
						}
						else{
							$k_2mer{$k}=int(($k_2mer{$k}+2)/3);
						}
					}
					my $zvalueref=&zvalue(\%kmer,\%k_1mer,\%k_2mer);
					$r->set("x",$zvalueref);
					$r->run(q'zvalue<-append(zvalue,x)');
				}
				$r->run(qq'mz<-matrix(zvalue,nrow=136,byrow=F)
						ZRF<-apply(mz,1,function(i){
							if(ncol(mz) %% 2 == 0){
								TNFmean<-mean(i,trim=(ncol(mz)-4)/ncol(mz))
							}else{
								TNFmean<-mean(i,trim=(ncol(mz)-3)/ncol(mz))
							}
							return(TNFmean)
						})');
				my $ZRFref=$r->get('ZRF');
				my @ZRF=@{$ZRFref};
				$ref=join("\t",@ZRF);
			}
		}
		$threshold =$LSC{$kb};
	}
	return($threshold,$ref,$kb);
}
__END__