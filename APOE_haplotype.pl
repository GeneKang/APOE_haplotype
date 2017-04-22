#!/usr/bin/perl 
use strict;

# APOE haplotype 
# www.genekang.com 
# qjmin@genekang.com  闵庆杰

#chr19   45411791
#GCGCTGATGGACGAGACCATGAAGGAGTTGAAGGCCTACAAATCGGAACTGGAGGAACAACTGACCCCGGTGGCGGAGGAGACGCGGGCACGGCTGTCCAAGGAGCTGCAGGCGGCGCAGGCCCGGCTGGGCGCGGACATGGAGGACGTGTGCGGCCGCCTGGTGCAGTACCGCGGCGAGGTGCAGGCCATGCTCGGCCAGAGCACCGAGGAGCTGCGGGTGCGCCTCGCCTCCCACCTGCGCAAGCTGCGTAAGCGGCTCCTCCGCGATGCCGATGACCTGCAGAAGCGCCTGGCAGTGTACCAGGCCGGGGCCCGCGAGGGCGCCGAGCGCGGCCTCAGCGCCATCCGCGAGCGCCTGGGGCCCCTGGTGGAACAGGGCCGCGTGCGGGCCGCCACTGTGGGCTCCCTGGCCGGCCAGCCGCTACAGGAGCGGGCCC
my $old = time();
my $seq="GCGCTGATGGACGAGACCATGAAGGAGTTGAAGGCCTACAAATCGGAACTGGAGGAACAACTGACCCCGGTGGCGGAGGAGACGCGGGCACGGCTGTCCAAGGAGCTGCAGGCGGCGCAGGCCCGGCTGGGCGCGGACATGGAGGACGTGTGCGGCCGCCTGGTGCAGTACCGCGGCGAGGTGCAGGCCATGCTCGGCCAGAGCACCGAGGAGCTGCGGGTGCGCCTCGCCTCCCACCTGCGCAAGCTGCGTAAGCGGCTCCTCCGCGATGCCGATGACCTGCAGAAGCGCCTGGCAGTGTACCAGGCCGGGGCCCGCGAGGGCGCCGAGCGCGGCCTCAGCGCCATCCGCGAGCGCCTGGGGCCCCTGGTGGAACAGGGCCGCGTGCGGGCCGCCACTGTGGGCTCCCTGGCCGGCCAGCCGCTACAGGAGCGGGCCC";

my $len=length($seq);
my $seed;
my $loc=45411791;
my %iterator;
#my %tmp;
for(my $i=1;$i<=($len-10+1);$i++){
	$seed= substr($seq,$i-1,10);
	$iterator{$seed}=$loc;
	$loc+=1;

	#if(exists $tmp{$seed}){
	#	$tmp{$seed}++;
	#}else{
	#	$tmp{$seed}=1;
	#}
}

my $seq2="GGGCCCGCTCCTGTAGCGGCTGGCCGGCCAGGGAGCCCACAGTGGCGGCCCGCACGCGGCCCTGTTCCACCAGGGGCCCCAGGCGCTCGCGGATGGCGCTGAGGCCGCGCTCGGCGCCCTCGCGGGCCCCGGCCTGGTACACTGCCAGGCGCTTCTGCAGGTCATCGGCATCGCGGAGGAGCCGCTTACGCAGCTTGCGCAGGTGGGAGGCGAGGCGCACCCGCAGCTCCTCGGTGCTCTGGCCGAGCATGGCCTGCACCTCGCCGCGGTACTGCACCAGGCGGCCGCACACGTCCTCCATGTCCGCGCCCAGCCGGGCCTGCGCCGCCTGCAGCTCCTTGGACAGCCGTGCCCGCGTCTCCTCCGCCACCGGGGTCAGTTGTTCCTCCAGTTCCGATTTGTAGGCCTTCAACTCCTTCATGGTCTCGTCCATCAGCGC";


my $len2=length($seq2);
my $seed2;
my $loc2=45412229;
my %iterator2;
for(my $i=1;$i<=($len2-10+1);$i++){
	$seed2= substr($seq2,$i-1,10);
	$iterator2{$seed2}=$loc2;
	$loc2=$loc2-1;

	#if(exists $tmp{$seed}){
	#	$tmp{$seed}++;
	#}else{
	#	$tmp{$seed}=1;
	#}
}

#foreach my $key (keys %tmp){
#	print "$key\t$tmp{$key}\n";
#}

my $inFile=$ARGV[0];
my $outdir=$ARGV[1];
my $name=$ARGV[2];
open OUT,">$outdir/$name";


my $line;
if( $inFile =~ /.*\.bz2/ ){ 
         open(IN, "bzcat $inFile |") or die "Can not open $inFile\n"; 
}elsif( $inFile =~ /.*\.gz/ ){ 
         open(IN, "gunzip -c $inFile |") or die "Can not open $inFile\n"; 
}else{
         open(IN, $inFile) or die "Can not open $inFile\n"; 
}

#my $i=1;
while($line=<IN>){
	#print "$i\n";
	chomp $line;
	my $fq2=<IN>; chomp $fq2;
	my $fq3=<IN>; chomp $fq3;
	my $fq4=<IN>; chomp $fq4;

	#CGCGGACATGGAGGACGTGT
	#GGAGGACGTGTGCGGCCGCCT
	#TGCGGCCGCCTGGTGCAGTA

	#TGCCGATGACCTGCAGAAGC
	#CCTGCAGAAGCGCCTGGCAGT
	#CGCCTGGCAGTGTACCAGGC
	if($fq2=~/CGCGGACATGGAGGACGTGT/ ||$fq2=~/GGAGGACGTGTGCGGCCGCCT/ ||$fq2=~/TGCGGCCGCCTGGTGCAGTA/ ||$fq2=~/TGCCGATGACCTGCAGAAGC/ ||$fq2=~/CCTGCAGAAGCGCCTGGCAGT/ ||$fq2=~/CGCCTGGCAGTGTACCAGGC/){
		 
		my $aa=system("perl smith-waterman_algorithm.pl $seq  $fq2 >$outdir/tmp");
		my $bb=system("perl smith-waterman_algorithm.pl $seq2 $fq2 >$outdir/tmp2");
		#print "$aa\n";
		
		if(-z "$outdir/tmp"){
			next;
		}else{
			open IN2,"$outdir/tmp";
			my $line2;
			while($line2=<IN2>){
				chomp $line2;
				if (length($line2)>0){
					my $ref=$line2;
					my $query=<IN2>;
					chomp $query;
					next if($ref=~/-/ || $query=~/-/);
					
					my $seedtmp= substr($ref,0,10);
					my $start=$iterator{$seedtmp};
					my $base1;
					my $base2;
					
					my $len2=length($query);
					for(my $i=0;$i<$len2;$i++){
						my $loctmp=$start+$i;
						my $base=uc(substr($query,$i,1));

						if($loctmp==45411941){
							$base1=$base;
						}
						if($loctmp==45412079){
							$base2=$base;
						}
					}
					
					#chr19   45411941        rs429358        T       C
					#chr19   45412079        rs7412          C       T

					#E1  rs429358-C rs7412-T,
					#E2  rs429358-T rs7412-T,
					#E3  rs429358-T rs7412-C, 
					#E4  rs429358-C rs7412-C,
					if ($base1 eq 'C' and $base2 eq 'T'){
						print OUT "$fq2\n";
						print OUT "APOE type: E1\n";
					}
					if ($base1 eq 'T' and $base2 eq 'T'){
						print OUT "$fq2\n";
						print OUT "APOE type: E2\n";
					}
					if ($base1 eq 'T' and $base2 eq 'C'){
						print OUT "$fq2\n";
						print OUT "APOE type: E3\n";
					}
					if ($base1 eq 'C' and $base2 eq 'C'){
						print OUT "$fq2\n";
						print OUT "APOE type: E4\n";
					}
					
				}
			}
			close IN2;
		}
		
		if(-z "$outdir/tmp2"){
			next;
		}else{
			open IN2,"$outdir/tmp2";
			my $line2;
			while($line2=<IN2>){
				chomp $line2;
				if (length($line2)>0){
					my $ref=$line2;
					my $query=<IN2>;
					chomp $query;
					next if($ref=~/-/ || $query=~/-/);
					
					my $seedtmp= substr($ref,0,10);
					my $start=$iterator2{$seedtmp};
					my $base1;
					my $base2;
					
					my $len2=length($query);
					for(my $i=0;$i<$len2;$i++){
						my $loctmp=$start-$i;
						my $base=uc(substr($query,$i,1));

						if($loctmp==45411941){
							$base1=$base;
						}
						if($loctmp==45412079){
							$base2=$base;
						}
					}
					
					#chr19   45411941        rs429358        T       C
					#chr19   45412079        rs7412          C       T

					#E1  rs429358-C rs7412-T,
					#E2  rs429358-T rs7412-T,
					#E3  rs429358-T rs7412-C, 
					#E4  rs429358-C rs7412-C,
					if ($base1 eq 'G' and $base2 eq 'A'){
						print OUT "$fq2\n";
						print OUT "APOE type: E1\n";
					}
					if ($base1 eq 'A' and $base2 eq 'A'){
						print OUT "$fq2\n";
						print OUT "APOE type: E2\n";
					}
					if ($base1 eq 'A' and $base2 eq 'G'){
						print OUT "$fq2\n";
						print OUT "APOE type: E3\n";
					}
					if ($base1 eq 'G' and $base2 eq 'G'){
						print OUT "$fq2\n";
						print OUT "APOE type: E4\n";
					}
					
				}
			}
			close IN2;
		}
	}
	#$i++;

}
close IN;
close OUT;

my $new = time();
my $delta = $new - $old;
print "Run time：$delta seconds\n";
