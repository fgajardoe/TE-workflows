#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper qw(Dumper);
use Getopt::Std;

our($opt_h,$opt_n,$opt_i,$opt_g);

# Usage
my $usage="
Usage:

   perl RepeatMaskerAln2KimuraBed.pl [-ngh] -i RepeatMaskerAln.cat[.gz] > RepeatsKimura.bed

Options:
   -i RepeatMaskerAln.cat[.gz]		RepeatMasker alignment file (as plain text).
   -g					Boolean. If present assumes compressed(.gz) input file.

   -n sep1, sep2,...,sepN		Split name by sep1 OR sep2 OR ... OR sepN. Useful when
   					repeat name contains TE family and order, and you want
					have them on separated columns.
        				Example: rnd-5_family-81#LTR/Ngaro 

   -h					Print this useful help.


Bye bye.

";

# Get options
getopts('i:hn:g');

# Check if there is input
if(($opt_h)||(!$opt_i)){
	print $usage;
	exit 0;
}



my $command;
if($opt_g){
        $command="zcat $opt_i |";
}
else{
        $command=$opt_i;
}

open(ALN, $command) or die $!;

my @lines=<ALN>;

my @content=();
my $defline="";
my %hash=();

my %kimuras=();
my $flag=0;

# data assigment
foreach my $line(@lines){
	chomp($line);
	if($line eq ""){
		next;
	}
	else{
		
		if($line=~m/^\d/){
			if(!exists $hash{$defline}){
				
				$hash{$defline}=1;

			}
			$defline=$line;
			$flag=0;
		}
		else{
			if($line=~/^Kimura/){

				my @tmp=split("=",$line);
				$tmp[1]=~s/^ //g;
				my $kimura=$tmp[1];
				$kimuras{$defline}=$kimura;
				$flag=1;
			}
			else{
				if($flag==0){
					$kimuras{$defline}="ND";
				}
			}
		}
	}
}

if(!exists $hash{$defline}){
	$hash{$defline}=1;
}


# print
foreach my $key(keys %hash){
	if($key ne ""){
		my @fields=split(" ",$key);
		my $scaff=$fields[4];

		my $start=$fields[5]-1;

		my $end=$fields[6];
		my $score=$fields[0];

		my $strand;
		my $name;
	
		if($fields[8] eq "C"){
			$strand="-";
			$name=$fields[9];

		}
		else{
			$strand="+";
			$name=$fields[8];
		}
		
		my $att;
		if($opt_n){
			                                     
			$opt_n=~s/,/\|/g;
			my @name_splitted=split(/$opt_n/,$name);
			$name=$name_splitted[0];
			shift(@name_splitted);
			$att=join("\t",@name_splitted);
		}

		my $kimura=$kimuras{$key};

		if($opt_n){
			print join("\t",$scaff,$start,$end,$name,$score,$strand,$kimura,$att)."\n";
		}
		else{
			print join("\t",$scaff,$start,$end,$name,$score,$strand,$kimura)."\n";
		}

	}
}

