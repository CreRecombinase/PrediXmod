#!/usr/bin/perl
use strict;
use warnings;
                                        #By NWK 20140202####
##Compute Expected dosage from impute2 data (also generates annotation RDS and zipped txt files on a per-chromosome level##

my $Infile = $ARGV[0];
my $chr = $ARGV[1];
my $Annofile = $ARGV[2];
my $Dosagefile = $ARGV[3];

open(IF,"<",$Infile);
open(AF,">",$Annofile);
open(DF,">",$Dosagefile);

while(<IF>){
    chomp;
    my ($tl,$rsid,$pos,$ref,$alt,@genos) = split(/\s/);
    my $al = join("\t",$tl,$rsid,$pos,$ref,$alt,$chr);
    print AF ($al."\n");
    while(my @three = splice @genos,0,3){
	my $dos = $three[1]+2*$three[2];
	print DF $dos."\t";
    }
    print DF "\n";
}
close(IF);
close(AF);
close(DF);
    
    
