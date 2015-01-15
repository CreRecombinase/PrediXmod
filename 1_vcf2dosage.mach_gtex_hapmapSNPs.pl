#!/usr/bin/perl
use strict;
use warnings;

####This perl script takes the GTEx imputed vcf file as input, removes ambiguous-strand SNPs (A/T and C/G)
#### and makes several output files:
#### .dosage.gz file for PrediXcan (compute.scores.py)
#### .mlinfo.gz and .mldose.gz MACH files for GCTA
#### .SNPxID matrix for quick scanning into R
#### .ID.list colnames of matrix
#### .SNP.list rownames of matrix
#### .bim plink bim file with SNP pos info



if( $#ARGV +1 !=4){
  die "\n Usage: GTEx_imputed_genotype_file.vcf.gz hapmapSNPsCEU.list output_directory projectname";
}


my $file = $ARGV[0]; #"GTEx_Analysis_2014-06-13_OMNI_2.5M_5M_451Indiv_allchr_genot_imput_info04_maf01_HEW1E6.vcf.gz";
my $snpfile = $ARGV[1];
my $outdir = $ARGV[2];
my $projectname = $ARGV[3];


open(VCF, "bgzip -cd ${file}|");
open(HAP, $snpfile); ##from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/hapmapSnpsCEU.txt.gz (only rsid column)

$snpfile =~ s/.+\/([^\/]+)\..+/\1/ 
my %hapmapsnps;
while(<HAP>){
    chomp;
    my $snp = (split '\t', $_)[4];
    $hapmapsnps{$snp} = 1;
}

#outfiles for downstream applications:
open(DOS, ">${outdir}/${projectname}.unamb.dosage");
open(ID, ">${outdir}/${projectname}.SNP.ID.list");
open(SNP, ">${outdir}/${projectname}.SNP.list");
open(SNPxID, ">${outdir}/${projectname}.SNPxID");
open(BIM, ">${outdir}/${projectname}.bim");
open(INTRO, ">${outdir}/intro");
open(MLINFO, ">${outdir}/${projectname}.mlinfo");

while(<VCF>){
  chomp;
  my ($first) = split(/\t/);
  if($first =~ m/##/){
    next;
  }
  my ($chr, $pos, $rs, $ref, $alt, $qual, $filter, $info, $format, @genos) = split(/\t/);
  my ($expfreq, $impinfo, $cert) = split(/;/,$info);
  my ($a, $freqalt) = split(/=/,$expfreq);
  my ($freqref) = 1 - $freqalt;
  my ($b, $quality) = split(/=/,$cert);
  my ($c, $rsq) = split(/=/,$impinfo);
  if($chr eq "#CHROM"){
    foreach my $id (@genos){
      my ($d,$e) = split(/-/,$id);
      my $shortid = $d . '-' .$e;	    
      print ID "$shortid\n";
      print INTRO "$shortid->$shortid MLDOSE\n";
    }
  }
  if($ref eq "A" && $alt eq "T"){ ##rm potentially ambiguous strand SNPs
    next;
  }elsif($ref eq "T" && $alt eq "A"){
    next;
  }elsif($ref eq "C" && $alt eq "G"){
    next;
  }elsif($ref eq "G" && $alt eq "C"){
    next;
  }elsif(defined($hapmapsnps{$rs}) && $pos =~ m/\d+/){ ###only pull rsid SNPs at this time & don't print header rows
    print DOS "$chr\t$rs\t$pos\t$ref\t$alt\t";
    print BIM "$chr\t$rs\t0\t$pos\t$ref\t$alt\n";
    print SNP "$rs\n";
    print MLINFO "$rs\t$ref\t$alt\t$freqref\t$quality\t$rsq\n";
    foreach my $geno (@genos){
      my ($probs, $gt, $dos) = split(/:/,$geno);
      print DOS "$dos\t";
      print SNPxID "$dos\t";
    }
    print DOS "\n";
    print SNPxID "\n";
  }
}




open(R, ">runR.R") or die "cant make runR.R\n";
print R "dat<-read.table(\"${outdir}/${projectname}.${snpfile}.SNPxID\")\n";
print R "dat<-t(dat)\n";
print R "write.table(dat,\"${outdir}/t.dos\",col.names=F,row.names=F,quote=F)\n";
close(R);
system("R --vanilla < runR.R");
system("paste -d\' \' ${outdir}/intro ${outdir}/t.dos > ${outdir}/${projectname}.${snpfile}.mldose");


system("gzip ${outdir}/${projectname}.${snpfile}.unamb.dosage");
system("gzip ${outdir}/${projectname}.${snpfile}.mldose");
system("gzip ${outdir}/${projectname}.${snpfile}.mlinfo");
system("rm ${outdir}/intro ${outdir}/t.dos");
