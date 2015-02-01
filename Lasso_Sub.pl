#!/usr/bin/perl
use strict;
use warnings;



my $tissuefile = $ARGV[0];
my $SCRIPTDIR = $ARGV[1];

open(TISFILE,$tissuefile);


while(<TISFILE>){
    chomp;
    my $line = $_;
    (my $tis = $_) =~ s/\s//g;
    my $dir = "/group/dolan-lab/nknoblauch/".$tis;
    if( -e $dir){
	opendir(DIR,$dir);
	my @files = grep(/.+EXPanno\.[0-9]+\.[0-9]+.RDS/,readdir(DIR));
	closedir(DIR);
	my $size = @files;
	if ($size >=1){
	    foreach my $file (@files){
		(my $chr= $file) =~ s/.+EXPanno\.([0-9]+)\.([0-9]+)\.RDS/$1/;
		(my $seg = $file)  =~ s/.+EXPanno\.([0-9]+)\.([0-9]+)\.RDS/$2/;
		my $filename = "/group/dolan-lab/nknoblauch/${tis}/GTEx_PrediXmod.${tis}.ResultsArray.${chr}.${seg}.RDS";
		$tis =~ s/\(|\)//g;
		my $shellfile = "${SCRIPTDIR}/LASSO.${tis}.${chr}.${seg}.sh";
		unless (-e $shellfile){
		    open(DOSCRIPT,">${shellfile}");
		    print DOSCRIPT "#!/bin/bash\n";
		    print DOSCRIPT "#PBS -N LASSO.${tis}.${chr}.${seg}\n";
		    print DOSCRIPT "#PBS -S /bin/bash\n";
		    print DOSCRIPT "#PBS -l mem=3gb\n";
		    print DOSCRIPT "#PBS -l walltime=36:00:00\n";
		    print DOSCRIPT "#PBS -o ${SCRIPTDIR}/LASSO.${tis}.${chr}.${seg}.out\n";
		    print DOSCRIPT "#PBS -e ${SCRIPTDIR}/LASSO.${tis}.${chr}.${seg}.err\n";
		    print DOSCRIPT "module load R/3.1.0\n";
		    print DOSCRIPT "Rscript /home/t.cri.nknoblauch/PrediXmod/4_CV_GTEx_lasso_adjusted.R GTEx_PrediXmod /group/dolan-lab/nknoblauch  \'${line}\' ${chr} ${seg}\n";
		    close(DOSCRIPT);
		    system("qsub ${SCRIPTDIR}/LASSO.${tis}.${chr}.${seg}.sh");
		    print "$line ${tis}.${chr}.${seg}\n";
		    sleep(1);
		}
	    }
	}
    } else{
	print("No such Directory ${dir}\n");
    }
}

