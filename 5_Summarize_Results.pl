#!/usr/bin/perl
use strict;
use warnings;

my $Projectname = $ARGV[0];
my $Basedir = $ARGV[1];
my $Outputdir = $ARGV[2];

print $Projectname."\n".$Basedir."\n".$Outputdir."\n";
opendir my $dh,$Basedir
    or die "$0: opendir: $!";

my @dirs = grep {-d "$Basedir/$_" && !/^\.{1,2}$/} readdir($dh);

foreach my $dir (@dirs){
    print $dir."\n";
    opendir(my $tdh,$Basedir."/".$dir);
    my @allfiles = readdir($tdh);
    my @annofiles = grep(/EXPanno/,@allfiles);
    my @finishedfiles = grep(/Finished/,@allfiles);
    my @betafiles = grep (/BetaTable/,@allfiles);
    my @resultfiles = grep (/ResultsArray/,@allfiles);
    print($#annofiles." annofiles and ".$#finishedfiles." finished files in ${dir}\n");
    if($#annofiles == $#finishedfiles and $#betafiles == $#resultfiles){
	if($#finishedfiles>0){
	    my $finishedBetafile = "${Outputdir}/${dir}.allBetas.txt";
	    my $finishedResultsfile = "${Outputdir}/${dir}.allResults.txt";
	    unless( -e $finishedBetafile or -e $finishedResultsfile){
		open(BB,">",$finishedBetafile) or die "Can't open BB\n";
		open(BR,">",$finishedResultsfile) or die "Can't open BR\n";
		my $fbetafile = $betafiles[0];
		my $fresultfile = $resultfiles[0];
		open(TB,"<","${Basedir}/${dir}/${fbetafile}") or die "Can't open TB ${fbetafile}\n";
		open(TR,"<","${Basedir}/${dir}/${fresultfile}") or die "Can't open TR\n";	
		my $tline = <TB>;
		print BB $tline;
		close TB or die "Can't close: ${fbetafile}!";
		$tline = <TR>;
		print BR $tline;
		close(TR);
		for my $i (0 .. $#betafiles){
		    unless( -e "${Basedir}/${dir}/".$betafiles[$i]){
			die "No such file {Basedir}/${dir}/".$betafiles[$i];
		    }
		    open(TB,"<","${Basedir}/${dir}/".$betafiles[$i]) or die "Can't open TB :".$betafiles[$i];
		    my $stline = <TB>;
		    while (<TB>){
			print BB $_;
		    }
		    close(TB);
		    unless( -e "${Basedir}/${dir}/".$resultfiles[$i]){
			die "No such file {Basedir}/${dir}/".$resultfiles[$i];
		    }
		    open(TR,"<","${Basedir}/${dir}/".$resultfiles[$i]) or die "Can't open TR :".$resultfiles[$i];
		    my $ttline = <TR>;
		    while (<TR>){
			print BR $_;
		    }
		    close(TR);
		}
		close(BB);
		close(BR);
	    }
	    print $dir."\n";
	}
	else{
	    print ("${dir} has no finished files\n");
	}
	
    }
    else{
	print ("${dir} has ".$#annofiles." annofiles, ".$#finishedfiles." finished files, ".$#resultfiles." resultfiles and ".$#betafiles." betafiles\n");
    }   
}


    
