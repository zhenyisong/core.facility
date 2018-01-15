#!/usr/bin/perl

# @author Yisog Zhen
# @since  2018-01-11
# @update 

#---
# code root
# /wa/zhenyisong/sourcecode/core.facility/qualiControl/QCparser.pl
#---

use Modern::Perl;
use autodie;
use FileHandle;
use File::Temp qw(tempfile);
use Switch;
use File::Basename;
use Cwd;

# I use anaconda to manager the Perl
# root env
# cpanm Modern::Perl
#--

# http://stackoverflow.com/questions/84932/how-do-i-get-the-full-path-to-a-perl-script-that-is-executing
# this method is the better way.
#---
# my $WORKING_DIR = dirname(__FILE__);
# the above discarded! 
# 
my $WORKING_DIR = getcwd();
chdir $WORKING_DIR;

my $picard_path = '/home/zhenyisong/data/results/chenlab/xiaoning/qc.step/hisat2';
chdir $picard_path;
my @picard_files = getAllPicardFiles();

for my $i (0..$#picard_files) {
    
    my $picard_dict_ref = getPicardMetrics($picard_files[$i]);
    my ($base_name) = ($picard_files[$i] =~ /(\w+)\.picard/);
    print $base_name, "\t", $picard_dict_ref->{'PCT_RIBOSOMAL_BASES'},"\n";
}




print "we complete the parsing procedure!\n";




#---
# @para
#    filename: the input is the picard file path (directories)
#              name which saves the 
#              picard  results
# @return
#   the hash reference which contains the picard file statistic
#---
sub getAllPicardFiles {
    my @picardFiles     = `ls *.picard`;
    return @picardFiles;
}

# @parameters
#    the file name from the picard analysis program
# @return
# the dictionary for all base statistic
# the specified column, in this case, is the PCT_RibosomeRNA
#----

sub getPicardMetrics {
    my $fileName = shift;
    my $fh       = FileHandle->new($fileName);
    my @metricsHeader = ();
    my @metricsNumber = ();
    my $hashRef       = undef;
    while(my $line = $fh->getline()) {
        chomp $line;
        if( $line =~ /^PF_BASES/) {
            
            @metricsHeader = ($line =~ /\S+/g);
            $line = $fh->getline();
            chomp $line;
            @metricsNumber = ($line =~ /\S+/g);
            last;
        }
    }
    #print @metricsHeader,"\n";
    # my $result = $metricsNumber[15];
    # initilize the hash ref using for-loop
    for my $i (0..$#metricsHeader) {
        $hashRef->{$metricsHeader[$i]} = $metricsNumber[$i];
    }
    $fh->close;
    return  $hashRef;
    
}

