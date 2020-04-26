#!/usr/bin/env perl

# core
use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use Log::Dispatch;
use Math::Round;
use Path::Tiny;
use Try::Tiny;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

our $dir_out = path('.');

# script usage
my $progname = path($0)->basename;
my $USAGE = <<"__USAGE__";
usage: $progname <cath_version> <domain_list>

Generates a feature file for a list of CATH domains for machine learning applications

required:

  <cath_version>      Version of CATH (eg 4.2/4.1)
  <DOMAIN list>       File path for the list of domains with FunFam ID in the format 1.10.8.10-ff-1234 in tsv format
  <result-dir>        Name of the result directory
  <feature_outfile>   Name of the feature file
  
options:

  --dir_out  <dir>    output directory
                      [$dir_out]
  -v|--verbose        More verbose logs

__USAGE__



# get input data
my ($cath_db_version, $DOMAINLIST, $DIR, $OUTFILE) = @ARGV;
  
# exit script unless input data is provided
if( scalar @ARGV != 4) {
    print $USAGE;
    exit;
}

my $SCRIPTS_DIR = path("$DIR/scripts");
my $FEATURE_FILE_DIR = path("$DIR/feature_files");

my $LIST = path("$DOMAINLIST");

my @lines = $LIST->lines;

my $dom_group=0;

foreach my $entry (@lines){
    
    unless( $entry =~ /^\#/){
        
        my ($DOMAIN, $SF, $FF, $DOPS) =split("\t", $entry);
        my $ff_arg = "$SF"."-ff-"."$FF";
        $dom_group++;
        
        my $DOMAIN_FEATURE_FILE = path("$FEATURE_FILE_DIR/$DOMAIN.features.csv");
        
        if( !-e "$DOMAIN_FEATURE_FILE" || -z "$DOMAIN_FEATURE_FILE"){
    
            print "Generating features for $DOMAIN, $ff_arg\n";
            system("perl $SCRIPTS_DIR/funsite-feature-generation.pl $cath_db_version $DOMAIN $ff_arg $DIR $DOMAIN_FEATURE_FILE $dom_group");
            
            system("cat $DOMAIN_FEATURE_FILE >> $OUTFILE");
            exit 0;
            
        }
    }
}


