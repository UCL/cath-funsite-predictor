#!/usr/bin/env perl

# core
use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use Path::Tiny;
use Try::Tiny;
use List::MoreUtils qw( zip );

# script usage
my $progname = path($0)->basename;
my $USAGE = <<"__USAGE__";
usage: $progname <cath_version> <domain_list>

example:

  $progname 4.2 /cath/people2/ucbtdas/domains.list

__USAGE__

# get input data
my ($cath_db_version, $DOMAIN, $funfam_scorecons_map) = @ARGV;

 #exit script unless input data is provided
if( scalar @ARGV != 2) {
    print $USAGE;
    exit;
}

#my $cath_db_version=4.2;

my $DIR_SC5=path("/export/sayoni/sc5_data_v4_2_0");
my $DIR_SCONS_MAP= path("/cath/people2/ucbtdas/git/cath-funsite/funsites-ml/generate_domain_features/results/scons_map");
my $SCRIPTS_DIR = path("/cath/people2/ucbtdas/git/cath-funsite/funsites-ml/generate_domain_features/scripts");
my $OUTDIR=path("/cath/people2/ucbtdas/git/cath-funsite/funsites-ml/generate_domain_features/results/sc5_map");

        if(-e "$funfam_scorecons_map"){

        my @scons_map_domain=$funfam_scorecons_map->lines;

        my $get_sc5=`fgrep "$DOMAIN" $DIR_SC5/domain_map`;
        chomp($get_sc5);

        if($get_sc5){

            my ($domain, $sc5_name)=split("\t",$get_sc5);
            my $SC5FILE= $DIR_SC5->path("$sc5_name.cathdom_sdp_scons.txt");

            my $sc5_groupsim_line = `fgrep "groupsim" $SC5FILE`;
            my $sc5_scons_line = `fgrep "scorecons" $SC5FILE`;
            my $sc5_domain_line = `fgrep "$DOMAIN" $SC5FILE`;
            chomp($sc5_domain_line);
            chomp($sc5_scons_line);
            chomp($sc5_groupsim_line);

            my @sc5_dom;
            my @sc5_gs;
            my @sc5_scons;

            if($sc5_domain_line & $sc5_groupsim_line){
                my ($sc5_domain_seq, $id)=split(" ",$sc5_domain_line);
                my ($sc5_domain_gs, $id2)=split(" ",$sc5_groupsim_line);
                my ($sc5_domain_scons, $id3)=split(" ",$sc5_scons_line);
                $sc5_domain_seq= uc($sc5_domain_seq);
                @sc5_dom = split("", $sc5_domain_seq);
                @sc5_gs= split("", $sc5_domain_gs);
                @sc5_scons= split("", $sc5_domain_scons);
            }

            my @sc5_dom_aa_gs;
            my @sc5_dom_aa_scons;
            my $index=0;
            my $map_count=1;
            foreach my $i (@sc5_dom){
                if($i ne "-"){
                    #print "$DOMAIN\t$i\t$sc5_gs[$index]\n";
                    push(@sc5_dom_aa_gs, $sc5_gs[$index]);
                    push(@sc5_dom_aa_scons, $sc5_scons[$index]);
                }
                $index++;
            }

            my $filename=$OUTDIR->path("$DOMAIN.sc5_map");
            open(SC5_MAP, '>', $filename) or die "Could not open file '$filename' $!";


            my $count=0;
            for my $line (@scons_map_domain){
                chomp($line);
                if($line=~/^\d/){
                    #print("$line\n");
                    my @tab=split("\t", $line);
                    my $pdb_res= $tab[6];
                    my $dom_aa = $tab[5];
                    if($dom_aa ne "-"){

                        my $gs_value=$sc5_dom_aa_gs[$count];
                        my $scons_value=$sc5_dom_aa_scons[$count];

                        if($gs_value && $scons_value && $pdb_res ne "-"){
                        #print($gs_value);
                        if($scons_value eq "-"){
                            $scons_value=-1;
                        }
                        if($gs_value eq "-"){
                            $gs_value=-1;
                        }
                        #print "$line\t$gs_value\t$scons_value\n";
                        print SC5_MAP "$line\t$gs_value\t$scons_value\n";
                    }
                        $count++;
                    }
                }

            }

            close(SC5_MAP);

        }

    }
