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

# local
use Path::Class;
use Cath::Schema::Biomap;
use Cath::Util;
use Cath::Data;
use Cath::SSGID;
use Cath::Stockholm;
use Bio::Structure::SecStr::DSSP::Res;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

our $dir_out = path('.');

# script usage
my $progname = path($0)->basename;
my $USAGE = <<"__USAGE__";
usage: $progname <cath_version> <domain_list>

Generates a features for one CATH domain for machine learning applications
  
options:

  --dir_out  <dir>    output directory
                      [$dir_out]
  -v|--verbose        More verbose logs
 

__USAGE__



# get input data
my ($cath_version, $DOMAIN, $ff_arg, $dir, $domain_feature_file, $dom_group) = @ARGV;
  
# exit script unless input data is provided
if( scalar @ARGV != 6) {
    print $USAGE;
    exit;
}

# set up logging
my $min_log_level = 'info';
my $LOG = Log::Dispatch->new(
  outputs => [
    [ 'Screen', min_level => $min_log_level, stderr => 0, newline => 0 ],
    [ 'Screen', min_level => 'warning', stderr => 1, newline => 0 ],
  ],
  callbacks => sub {
    my %p = @_;
    my @lines = split (/\n/, $p{message});
    my $msg = join( "", map { sprintf( "%s [%s] %s\n", DateTime->now()->datetime(), $p{level}, $_ || '' ) } @lines );
    return $msg;
  },
);

            
my $cath_db_version = try { Cath::Util::get_cath_version( $cath_version ) }
catch {
    $LOG->error( "Error: failed to parse CATH version '$cath_version': $_" );
    exit;
};


my $funfam_id = try { Cath::SSGID->new( id => $ff_arg ) }
catch {
    $LOG->error( "Error: failed to parse FunFam ID from string '$ff_arg': $_" );
    exit;
};

my $APPS_DIR = path("$dir/apps");
my $FEATURE_DIR = path("$dir/results");
my $DOMPDB = path("$dir/data/pdbs/$DOMAIN.pdb");
my $DOM_FEATURE_FILE = path("$domain_feature_file");
my $FH = $DOM_FEATURE_FILE->openw;

#################################################################################
# Run and load FunFam Alignment features incl Scorecons
#################################################################################

# define cath data directory
my $funfam_id_str = $funfam_id->to_string( join_char => '-' ); 
my $funfam_aln_file = Cath::Data->get_data_file( 
    file_type => 'FunfamStockholm', 
    id => $funfam_id, 
    version => $cath_db_version 
)->to_string;

my $aln = Cath::Stockholm->parse_file( $funfam_aln_file );
    
our %dom_feature_data=();
    
my $aln_dops = $aln->dops; #*#
    
# run the scorecons mapper on the $funfam_aln_file
my $funfam_scorecons_map = $dir_out->path("$ff_arg-scorecons.map");

unless(-e "$funfam_scorecons_map"){
    
    system("/usr/local/svn/source/cpan/trunk/Cath/script/cath-funfam-scorecons-mapper.pl -c $cath_db_version -f $funfam_id -r $DOMAIN -o $funfam_scorecons_map");

}

my @scons_lines = $funfam_scorecons_map->lines;
my $header = shift @scons_lines;
    
our %aln_pos_map=();
my $domain_length = 0;

foreach my $line (@scons_lines){
		
    chomp($line);
	my @array = split("\t", $line);
	my $aln_pos = $array[1];
	my $scons =$array[2];
	my $residue_aa = $array[5];
	my $residue_num = $array[6];
	
	unless($residue_num eq "-" || $residue_aa eq "X"){
            

		$aln_pos_map{$aln_pos}=$residue_num;
        
		$dom_feature_data{scons}{$residue_num} = $scons;            	#*#
		$dom_feature_data{dops}{$residue_num} = $aln_dops;            	#*#
		$dom_feature_data{residue_aa}{$residue_num} = $residue_aa;  	#*#
		$dom_feature_data{dom_group}{$residue_num} = $dom_group;    	#*#
		$domain_length++;
	}
	
}

#################################################################################
# Load CATH, FOLD INFO
#################################################################################

my $SF = $funfam_id->superfamily_id;
$SF=~ /^(\w\.\w+\.\w+)\.\w+/;
my $FOLD = $1;
$FOLD=~ /^(\w\.\w+)\.\w+/;
my $CAth = $1;

#print "domain: $DOMAIN, FF: $ff_arg, SF: $SF, FOLD: $FOLD\n\n";
$DOMAIN=~ /^(\w\w\w\w)(\w)/;
my $PDB = $1;
my $CHAIN = $2;

# connect to the database
my $db = Cath::Schema::Biomap->connect_by_version("$cath_version");

my $pdbchain_domain_count=0;
# search for entries that map to superfamily '1.10.8.10'
my $rs = $db->resultset("Domain")->search( { pdb_code => "$PDB" , chain_code => "$CHAIN" } );
# go through the results row by row
while ( my $row = $rs->next ) {
         my $dom = $row->domain_id;
        $pdbchain_domain_count++;
}

foreach my $line (@scons_lines){
	chomp($line);
	my @array = split("\t", $line);
	my $residue_aa = $array[5];	
	my $residue_num = $array[6];
	unless($residue_num eq "-" || $residue_aa eq "X"){
		$dom_feature_data{domain_length}{$residue_num} = $domain_length;   	                #*#
		$dom_feature_data{domain_num_in_chain}{$residue_num} = $pdbchain_domain_count;  	#*#
        $dom_feature_data{fold}{$residue_num} = $FOLD;   	                                #*#
        $dom_feature_data{CAth}{$residue_num} = $CAth;   	                            #*#
        $dom_feature_data{SF}{$residue_num} = $SF;
	}
}

#################################################################################
# Load PSSM INFO
#################################################################################

my $pssm_file = path("$FEATURE_DIR/pssm/$DOMAIN.pssm");
my $OUTDIR_PSSMFEATURE=path("/cath/people2/ucbtdas/git/cath-funsite/funsites-ml/generate_domain_features/feature_files/pssm_features");
my $SCRIPTS_DIR = path("/cath/people2/ucbtdas/git/cath-funsite/funsites-ml/generate_domain_features/scripts");

system("python $SCRIPTS_DIR/generate_pssm_features_cath_domain.py $DOMAIN $pssm_file $funfam_scorecons_map > $OUTDIR_PSSMFEATURE/$DOMAIN.pssmfeatures");

#################################################################################
# Load SC5 INFO
#################################################################################

my $sc5_file = path("$FEATURE_DIR/sc5_map/$DOMAIN.sc5_map");
system("perl $SCRIPTS_DIR/map_sc5_sdp_scons_to_domain.pl 4.2 $DOMAIN $funfam_scorecons_map");

#################################################################################
# Load Naccess RSA & Distance of residue to Center of Mass
#################################################################################

my $naccess_file = path("$FEATURE_DIR/rsa/$DOMAIN.rsa");
my $CM_file = path("$FEATURE_DIR/CM/$DOMAIN.center_of_mass");
&load_naccess_cm_features($naccess_file, $APPS_DIR, $DOMPDB, $CM_file);



#################################################################################
# Load Speedfill cleft features and calculate distance to center of mass
#################################################################################

my $cleft_file = path("$FEATURE_DIR/speedfill/$DOMAIN.clefts.depth");
&load_cleft_features($cleft_file);




#################################################################################
# Load AA Indices -> map to {residue_aa}
#################################################################################

our %aa_indices=();
my $AAindices_file = path("$APPS_DIR/aa_indices/AAindex.tsv");
my @aa_index_names = &load_aaindices($AAindices_file);


for my $residue_num ( sort keys %{ $dom_feature_data{residue_aa} } ) {

    my $amino_acid = $dom_feature_data{residue_aa}{$residue_num};
    foreach my $aa_index (@aa_index_names){
        
        $dom_feature_data{$aa_index}{$residue_num} = $aa_indices{$amino_acid}{$aa_index};    #*#
        
    }
}




#################################################################################
# Load distmap
#################################################################################

our %distmap=();
my $DISTMAP = path("$FEATURE_DIR/distmap/$DOMAIN.distmap");
&load_distmap($DISTMAP);




#################################################################################
# Load min. distance to clefts [123] and surface residues [25% RSA cutoff]
#################################################################################

&load_min_distance_to_cleft();
&load_min_distance_to_surface();




#################################################################################
# Load structural neighborhood
#################################################################################

# Calculate average RSA, Scorecons, no. of pockets residue in the structural neighbourhood
# Calculate average AA properties of residues in the structural neighbourhood

my $distance_cutoff=5;
&load_structural_neighborhood_properties($distance_cutoff);





#################################################################################
# Load DSSP features
#################################################################################

my $dssp_file = path("$FEATURE_DIR/dssp/$DOMAIN.dssp");
&load_dssp_features($dssp_file);




#################################################################################
# Load FoldX features
#################################################################################

my $foldx_filename = "$FEATURE_DIR/foldx/$DOMAIN"."_AS.fxout";
my $foldx_file = path("$foldx_filename");
&load_foldx_features($foldx_file);




#################################################################################
# Load PSAIA features
#################################################################################

my $psaia_file = path("$FEATURE_DIR/psaia/$DOMAIN.psaia");
&load_psaia_features($psaia_file);




#################################################################################
# Load centrality measures
#################################################################################

my $betweenness_centrality = path("$FEATURE_DIR/centrality/$DOMAIN.pdb.betweenness");
&load_centrality_feature("betweenness", $betweenness_centrality);

my $closeness_centrality = path("$FEATURE_DIR/centrality/$DOMAIN.pdb.closeness");
&load_centrality_feature("closeness", $closeness_centrality);

my $degree_centrality = path("$FEATURE_DIR/centrality/$DOMAIN.pdb.degree");
&load_centrality_feature("degree", $degree_centrality);




#################################################################################
# Load normalised B-factor
#################################################################################

&load_norm_bfactors($DOMPDB);





#################################################################################
# Load IBIS LIG/PPI freq scores - > map to {ff_aln_pos}
#################################################################################

my $ibis_freq_score_file = path("/export/sayoni/funfam_ibis_freq_scores/4.2/$ff_arg.ibis.anno.freq");
&load_ff_ibis_freq($ibis_freq_score_file);




#################################################################################
# Load Annotations
#################################################################################

&load_annotations($DOMAIN);

# Load MDA (SFNUM AND DOMAIN NUM IN A CHAIN)#????

#print Dumper(\%dom_feature_data);


#################################################################################
# Print feature names and features and sliding features +-3 residues
#################################################################################

if ($dom_group != 0) {
    
    my @feature_names;
    
    push(@feature_names, "residue_string");
    
    for my $feature_name ( sort keys %dom_feature_data ) {
    
        if ($feature_name =~ /^annotation/ || $feature_name =~ /residue_aa/) {
            push(@feature_names,$feature_name);
        }
        else{
            push(@feature_names,$feature_name);
         
            #my $feature_before1_string = "$feature_name"."_before1";
            #my $feature_before2_string = "$feature_name"."_before2";
            #my $feature_after1_string = "$feature_name"."_after1";
            #my $feature_after2_string = "$feature_name"."_after2";
            
            #push(@feature_names,$feature_before1_string);
            #push(@feature_names,$feature_before2_string);
            #push(@feature_names,$feature_after1_string);
            #push(@feature_names,$feature_after2_string);
    
        }
           
    }

    foreach my $feature (@feature_names){
        print $FH "$feature,";
    }
    print $FH "\n";
}

for my $residue_num ( sort keys %{ $dom_feature_data{scons} } ) {
        
        my $residue_string = "$DOMAIN"."_"."$residue_num";
        print $FH "$residue_string,";
        
        for my $feature_name ( sort keys %dom_feature_data ) {
    
            if ($feature_name =~ /^annotation/ || $feature_name =~ /residue_aa/) {
                
                # print annotations
                my $feature_value = 0;
                
                if ($dom_feature_data{$feature_name}{$residue_num}) {
                    $feature_value =  $dom_feature_data{$feature_name}{$residue_num};
                    print $FH "$feature_value,";
                }
                else{
                    print $FH "$feature_value,";
                }
                
            }
            else{
                
                # print non-annotations with sliding features (first feature, then its before1,before2, then after1,after2)
                
                my $feature_value = "";
                if (defined($dom_feature_data{$feature_name}{$residue_num})) {
                    my $feature_value =  $dom_feature_data{$feature_name}{$residue_num};
                    print $FH "$feature_value,";
                }
                else{
                    print $FH "$feature_value,";
                }
                
=head
                $feature_value = "";
                my $res_before1 = $residue_num - 1;
                if ($dom_feature_data{$feature_name}{$res_before1}) {
                    $feature_value =  $dom_feature_data{$feature_name}{$res_before1};
                    print $FH "$feature_value,";
                }
                else{
                    print $FH "$feature_value,";
                }


                $feature_value = "";
                my $res_before2 = $residue_num - 2;
                if ($dom_feature_data{$feature_name}{$res_before2}) {
                    $feature_value =  $dom_feature_data{$feature_name}{$res_before2};
                    print $FH "$feature_value,";
                }
                else{
                    print $FH "$feature_value,";
                }           
            
            
                $feature_value = "";
                my $res_after1 = $residue_num + 1 ;
                if ($dom_feature_data{$feature_name}{$res_after1}) {
                    $feature_value =  $dom_feature_data{$feature_name}{$res_after1};
                    print $FH "$feature_value,";
                }
                else{
                    print $FH "$feature_value,";
                }
                
             
                $feature_value = "";
                my $res_after2 = $residue_num + 2;
                if ($dom_feature_data{$feature_name}{$res_after2}) {
                    $feature_value =  $dom_feature_data{$feature_name}{$res_after2};
                    print $FH "$feature_value,";
                }
                else{
                    print $FH "$feature_value,";
                }
=cut

            }
            
        }
        print $FH "\n";
        #exit 0;

}

unlink("$funfam_scorecons_map");


###########################
# SUBROUTINES
###########################

sub load_naccess_cm_features{
    

    my ($naccess_file, $APPS_DIR, $DOMPDB, $CM_file) = @_;
    
    if (!-e "$naccess_file" || -z "$naccess_file") {
        # generate naccess file
        print "ERROR: $naccess_file is empty or does not exist. Exiting!\n";
        exit;
    }
    else{
        # load naccess file
        my @rsalines = $naccess_file->lines;
        
        foreach my $l (@rsalines){
            
            chomp($l);
            if($l=~ /^RES/){
                
                my $residue_num = substr($l, 9, 5);
                $residue_num=~ s/^\s+|\s+$//g;
                chomp($residue_num);
                my $rsa_allatoms = substr($l, 22, 6)+0;
                my $rsa_totside = substr($l, 35, 6)+0;
                my $rsa_mainchain = substr($l, 48, 6)+0;
                my $rsa_nonpolar = substr($l, 61, 6)+0;
                my $rsa_polar = substr($l, 74, 6)+0;
                
                if($rsa_allatoms == -99.9) {
                    print "$rsa_allatoms in $DOMAIN is -99.9! Exiting!\n";
                    exit;
                }

                $dom_feature_data{rsa_allatoms}{$residue_num} = $rsa_allatoms;    #*#
                $dom_feature_data{rsa_totside}{$residue_num} = $rsa_totside;    #*#
                $dom_feature_data{rsa_mainchain}{$residue_num} = $rsa_mainchain;    #*#
                $dom_feature_data{rsa_nonpolar}{$residue_num} = $rsa_nonpolar;    #*#
                $dom_feature_data{rsa_polar}{$residue_num} = $rsa_polar;    #*#
                
                # Find distance of residue to CM
                my $res_dist_cm_file = $dir_out->path("residue_distance_to_CM");
                system("perl $APPS_DIR/find_pdb_distance_pairs.pl $DOMPDB $residue_num $CM_file > $res_dist_cm_file");
                my $res_dist_cm = $res_dist_cm_file->lines;
                chomp($res_dist_cm);
                $res_dist_cm = nearest(0.01, $res_dist_cm);
                $dom_feature_data{dist_to_CM}{$residue_num} = $res_dist_cm;    #*#
                
            }
            
        }
    }
    
}
    

sub load_cleft_features{

    my ($speedfill_file) = shift;
    
    if (!-e "$speedfill_file" || -z "$speedfill_file") {
        # generate speedfill file
        print "ERROR: $speedfill_file is empty or does not exist. Exiting!\n";
        exit;
    }
    else{
        # load speedfill file
        my @cleftlines = $speedfill_file->lines;
        
        foreach my $line (@cleftlines){
            
            chomp($line);
            if($line=~ /^ATOM/){
                my $residue_num = substr($line, 22, 5);
                $residue_num=~ s/^\s+|\s+$//g;
                chomp($residue_num);
                my $cleft_num = substr($line, 56, 4)+0;
                my $cleft_depth = substr($line, 60, 6)+0;
                
                $dom_feature_data{cleft_num}{$residue_num} = 0;
                $dom_feature_data{cleft_depth}{$residue_num} = 0;
                
                if($dom_feature_data{cleft_num}{$residue_num}){
                    
                    my $temp = $dom_feature_data{cleft_num}{$residue_num};
                    if ($cleft_num !=0 && $cleft_num < $temp) {
                        $dom_feature_data{cleft_num}{$residue_num} = $cleft_num;        #*#
                        $dom_feature_data{cleft_depth}{$residue_num} = $cleft_depth;    #*#
                    }
                    
                }
                else{
                    
                    $dom_feature_data{cleft_num}{$residue_num} = $cleft_num;        #*#
                    $dom_feature_data{cleft_depth}{$residue_num} = $cleft_depth;    #*#
                    
                }
                
            }
            
        }
    }
    
}
    



sub load_dssp_features{

    my $dssp_file = shift;
    my $dssp_obj = Bio::Structure::SecStr::DSSP::Res->new('-file'=>"$dssp_file");

    my @residues_ids = $dssp_obj->residues();
  	my %dssp_res=();
    
    foreach my $id (@residues_ids){
        
        chomp($id);
        $id=~ /^(\S+)\:/;
        my $residue_num = $1;
		
        my $sec_str = $dssp_obj->resSecStrSum($residue_num);
		$sec_str=~ s/^\s+|\s+$//g;
									
        #my $dssp_acc = $dssp_obj->resSolvAcc($residue_num);
        #my $dssp_surfarea = $dssp_obj->resSurfArea($residue_num);
        
		my $phi = $dssp_obj->resPhi($residue_num);						
        my $psi = $dssp_obj->resPsi($residue_num);						
        my $oBonds_ptr = $dssp_obj->resHB_O_HN($residue_num);
        my $nhBonds_ptr = $dssp_obj->resHB_NH_O($residue_num);
        my $resTco = $dssp_obj->resTco($residue_num);
        my $kappa = $dssp_obj->resKappa($residue_num);
        my $alpha = $dssp_obj->resAlpha($residue_num);
        
        $dom_feature_data{dssp_type}{$residue_num} = $sec_str;
        $dom_feature_data{phi}{$residue_num} = $phi;
        $dom_feature_data{psi}{$residue_num} = $psi;
        $dom_feature_data{oBonds_ptr}{$residue_num} = $oBonds_ptr;
        $dom_feature_data{nhBonds_ptr}{$residue_num} = $nhBonds_ptr;
        $dom_feature_data{resTco}{$residue_num} = $resTco;
        $dom_feature_data{kappa}{$residue_num} = $kappa;
        $dom_feature_data{alpha}{$residue_num} = $alpha;
        
	}

}


sub load_foldx_features{
    
    my $foldx_file = shift;
    
    if(-e "$foldx_file" && !-z "$foldx_file"){
			
        my @info = $foldx_file->lines;
        
        foreach my $l (@info){
        
            $l=~ /^\w+\s(\-*\w+)\sto ALA energy change is\s(.*)$/;
            my $residue_num = $1;
            my $foldx_value = "$2";
            $foldx_value = nearest(0.001,$foldx_value);
            $dom_feature_data{foldx_alascan}{$residue_num} = $foldx_value;
            
        }
    }
    else{
        print "ERROR: $foldx_file is empty or does not exist. Exiting!\n";
        exit;
    }
    
}


sub load_psaia_features{
    
    my $psaia_file = shift;
    
    if(-e "$psaia_file" && !-z "$psaia_file"){
		
        my @psais_lines= $psaia_file->lines;
		my $count=0;
		
        foreach my $line(@psais_lines){
				
            if($count >= 8){
                    
		chomp($line);
		my @tab = split /\s+/,$line;
		my $residue_num = $tab[7];
                
		my $avg_dpx = $tab[19]; $avg_dpx = nearest(0.01, $avg_dpx);
                my $max_dpx = $tab[23]; $max_dpx = nearest(0.01, $max_dpx);
                my $min_dpx = $tab[24]; $min_dpx = nearest(0.01, $min_dpx);
        
				my $avg_cx = $tab[25]; $avg_cx = nearest(0.01, $avg_cx);
                my $max_cx = $tab[29]; $max_cx = nearest(0.01, $max_cx);
                my $min_cx = $tab[30]; $min_cx = nearest(0.01, $min_cx);
				my $hydrophobicity_psaia = $tab[31];
		
                $dom_feature_data{avg_dpx}{$residue_num} = $avg_dpx;
                $dom_feature_data{max_dpx}{$residue_num} = $max_dpx;
                $dom_feature_data{min_dpx}{$residue_num} = $min_dpx;
                
                $dom_feature_data{avg_cx}{$residue_num} = $avg_cx;
                $dom_feature_data{max_cx}{$residue_num} = $max_cx;
                $dom_feature_data{min_cx}{$residue_num} = $min_cx;
                
                $dom_feature_data{hydrophobicity_psaia}{$residue_num} = $hydrophobicity_psaia;
                
		}
                
		$count++;
                
		}   
    }
    else{
        print "ERROR: $psaia_file is empty or does not exist. Exiting!\n";
        exit;
    }
}


sub load_centrality_feature{
    
    my ($centrality_feature, $centrality_file) = @_;
    
    if (-e "$centrality_file" && !-z "$centrality_file") {
        
        my %centrality=();
            my @lines = $centrality_file->lines;
            foreach my $line (@lines){
                chomp($line);
                my ($num, $centrality_measure) = split("\t", $line);
            $centrality{$num} = $centrality_measure;
        }
        my $count=0;
             for my $residue_num ( sort keys %{ $dom_feature_data{scons} } ) {
            $count++;
            $dom_feature_data{$centrality_feature}{$residue_num} = $centrality{$count};
            
        }
        
    }
    else{
        print "ERROR: $centrality_file is empty or does not exist. Exiting!\n";
        exit;
    }
    
}


sub load_distmap{
    
    my ($DISTMAP) = @_;
    
    if (-e $DISTMAP && !-z $DISTMAP) {
        my @dist_lines = $DISTMAP->lines;
        foreach my $line (@dist_lines){
            chomp($line);
            my ($r1, $r2, $d) =split("\t", $line);
            $distmap{$r1}{$r2} = $d;
        }
    }
    else{
        print "ERROR: $DISTMAP is empty or does not exist. Exiting!\n";
        exit;
    }
}


sub find_distance{
    
    my ($residue1, $residue2) = @_;
    
    my $min_distance;
    
    if ($residue1 eq $residue2) {
        $min_distance = 0;
    }
    elsif($distmap{$residue1}{$residue2}){
        $min_distance = $distmap{$residue1}{$residue2};
    }
    elsif($distmap{$residue2}{$residue1}){
        $min_distance = $distmap{$residue2}{$residue1};
    }
    else{
        print "ERROR: Cannot find distance $residue1-$residue2 in $DISTMAP! Exiting!\n";
        exit;
    }
    
    return $min_distance;
}


sub load_min_distance_to_cleft{

    # for each residue in the pdb look at the cleft residues and find min distance
    for my $residue_num ( sort keys %{ $dom_feature_data{scons} } ) {
        
        my %min_cleft_dist = ();

        for my $residue_iter ( sort keys %{ $dom_feature_data{scons} } ) {
            
            my $cleft_iter = 1;
            while ($cleft_iter <= 3) {
            
                if ($dom_feature_data{cleft_num}{$residue_iter} == $cleft_iter) {
                    
                    my $min_distance = &find_distance($residue_num,$residue_iter);
                    
                    if($min_cleft_dist{$cleft_iter}){
                        my $temp = $min_cleft_dist{$cleft_iter};
                        if ($min_distance < $temp) {
                            $min_cleft_dist{$cleft_iter} = $min_distance;
                        }
                    }
                    else{
                        
                        $min_cleft_dist{$cleft_iter} = $min_distance;
                    }
                    
                }
            
                $cleft_iter++;
            }
                
         }
        
        # add the min distance of the residue to the clefts to feature hash
        foreach my $cleft (sort keys %min_cleft_dist){
            
            chomp($cleft);
            my $cleft_string = "min_dist_to_cleft_"."$cleft";
            $dom_feature_data{$cleft_string}{$residue_num} = $min_cleft_dist{$cleft};
            #print "$residue_num\t$cleft_string\t$min_cleft_dist{$cleft}\n";
            
        }  
    }

}


sub load_min_distance_to_surface{
    
    #suraface is regarded as any residue with >= 25% RSA (rsa_allatoms)
    
    # for each residue in the pdb look at the surface residues and find min distance
    for my $residue_num ( sort keys %{ $dom_feature_data{scons} } ) {
        
        my $min_rsa_dist;

        for my $residue_iter ( sort keys %{ $dom_feature_data{scons} } ) {
            
            if ($dom_feature_data{rsa_allatoms}{$residue_iter} >= 25) {
                    
                my $min_distance = &find_distance($residue_num,$residue_iter);
                    
                if ($min_rsa_dist) {
                    
                    if($min_distance < $min_rsa_dist){
                        $min_rsa_dist = $min_distance;
                    }
                }
                else{
                    $min_rsa_dist = $min_distance;
                }
                
            }
                
         }
        
        # add the min distance of the residue to the surface to feature hash
        $dom_feature_data{dist_to_surface}{$residue_num} = $min_rsa_dist;
    }

}



sub load_structural_neighborhood_properties{
    
    my $distance_cutoff = shift;
    
    # for each residue, find residues within distance cutoff and get avg features for amino acid properties
    
    for my $residue_num ( sort keys %{ $dom_feature_data{scons} } ) {
        
        my %neighbours=();
        
        # add itself to the structural neighbourhood
        $neighbours{$residue_num}=1;
        
        my $scons_sum = 0;
        my $surface_residues=0;
        my $polarity_sum=0;
        my $hydrophobicity_sum=0;
        my $hydropathicity_sum=0;
        my $charged_sum=0;
        my $electric_effect_sum=0;
        
        for my $residue_iter ( sort keys %{ $dom_feature_data{scons} } ) {
         
            my $distance = &find_distance($residue_num,$residue_iter);
            
            if ($distance <= $distance_cutoff) {
                
                $neighbours{$residue_iter}=1;
                $scons_sum = $scons_sum + $dom_feature_data{scons}{$residue_iter};
                if ($dom_feature_data{rsa_allatoms}{$residue_iter} >= 25) {
                    $surface_residues++;
                }
                my $aminoacid = $dom_feature_data{residue_aa}{$residue_iter};
                if ($aminoacid ne "X") {
                
                    $polarity_sum = $polarity_sum + $aa_indices{$aminoacid}{polarity};
                    $hydrophobicity_sum = $hydrophobicity_sum + $aa_indices{$aminoacid}{hydrophobicity};
                    $hydropathicity_sum = $hydropathicity_sum + $aa_indices{$aminoacid}{hydropathicity};
                    $charged_sum = $charged_sum + $aa_indices{$aminoacid}{charge};
                    $electric_effect_sum = $electric_effect_sum + $aa_indices{$aminoacid}{localised_electrical_effect};
                }
            }
              
        }
        
        my $neighbour_count = keys %neighbours;
        
        my $avg_scons = $scons_sum/$neighbour_count;
        $avg_scons = nearest(0.01, $avg_scons);
        $dom_feature_data{avg_scons}{$residue_num}= $avg_scons;
        
        my $avg_surface_residues = $surface_residues/$neighbour_count;
        $avg_surface_residues = nearest(0.01, $avg_surface_residues);
        $dom_feature_data{avg_surface_residues}{$residue_num}= $avg_surface_residues;
        
        my $avg_polarity = $polarity_sum/$neighbour_count;
        $avg_polarity = nearest(0.01,$avg_polarity);
        $dom_feature_data{avg_polarity}{$residue_num}= $avg_polarity;
        
        my $avg_hydrophobicity = $hydrophobicity_sum/$neighbour_count;
        $avg_hydrophobicity = nearest(0.01,$avg_hydrophobicity);
        $dom_feature_data{avg_hydrophobicity}{$residue_num}= $avg_hydrophobicity;
        
        my $avg_hydropathicity = $hydropathicity_sum/$neighbour_count;
        $avg_hydropathicity = nearest(0.01,$avg_hydropathicity);
        $dom_feature_data{avg_hydropathicity}{$residue_num}= $avg_hydropathicity;
        
        my $avg_charged = $charged_sum/$neighbour_count;
        $avg_charged = nearest(0.01,$avg_charged);
        $dom_feature_data{avg_charged}{$residue_num}= $avg_charged;
        
        my $avg_electric_effect = $electric_effect_sum/$neighbour_count;
        $avg_electric_effect = nearest(0.01,$avg_electric_effect);
        $dom_feature_data{avg_electric_effect}{$residue_num}= $avg_electric_effect;
        
    }
    
}


sub load_norm_bfactors{
    
    my $dompdb = shift;
    
    my %bfactors=();
    
    if (-e $dompdb){
    
        my @pdblines = $dompdb->lines;
        
        foreach my $pdbline (@pdblines){
            chomp($pdbline);
            my $r = substr($pdbline, 22, 5);
            chomp($r);
            $r=~ s/^\s+|\s+$//g;
            my $value = substr($pdbline, 60, 6);
            chomp($value);
            
            if($bfactors{$r}){
                my $temp = $bfactors{$r};
                if ($value < $temp) {
                    $bfactors{$r}=$value;
                }
            }
            else{
                $bfactors{$r}=$value;
            }
            
        }
        
        my $maxbfactor = max values %bfactors;
        my $minbfactor = min values %bfactors;
        
        my $chkbfactorrange = $maxbfactor - $minbfactor;
        
        for my $residue_num ( sort keys %{ $dom_feature_data{scons} } ) {
            if($chkbfactorrange != 0){
                
                foreach my $residue_num (keys %bfactors){
                    
                    my $res_bfactor = $bfactors{$residue_num};
                    my $res_bfactor_n = ($res_bfactor - $minbfactor)/($maxbfactor - $minbfactor);
                    $res_bfactor_n = nearest(0.01, $res_bfactor_n);
                    $dom_feature_data{res_bfactor_n}{$residue_num} = $res_bfactor_n;
                    
                }
                
            }
            else{
                $dom_feature_data{res_bfactor_n}{$residue_num} = 0;
            }
        }
    
    }
    undef(%bfactors);
}


sub load_ff_ibis_freq{
    
    my $ibis_freq_file = shift;
    
    if (-e "$ibis_freq_file") {
        
        my @lines = $ibis_freq_file->lines;
        
        foreach my $line (@lines){

            my ($ff_name, $ibis_type, $aln_pos, $freq_score, $num_pdbs, $tot_pdbs) = split /\s{2,}/, $line;
            if ($freq_score > 1) {
                $freq_score = 1;
            }
            my $ibis_type_freq = "$ibis_type"."_freq";
            my $residue_num = $aln_pos_map{$aln_pos};
            
            $dom_feature_data{$ibis_type_freq}{$residue_num}= $freq_score;  #*#
            
            
      
        }
        
    }
    
    #assign zero
    for my $residue_num ( sort keys %{ $dom_feature_data{scons} } ) {
        unless($dom_feature_data{ibis_lig_freq}{$residue_num}){
               $dom_feature_data{ibis_lig_freq}{$residue_num}=0;
        }
        unless($dom_feature_data{ibis_ppi_freq}{$residue_num}){
            $dom_feature_data{ibis_ppi_freq}{$residue_num}=0;
        }
    } 
    
}


sub load_aaindices{
	
    my $aa_indices_file = shift;
    my @AAindices_lines = $aa_indices_file->lines;

    my @headers;
    
    foreach my $line (@AAindices_lines){
        
        chomp($line);

        my @indices = split("\t", $line);
        
        if($line =~ /^#/){
            
            @headers= @indices;
            shift @headers;
            #print "@headers\n";
            
        }
        else{
            
            my $i=0;
            my $AA = shift @indices;
        
            foreach my $col (@indices){
            
                $aa_indices{$AA}{$headers[$i]}=$col;
                $i++;
            }
        }

    }
    
    return @headers;

}
sub load_annotations{
    
    my $DOMAIN = shift;
    $DOMAIN=~ /^(\w\w\w\w)(\w)/;
    my $PDB = $1;
    my $CHAIN = $2;
    my $STRUC = "$PDB"."$CHAIN";
    
    my $domcsa_info_file = path("/cath/people2/ucbtdas/git/cath-funsite/funsites-ml/generate_domain_features/data/M-CSA/cath_domain_MCSA_dataset_good_DOPS.tab");
    my $domlig_info_file = path("/cath/people2/ucbtdas/FUNSITE-ML/dataset/IBIS/cath_domain_LIG_dataset_good_DOPS_3Apdbs_NR.tab");
    my $domppi_info_file = path("/cath/people2/ucbtdas/FUNSITE-ML/dataset/IBIS/cath_domain_PPI_dataset_good_DOPS_3Apdbs_NR.tab");
    my $dom_othersite_info_file = path("/cath/people2/ucbtdas/FUNSITE-ML/dataset/IBIS/cath_domain_OTHER_dataset_good_DOPS_3Apdbs.tab");
    my $dom_biolip_info_file = path("/cath/people2/ucbtdas/git/cath-funsite/funsites-ml/generate_domain_features/data/BioLip/cath_domain_BIOLIP_LIG_small_molecules_dataset_good_DOPS.tab");
    my $did3_file_interchain = path("/export/sayoni/3did/3did_flat_interchain");
    my $did3_file_intrachain = path("/export/sayoni/3did/3did_flat_intrachain");
    my $proteindb_annofile= path("/export/sayoni/protindb/DS4000/parsed/$STRUC.protindb.anno");
    
    my @domain_csa_lines = `LC_ALL=C fgrep -w "$DOMAIN" $domcsa_info_file`;
    foreach my $domline (@domain_csa_lines){
        chomp($domline);
        my @tab = split("\t", $domline);
        if($tab[0] eq $DOMAIN){
            my $csa_res = $tab[6];
            my $csa_roletype = $tab[8];
            my $csa_role = $tab[9];
            $dom_feature_data{annotation_MCSA}{$csa_res} = 1;    #*#
            $dom_feature_data{annotation_MCSA_roletype}{$csa_roletype} = 1;    #*#
            $dom_feature_data{annotation_MCSA_role}{$csa_role} = 1;    #*#
        }
    }
    
    
    my @domain_lig_lines = `LC_ALL=C fgrep -w "$DOMAIN" $domlig_info_file`;
    foreach my $line (@domain_lig_lines){
        chomp($line);
        my @tab = split("\t", $line);
        if($tab[0] eq $DOMAIN){
            my $ibis_lig_res = $tab[6];
            $dom_feature_data{annotation_IBIS_LIG}{$ibis_lig_res} = 1;    #*#
        }
    }
    
    my @domain_ppi_lines = `LC_ALL=C fgrep -w "$DOMAIN" $domppi_info_file`;
    foreach my $line (@domain_ppi_lines){
        chomp($line);
        my @tab = split("\t", $line);
        if($tab[0] eq $DOMAIN){
            my $ibis_ppi_res = $tab[6];
            my $intrachain = $tab[9];
            
            if ($intrachain == 1) {
                $dom_feature_data{annotation_IBIS_PPI_INTRACHAIN}{$ibis_ppi_res} = 1;    #*#
            }
            else{
                $dom_feature_data{annotation_IBIS_PPI_INTERCHAIN}{$ibis_ppi_res} = 1;    #*#
            }
        }
    }
    
    my @biolip_lines = `LC_ALL=C fgrep -w "$DOMAIN" $dom_biolip_info_file`;
    foreach my $line (@biolip_lines){
        chomp($line);
        my @tab = split("\t", $line);
        if($tab[0] eq $DOMAIN){
            my $biolip_lig = $tab[6];
            my $biolip_ligand_name = $tab[7];
            $dom_feature_data{annotation_BIOLIP}{$biolip_lig} = 1;    #*#
            $dom_feature_data{annotation_BIOLIP_ligand}{$biolip_ligand_name} = 1;    #*#
        }
    }
    
    my @did3_info_interchain = `LC_ALL=C fgrep -w "$PDB" $did3_file_interchain`;
	foreach my $line (@did3_info_interchain){
	 	chomp($line);
	 	my @array = split("\t",$line);
	 	my $did3_res_interchain = $array[2];
	 	$did3_res_interchain=~ s/\s+//g;
	 	if($array[0] eq $PDB && $array[1] eq $CHAIN){
	 		$dom_feature_data{annotation_3DID_INTERCHAIN}{$did3_res_interchain} = 1;    #*#
	 	}
	}
	 		
    my @did3_info_intrachain = `LC_ALL=C fgrep -w "$PDB" $did3_file_intrachain`;
	foreach my $line (@did3_info_intrachain){
	 	chomp($line);
	 	my @array = split("\t",$line);
	 	my $did3_res_intrachain = $array[2];
	 	$did3_res_intrachain=~ s/\s+//g;
	 	if($array[0] eq $PDB && $array[1] eq $CHAIN){
	 		$dom_feature_data{annotation_3DID_INTRACHAIN}{$did3_res_intrachain} = 1;    #*#
	 	}
	}
    
	if(-e "$proteindb_annofile"){
	 	my @proteindb_info = `LC_ALL=C fgrep -w "$STRUC" $proteindb_annofile`;
	 	if(@proteindb_info){
	 		foreach my $anno (@proteindb_info){
	 			chomp($anno);
	 			my @t = split("\t", $anno);
	 			if($t[1] eq $PDB && $t[2] eq $CHAIN){
                    my $res = $t[3];
                    $dom_feature_data{annotation_PROTINDB}{$res} = 1;    #*#
	 			}
	 		}
	 	}
	}
    
    #assign zero
    for my $residue_num ( sort keys %{ $dom_feature_data{scons} } ) {
        unless($dom_feature_data{annotation_MCSA}{$residue_num}){
               $dom_feature_data{annotation_MCSA}{$residue_num}=0;
        }
        unless($dom_feature_data{annotation_MCSA_role}{$residue_num}){
               $dom_feature_data{annotation_MCSA_role}{$residue_num}="NO_ROLE";
        }
        unless($dom_feature_data{annotation_MCSA_roletype}{$residue_num}){
               $dom_feature_data{annotation_MCSA_roletype}{$residue_num}="NO_ROLETYPE";
        }
        unless($dom_feature_data{annotation_IBIS_LIG}{$residue_num}){
            $dom_feature_data{annotation_IBIS_LIG}{$residue_num}=0;
        }
        unless($dom_feature_data{annotation_PROTINDB}{$residue_num}){
            $dom_feature_data{annotation_PROTINDB}{$residue_num}=0;
        }
        unless($dom_feature_data{annotation_3DID_INTRACHAIN}{$residue_num}){
            $dom_feature_data{annotation_3DID_INTRACHAIN}{$residue_num}=0;
        }
        unless($dom_feature_data{annotation_3DID_INTERCHAIN}{$residue_num}){
            $dom_feature_data{annotation_3DID_INTERCHAIN}{$residue_num}=0;
        }
        unless($dom_feature_data{annotation_BIOLIP}{$residue_num}){
            $dom_feature_data{annotation_BIOLIP}{$residue_num}=0;
        }
        unless($dom_feature_data{annotation_BIOLIP_ligand}{$residue_num}){
            $dom_feature_data{annotation_BIOLIP_ligand}{$residue_num}="NO_LIGAND";
        }
        unless($dom_feature_data{annotation_IBIS_PPI_INTERCHAIN}{$residue_num}){
            $dom_feature_data{annotation_IBIS_PPI_INTERCHAIN}{$residue_num}=0;
        }
        unless($dom_feature_data{annotation_IBIS_PPI_INTRACHAIN}{$residue_num}){
            $dom_feature_data{annotation_IBIS_PPI_INTRACHAIN}{$residue_num}=0;
        }
    }
            
}
