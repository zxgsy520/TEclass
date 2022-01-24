#!/usr/bin/env perl

# A script to create decision trees to classify unknown repetitive elements. 

use warnings;
use Getopt::Std;
use Dumpvalue;
use lib X 
use Configuration;
use TEclass;

getopts('o:O:W:'); 
our ($opt_o, $opt_O, $opt_W);


# Sets the working directory
my $workdir;
if ($opt_W) {
    $workdir = makedir($opt_W); 
}
else {
    $workdir = makedir("workdir");
} 

# Sets the output dir
my $path_to_RF;
if ($opt_O) {
    $path_to_RF = makedir("$opt_O/RandomForests", 1);   
}
else {
    $path_to_RF = makedir("$path_to_TEclass/classifiers/RandomForests", 1);
} 


# Builds a hash of oligos depending on the specified oligomer size, 
# the default is tetramers
if (defined $opt_o && $opt_o == 5) {
    build_pentamer_hash();
    our $features_no = 1024;
    $path_to_RF  .= "5";
}
else {
    build_tetramer_hash();
    our $features_no = 256;   
}


my @files = qw/ short medium long xlong /;
my @analysis = qw/ forward_vs_reverse 
                   DNA_vs_Retro
                   LTR_vs_nonLTR
                   LINE_vs_SINE /;

#-------------------------------------------------------------------------------
#                               MAIN - Builds decision trees
#-------------------------------------------------------------------------------

# Loops through the files and the different comparisons
foreach my $j (0..$#analysis) {  
    print STDERR ($j+1), ". $analysis[$j]\n";  
    foreach my $i (0..$#files) { 
        # Does not build test set, because there are no long SINEs
        if ($files[$i] =~ m/long|xlong/ && $analysis[$j] eq 'LINE_vs_SINE') {
           last;
        }    
        
        build_repeat_vectors("$workdir/$files[$i]", $analysis[$j], 1, 0);
        standardise_vectors_train("$workdir/$files[$i]", $analysis[$j], $path_to_svm, 0.0001);
        translate_from_SVM_to_CSV("$workdir/$files[$i]", $analysis[$j]);
        train_RF("$workdir/$files[$i]", $analysis[$j], $features_no, $path_to_RF);        
    }
}


#-------------------------------------------------------------------------------    
#                               - Subroutines -   
#    train_RF     
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
# Purpose   : Create a random forest model      
# Usage     : &translate_from_SVM_to_LVQ(file, analysis_type,  dimensions,)
# Arguments : string, string, integer, 
# Returns   : nothing, prints
# Globals   : none
#******************
sub train_RF {
    my ($file, $analysis_type, $features, $path_to_RF) = @_;
    print STDERR "  - Training decision trees ... ";
    my $cmd =   "$path_to_librf/rftrain " 
           . "-t 500 "
           . "-m $file.$analysis_type.rf.model "
           . "-d $file.$analysis_type.csv " 
           . "-f $features --csv "
           . "-l $file.$analysis_type.lab " 
           . " &>$file.$analysis_type.rf.log";    
    system("$cmd"); 
    
    $cmd = "mv $file.$analysis_type.csv " 
           . " $file.$analysis_type.lab "
           . " $file.$analysis_type.rf.model "
           . " $file.$analysis_type.rf.log "
           . " $file.$analysis_type.range "
           . " $path_to_RF";
    system("$cmd");    
    
    #system("mv $file.$analysis_type.range $path_to_RF");
    print STDERR "done\n";
}

