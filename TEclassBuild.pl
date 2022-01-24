#!/usr/bin/env perl
# A script to automate the creation of SVM models that are used to classify 
# putative TE consensi.  

use warnings;
use Getopt::Std;
use lib X    
use Configuration;
use TEclass;

getopts('hs:m:p:o:w:x:C:g:');
our($opt_h, $opt_s, $opt_m, $opt_p, $opt_o, $opt_w, $opt_x, $opt_C, $opt_g);

# Defaults
if (!defined $opt_m) {  $opt_m = 1024; }
if (!defined $opt_p) {  $opt_p = 1; }
if (!defined $opt_C) {  $opt_C = 17; }
if (!defined $opt_g) {  $opt_g = -11; }


# Usage
if ($opt_h) {
    print STDERR "\n"
       . "Usage:\n" 
       . " ./TEclassBuild.pl -x <int> [+ optional flags]\n\n"
       . "Flags:\n"  
       . "-x <int> : the parameter selection mode of libSVM \n"
       . "     0   : parameter search = grid   \n"
       . "     1   : parameter search = list  \n"
       . "     2   : user defined C and g (defaults 17 and -11 respectively)\n"                       
       . "-C <int> : user defined C (only if -x is set to 2) \n"
       . "-g <int> : user defined g (only if -x is set to 2) \n\n"   
       . "-w 'workdir': Sets the working directory for tmp files. Default is ./workdir\n"
       . '-o \'outdir\' : Sets the output directory. Default is $TEclassdir/classifiers', "\n"       
       . "-s <int> : builds SVM models on samples of specified size,\n"
       . "-m <int> : sets the memory allocated to each CPU for SVM (default is 1024 MB)\n"
       . "-p <int> : sets the number of CPUs used in SVM training (default is 1)\n\n"       
       . "-h : this help\n";
    exit;
}

# Checks for the necessary params 
if (!defined $opt_x) {
    print STDERR "You must specify the parameter search options (-x).\n";
    exit;
}
elsif ($opt_x != 0 && $opt_x != 1 && $opt_x != 2) {
    print STDERR "The value of -x must be 0, 1 or 2. \n",
                 "Please specify a correct value.\n";
    exit;             
}

# Warns if a user provides too many processors
if ($opt_x == 0 && $opt_p > 99) {
    print STDERR "\nToo many processors specified, TEclass will use only 99 for SVM training.\n"; 
}
if ($opt_x == 1 && $opt_p > 12) {
    print STDERR "\nToo many processors specified, TEclass will use only 12 for SVM training.\n"; 
}
if ($opt_x == 2 && $opt_p > 1) {
    print STDERR "\nThe -x 2 option is not parallelized, only 1 CPU will be used for SVM training.\n"; 
}


# Sets the working directory
my $workdir;
if ($opt_w) { 
    $workdir = makedir($opt_w);
}
else {        
    $workdir = makedir("workdir"); 
} 


# Sets the output dir
my $outdir;
if ($opt_o) {
    $outdir = makedir($opt_o, 1);  
}
else {
    $outdir = makedir("$path_to_TEclass/classifiers", 1);
} 

#--------------------------- Building the classifiers --------------------------

# Preprocessing - creates .lib files from the GIRI and Repeatmasker editions of 
# RepBase 
print "\n### Preprocessing\n";
system("perl $path_to_TEclass/Preprocessing.pl -W $workdir");

# Builds ICM models with Glimmer 
print "\n### Building gene models\n";
system("perl $path_to_TEclass/Build_ICM_model.pl -W $workdir -O $outdir ");

# Builds codebooks 
print "\n### Building LVQ codebooks\n";
system("perl $path_to_TEclass/Build_LVQ_sets.pl -W $workdir -O $outdir ");

# Builds decision trees
print "\n### Building Random Forests\n";
system("perl $path_to_TEclass/Build_Random_Forests.pl -W $workdir -O $outdir ");

# Builds the SVM models
if ($opt_s) {  
    print "\n### Building tetramer SVM classifiers\n";
    system("perl $path_to_TEclass/Build_SVM_models.pl -s $opt_s -o 4 -m $opt_m -p $opt_p -W $workdir -O $outdir -x $opt_x -C $opt_C -g $opt_g"); 
    print "\n### Building pentamer SVM classifiers\n";
    system("perl $path_to_TEclass/Build_SVM_models.pl -s $opt_s -o 5 -m $opt_m -p $opt_p -W $workdir -O $outdir -x $opt_x -C $opt_C -g $opt_g");
}
else {
    print "\n### Building tetramer SVM classifiers\n";
    system("perl $path_to_TEclass/Build_SVM_models.pl -o 4 -m $opt_m -p $opt_p -W $workdir -O $outdir -x $opt_x -C $opt_C -g $opt_g"); 
    print "\n### Building pentamer SVM classifiers\n";
    system("perl $path_to_TEclass/Build_SVM_models.pl -o 5 -m $opt_m -p $opt_p -W $workdir -O $outdir -x $opt_x -C $opt_C -g $opt_g");
}

# Cleanup
`rm -r $workdir`;

