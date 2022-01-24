#!/usr/bin/env perl

# A script to create Learning vectors to classify unknown repetitive elements. 
 
use warnings;
use Getopt::Std;
use Dumpvalue;
use lib X 
use Configuration;
use TEclass;

getopts('hO:W:'); 
our($opt_h, $opt_O, $opt_W);


# Pints help
if ($opt_h) {
    print STDERR "\n"
       . "Usage:\n" 
       . " ./Build_LVQ_sets.pl [-h] -W 'workdir' -O 'outputdir' \n\n"
       . "flags:\n" 
       . "-W 'workdir' : Sets the working directory for temp. files.\n"
       . "-O 'outdir' : Sets the output directory.\n"      
       . "-h : this help\n";
    exit;
}

# Sets the working directory
my $workdir;
if ($opt_W) {
    $workdir = makedir($opt_W); # TEclass_build checks the path for correctness  
}
else {
    $workdir = makedir("workdir");
} 


# Sets the output dir
my $path_to_codebooks;
if ($opt_O) {
    $path_to_codebooks = makedir("$opt_O/LVQ_codebooks", 1); # TEclass_build checks the path for correctness  
}
else {
    $path_to_codebooks = makedir("$path_to_TEclass/classifiers/LVQ_codebooks", 1);
} 
 


#-------------------------------------------------------------------------------
#                                 - Making the LVQ -
#-------------------------------------------------------------------------------

#Builds a hash of oligos
build_tetramer_hash();
my @files = qw/ short medium long xlong /;                       
my $features = 256;

#-------------------------- 1. Forward vs. Reverse strand ----------------------
print STDERR "1. Reverse vs. Forward strand LVQ.\n";
foreach my $i(0..$#files) { 
    build_LVQ_repeat_vectors("$workdir/$files[$i]", 'forward_vs_reverse');

    standardise_vectors_train("$workdir/$files[$i]", 'forward_vs_reverse', $path_to_svm, 0.0001);
    system("mv $workdir/*.range $path_to_codebooks");
    
    translate_from_SVM_to_LVQ("$workdir/$files[$i]", 'forward_vs_reverse', $features);
    system("mv $workdir/*.codebook $path_to_codebooks");
}


#------------------------ 2. DNA vs. LTR vs. LINEs vs. SINEs  ------------------
print STDERR "2. Repeat LVQ\n";
foreach my $i (0..$#files) { 
    build_LVQ_repeat_vectors("$workdir/$files[$i]", 'repeats');
    
    standardise_vectors_train("$workdir/$files[$i]", 'repeats', $path_to_svm, 0.0001);
    system("mv $workdir/*.range $path_to_codebooks");
    
    translate_from_SVM_to_LVQ("$workdir/$files[$i]", 'repeats', $features);
    system("mv $workdir/*.codebook $path_to_codebooks");
}


#-------------------------------------------------------------------------------    
#                               - Subroutines -   
#    build_repeat_vectors     
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Purpose   : To build a vector of oligomer frequencies for every repeat           
# Usage     : build_repeat_vectors($string, $string);
# Arguments : $file, $analysis type
# Returns   : nothing, prints  
# Globals   : none 
#******************
sub build_LVQ_repeat_vectors {
    my ($file, $analysis_type) = @_;
    
    open(IN, "<", "$file.lib") or die "cannot open $file.lib $!";
    open(OUT, ">", "$file.$analysis_type.vect") 
                             or die "cannot open $file.$analysis_type.vect $!";
    my $outfile = *OUT;
    print STDERR "$file repeats:\n"; 
    print STDERR "  - Writing training set ... ";
    # Processes the linearised .lib file
    while (1) {
        my $repeat = <IN>;
        $repeat =~ s/>//;
        chomp($repeat); 
        my $sequence = <IN>;
        if ($analysis_type eq 'forward_vs_reverse') {
            # Sequence in forward direction
            process_sequence($sequence, 1, $outfile);
            # Reverse complemented sequence
            $sequence = reverse_complement($sequence);
            process_sequence($sequence, 0, $outfile);
        }
        elsif ($analysis_type eq 'repeats') {
            my @tags = split/\s+/, $repeat;
            my @prefix = split/_/, $tags[0];  
            # Assigns a numeric value for the different repeat types 
            # (to make svm-scale happy)
            my $label = 0;
            if    ($prefix[-1] eq DNA)  { $label = 1; } 
            elsif ($prefix[-1] eq LTR)  { $label = 2; }
            elsif ($prefix[-1] eq LINE) { $label = 3; }
            elsif ($prefix[-1] eq SINE) { $label = 4; }

            process_sequence($sequence, $label, $outfile);
        }
        last if eof(IN);
    }
    print STDERR "done.\n";
    close OUT; 
    close IN;
}
