#!/usr/bin/env perl

# A script to create SVM databases for clasification of unknown repetitive 
# elements.

# TODO : Rewrite the threading part, the current version is temporary. 

use warnings;
use Getopt::Std;
use lib X 
use Configuration;
use TEclass;
use Dumpvalue;
use threads;

getopts('ho:s:m:p:O:W:x:C:g:'); 
our($opt_h,$opt_o,$opt_s,$opt_m,$opt_p,$opt_O,$opt_W,$opt_x,$opt_C,$opt_g);

#-------------------------------------------------------------------------------
#                                   Defaults
#-------------------------------------------------------------------------------

if ($opt_m) { our $memory = $opt_m; }
else {        our $memory = 1024; }

if ($opt_p) { our $processors = $opt_p; }
else {        our $processors = 1; }

if (! defined $opt_C) { $opt_C = 17; }
if (! defined $opt_g) { $opt_g =-11; }


# Pints help
if ($opt_h) {
    print STDERR "\n"
       . "Usage:\n" 
       . " ./Build_SVM_models.pl [-flags]\n\n"
       . "flags:\n"
       . "-W 'workdir': Sets the working directory for temp. files.\n"
       . "-O 'outdir' : Sets the output directory.\n"                
       . "-s <int>: Build models using samples of specified size [e.g. -s 300]\n"
       . "-o <int>: sets the oligo size [e.g. -o 4], 4 or 5-mers are allowed\n"
       . "-m <int>: sets the memory allocated to each CPU (default 1024 MB)\n"
       . "-p <int>: sets the number of CPUs used (default 1)\n\n"
       . "-x <int> : the parameter selection mode of the SVM \n"
       . "   = 0   : parameter search = grid   \n"
       . "   = 1   : parameter search = list  \n"
       . "   = 2   : user defined C and g (defaults are 17 and -11 respectively)\n"                       
       . "-C <int>: user defined C (only if -x is set 2) \n"
       . "-g <int>: user defined g (only if -x is set 2) \n\n"                     
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


if (!defined $opt_o) {
    print STDERR "You must specify the oligomer type (4 or 5).\n";
    exit;
}

# Sets the working directory
our $workdir;
if ($opt_W) {
    $workdir = makedir($opt_W);   
}
else {
    $workdir = makedir("workdir");
} 

 
# Sets the output depending on oligo 
our $path_to_oligomers;
if ($opt_o == 4) {
    build_tetramer_hash();
    if ($opt_O) {
        $path_to_oligomers = makedir("$opt_O/tetramer_models", 1);  
    }
    else {
        $path_to_oligomers = makedir("$path_to_TEclass/classifiers/tetramer_models", 1);
    } 
}
elsif ($opt_o == 5) {
    build_pentamer_hash();
    if ($opt_O) {
        $path_to_oligomers = makedir("$opt_O/pentamer_models", 1);  
    }
    else {
        $path_to_oligomers = makedir("$path_to_TEclass/classifiers/pentamer_models", 1);
    } 
}


#-------------------------------------------------------------------------------
#                               MAIN - Builds SVM models
#-------------------------------------------------------------------------------

my @files = qw/ short medium long xlong /;
my @analysis = qw/ forward_vs_reverse 
                   DNA_vs_Retro
                   LTR_vs_nonLTR
                   LINE_vs_SINE /;

open (OUT, ">", "$workdir/CrossValidations.txt") or die "cant open OUT $!";

# Loops through the files and the different comparisons
foreach my $j (0..$#analysis) {  
    print STDERR ($j+1), ". $analysis[$j]\n";  
    foreach my $i (0..$#files) { 
        # Does not build test set, because there are no long SINEs
        last if ($files[$i] =~ m/long|xlong/ && $analysis[$j] eq 'LINE_vs_SINE');    
        build_repeat_vectors("$workdir/$files[$i]", $analysis[$j], '-1', '1');
        sample_scale("$workdir/$files[$i]", $analysis[$j]); 
        
        # Grid search for best C and g
        if ($opt_x == 0) {
            my @params = ();
            for (my $k=-3; $k<=17; $k+=2) {        
                for (my $l=1; $l>=-15; $l-=2) {
                    push @params, [$k, $l];
                }
            }   
            my ($C, $g) = best_Cg("$workdir/$files[$i]", $analysis[$j], \@params);
            train_selection("$workdir/$files[$i]", $analysis[$j], $C, $g);
        }
        # Prededfined list
        elsif ($opt_x == 1) {
            my @params = ([ 1, -5], [ 1, -3], [ 3, -5], [ 5, -3], [11, -9], [13, -9], 
                          [15,-11], [17,-13], [17,-11], [17, -7], [17, -5], [17, -3]);        
            my ($C, $g) = best_Cg("$workdir/$files[$i]", $analysis[$j], \@params);            
            train_selection("$workdir/$files[$i]", $analysis[$j], $C, $g);
        } 
        # User defined C and g
        elsif ($opt_x == 2) {
            train_selection("$workdir/$files[$i]", $analysis[$j], $opt_C, $opt_g);
        }              
    }
}
close OUT;
system("mv $workdir/CrossValidations.txt $path_to_oligomers");


#--------------------------------- Subroutines ---------------------------------   
#    sample_scale
#    best_Cg  
#    train_selection
#    system_call
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Purpose   : takes samples fron,  and scales the inplut data    
# Usage     : sample_scale($file, $analysis_type); 
# Arguments : the file vith oligomer vectors, the type of the analysis 
# Returns   : nothing, prints
# Globals   : none
#********************** 
sub sample_scale { 
    my ($file, $analysis_type) = @_; 
    
    # The -s flag is on: take samples. This should be used mainly for testing 
    if ($opt_s) {
        print "  - Taking sample of $opt_s families ... ";
        my $cmd = "python $path_to_scripts/subset.py "
                . "-s 0   $file.$analysis_type.vect " 
                . "$opt_s $file.$analysis_type.vect.sample"; 
        system($cmd);
        print "done.\n";                
        # Cleanup - unlink unnecesary files
        unlink "$file.$analysis_type.vect";
        rename "$file.$analysis_type.vect.sample", "$file.$analysis_type.vect";
    }


    # Scales the training data
    print STDERR '  - Scaling training data ... ';
    my $cmd = "$path_to_svm/svm-scale -l -1 -u 1"
           . " -s $file.$analysis_type.vect.range"
           . " $file.$analysis_type.vect "
           . ">$file.$analysis_type.vect.scale";    
    system($cmd);
    print STDERR  "done.\n";
} 
 
 
#-------------------------------------------------------------------------------
# Purpose   : Select best C and g for SVM    
# Usage     : best_Cg($file, $analysis_type, $param_list_reference); 
# Arguments : the file vith oligomer vectors, the type of the analysis, reference to param. list 
# Returns   : best C and g
# Globals   : none
#**********************
sub best_Cg {
    my ($file, $analysis_type, $param_list) = @_;
    my @params = @{$param_list};          
    my @commands = ();    
    
    # Build the SVM model        
    print STDERR "  - Building SVM model ... ";
    # writes all the system calls into an array     
    
    # Command list
    foreach my $i (0..$#params) {
        my $C = 2**$params[$i][0];       
        my $g = 2**$params[$i][1];              
        my $cmd = "$path_to_svm/svm-train -c $C -g $g -m $memory -v 5 " 
                . "$file.$analysis_type.vect.scale "
                . "&>$workdir/output_$C.$g.txt";        
        push(@commands, $cmd);        
    }

    #Dumpvalue->new->dumpValue(\@commands);
    
    # Calculates the number of rounds (it depends on the number of availabe 
    # processors)
    my $rounds;
    if (@commands % $processors == 0 ) {
        $rounds = @commands/$processors;
    }
    elsif (@commands % $processors != 0 ) {
        $rounds = int(@commands/$processors) + 1;    
    }

    
    ### Fires off the threads
    # The number of rounds
    foreach my $k (0..$rounds-1) {
        # The number of threads in the last round, where the number of threads 
        # may be lower than the number of specified processors
        if ($k == $rounds - 1) {
            my $last_round = @commands - ( ($rounds-1) * $processors );
            foreach my $l (0..$last_round-1) {
                my $position_in_array = ($k * $processors) + $l;
                my $thr = threads->new(\&system_call, $commands[$position_in_array]);
            } 
            # Joins them 
            foreach my $threads (threads->list()) {
                $threads->join();
            }    
        }
        # Fires out threads in the other rounds, where the number of processors 
        # equals the specified   
        else {
            foreach my $l (0..$processors-1) {
                my $position_in_array = ($k * $processors) + $l;
                my $thr = threads->new(\&system_call, $commands[$position_in_array]);
            }
            # Joins them 
            foreach my $threads (threads->list()) {
                $threads->join();
            }    
        }
    }
    
    # Processes the output data     
    my %CVs = ();            
    foreach my $i (0..$#params) {
        my $C = 2**$params[$i][0];       
        my $g = 2**$params[$i][1];             
        open (IN, "<", "$workdir/output_$C.$g.txt") or die "can't open IN output_$C.$g.txt $!";
        my @lines = <IN>;
        close IN;            
        unlink "$workdir/output_$C.$g.txt";
        my $cross_validation = (split/\s+/, $lines[-1])[-1]; 
        $cross_validation =~ s/\%//;
        $CVs{$cross_validation} = [ $C, $g ];   
    }
    
    my $best_CV = (sort {$a <=> $b} keys %CVs)[-1];
    my $best_C = log($CVs{$best_CV}[0])/log(2);
    my $best_g = log($CVs{$best_CV}[1])/log(2);
    
    print STDERR "done. (Cross Validation: $best_CV%) \n";
    my $type = (split/\//, $file)[-1];
    print OUT "$type.$analysis_type.vect.scale   Cross Validation= $best_CV%  ", 
              "log2(best_C)= $best_C  log2(best_g)= $best_g\n";
    
    #Dumpvalue->new->dumpValue(\%CVs);
    return $best_C, $best_g;
}

#-------------------------------------------------------------------------------
# Purpose   : runs the final SVM training    
# Usage     : train_selection($file, $analysis_type, $C, $g); 
# Arguments : scalars  
# Returns   : nothing, executes
# Globals   : path vars
#********************** 
sub train_selection {
    my $t1 = time;
    my ($file, $analysis_type, $C, $g) = @_;
    $C = 2**$C;
    $g = 2**$g;
    print STDERR '  - Training ... ';
    my $cmd = "$path_to_svm/svm-train -c $C -g $g " 
            . "$file.$analysis_type.vect.scale " 
            . "$file.$analysis_type.vect.scale.model &>/dev/null";
    system($cmd);
    print STDERR "done.\n\n";
    
    # Cleanup - unlink unnecesary files
    unlink "$file.$analysis_type.vect";        
    unlink "$file.$analysis_type.vect.scale";

    # Moves the results to the right folder
    rename("$file.$analysis_type.vect.scale.model", "$file.$analysis_type.vect.model");
    system("mv $file.$analysis_type.vect.range $path_to_oligomers");              
    system("mv $file.$analysis_type.vect.model $path_to_oligomers");    
}


#-------------------------------------------------------------------------------
# Purpose   : implements the system function as a subroutine    
# Usage     : system_call($cmd); 
# Arguments : string  
# Returns   : nothing, executes
# Globals   : none
#********************** 
sub system_call {
    my $cmd = $_[0];
    system($cmd);
}
