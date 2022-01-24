#!/usr/bin/env perl

# A script to classify unknown transposn consensi, using the SVM, LVQ, RF, ICM 
# databases. The input file is in fasta format.

use warnings;
use lib X  
use TEclass;
use Getopt::Std;
use Configuration;

getopts('c:ho:tr');
our ($opt_c, $opt_h, $opt_o, $opt_t, $opt_r); 


#-------------------------------- Prints help ----------------------------------

if ($opt_h) {
    print STDERR "\nUsage:\n"
               . "perl TEclassTest.pl [flags] <file.fa>\n\n"
               . "flags:\n"
               . " -c '/some/directory' : uses a user provided path to the classifiers.\n" 
               . "      Without the flag the default is the TEclass-2.1/classifiers directory.\n"
               . " -o '/some/directory' : writes the output to a user provided directory.\n" 
               . "      Without the flag a directory is created for the output in the current \n" 
               . "      folder using the input filename and a timestamp.\n"               
               . " -t : uses a decision algorithm that performs best when the\n" 
               . "      SVM models were built only on a subset of the entire\n"
               . "      RepBase.\n"
               . " -r : reverse complements the input sequences that are \n" 
               . "      predicted to be in the reverse orientation\n"
               . " -h : this help\n\n";
    exit;
}


#------------------------- Path and other global variables ---------------------

# Uses other than default classifiers
$path_to_classifiers = $opt_c if ($opt_c);
if (! -d "$path_to_classifiers/ICM_model" 
 || ! -d "$path_to_classifiers/LVQ_codebooks" 
 || ! -d "$path_to_classifiers/RandomForests"
 || ! -d "$path_to_classifiers/tetramer_models"
 || ! -d "$path_to_classifiers/pentamer_models") {
    print STDERR "Some of your classifiers appear to be missing.\n";
    exit;
}
 
# Paths to model folders 
our $path_to_tetramers     = $path_to_classifiers . '/tetramer_models';
our $path_to_pentamers     = $path_to_classifiers . '/pentamer_models';
our $path_to_LVQ_codebooks = $path_to_classifiers . '/LVQ_codebooks';
our $path_to_RFmodels      = $path_to_classifiers . '/RandomForests';
our $path_to_ICM_model     = $path_to_classifiers . '/ICM_model';

# Global hash variables
our %ORFs     = ();  # repeat => ORF coordinates (a..b c..d) 
our %SVM4     = ();  # repeat => forward Retro Non-LTR, LINE
our %SVM5     = ();  # repeat => forward Retro Non-LTR, LINE
our %LVQ      = ();  # repeat => forward LINE
our %RFs      = ();  # repeat => forward Retro Non-LTR, LINE
our %rev_comp = ();  # repeat => 1     - only reverse complemented repeats
our %RESULT   = ();  # repeat => LINE
our %OUT      = ();  # repeat => sequence


#---------------------------------- Preprocessing ------------------------------
#
#  - Linearises and classifies the input sequences according to their length -                              
#
#-------------------------------------------------------------------------------

# No input
if (! defined $ARGV[0]) {
    print STDERR "Please provide input data.\n";
    exit;
}

# Gets out the actual filename from the provided path
my $inputfile;  
if ($ARGV[0] =~ m/\//) {
    $inputfile = (split/\//, $ARGV[0])[-1];
}
else {
    $inputfile = $ARGV[0];
}

# Makes the output directory if provided
our $path_to_outdir;
if ($opt_o) {
    $path_to_outdir = $opt_o;
}
else {
    my $t = time;
    $path_to_outdir = $inputfile . '_' . $t;
}

if (-d $path_to_outdir) {
    my $t = time;
    $path_to_outdir .= '_' . $t;  
}

# Makes the output dir
if (! -d $path_to_outdir) {
    mkdir $path_to_outdir or die "Cannot make output directory $path_to_outdir. \n";  
}

system("cp $ARGV[0] $path_to_outdir");


# Opens temporary files
open( SHORT,  ">", "$path_to_outdir/short.tmp")  or die "Cannot open SHORT  $!";
open( MEDIUM, ">", "$path_to_outdir/medium.tmp") or die "Cannot open MEDIUM $!";
open( LONG,   ">", "$path_to_outdir/long.tmp")   or die "Cannot open LONG   $!";
open( XLONG,  ">", "$path_to_outdir/xlong.tmp")  or die "Cannot open XLONG  $!";

my @filehandles = ( *SHORT, *MEDIUM, *LONG, *XLONG );

# In case the sequences are wrapped it linearises them
linearize("$inputfile");

my $name;
my $sequence;
# Splits the uploaded file int four files, according to their length
# (short, medium, long, xlong)
print STDERR "\nBuilding .tmp files from .fa\n";
open( UPLOADFILE, "<", "$path_to_outdir/$inputfile" ) or die "UPLOADFILE $!";
while (<UPLOADFILE>) {
    chomp;
    s/\r//;
    # Identifies a new repeat
    if ( $_ =~ m/^>/ ) {
        if ( defined($sequence) ) {
            classify_length( $name, $sequence, @filehandles );
            $sequence = ();
        }
        $name = $_;
    }
    else {
        $sequence .= uc($_);
    }
}
classify_length( $name, $sequence, @filehandles );
close UPLOADFILE;
close SHORT;
close MEDIUM;
close LONG;
close XLONG;


# The list of files to be processed (with the right extension)
our @files = qw/ short medium long xlong /;


#-------------------------------------- MAIN -----------------------------------
#
#                               - Tests the repeats -
#      
#-------------------------------------------------------------------------------


#------------------------ Forward vs. Reverse classification -------------------

{
    print STDERR "\n1.Determining repeat orientation (forward vs. reverse)\n";

    # 1st SVM classification (pentamers)
    print STDERR "\n  First SVM classification:\n";    
    build_vectors(5);
    make_classification( \&test_with_SVM, 'forward_vs_reverse');

    # 2nd SVM classification (tetramers)
    print STDERR "\n  Second SVM classification:\n";
    build_vectors(4);
    make_classification( \&test_with_SVM, 'forward_vs_reverse');

    # LVQ classification
    print STDERR "\n  LVQ classification:\n"; 
    make_classification( \&test_with_LVQ, 'forward_vs_reverse');
    
    # RF classification
    print STDERR "\n  RF classification:\n";    
    make_classification( \&test_with_RF,  'forward_vs_reverse');

    # Decides about the orientation of the sequences 
    forward_vs_reverse();
}

#--------------------------- Repeat type classification ------------------------

{
    print STDERR "\n\n2. Determining the repeat type. \n";

    # Determines the ORFs
    print STDERR "  Search for ORFs ... ";
    search_for_ORFs();   # this is to be changed 
    print STDERR "done\n";

    # 1st SVM classification (pentamers)
    print STDERR "\n  First SVM classification:\n";
    build_vectors(5);
    make_classification( \&test_with_SVM, 'repeats'); 

    # 2nd SVM classification (tetramers)
    print STDERR "\n  Second SVM classification:\n";
    build_vectors(4);
    make_classification( \&test_with_SVM, 'repeats');

    # LVQ classification
    print STDERR "\n  LVQ classification:\n"; 
    make_classification( \&test_with_LVQ, 'repeats');

    # RF classification
    print STDERR "\n  RF classification:\n";
    make_classification( \&test_with_RF,  'repeats');

    # Prints the result files
    repeat_type();
    write_lib_and_stat_files();
    write_html_file();
    
    # Deletes the unnecesary files
    cleanup();
    print STDERR "\nClassification ready.\n"
}


#----------------------------------- SUBROUTINES -------------------------------
#
#    make_classification
#    linearize
#    test_with_SVM
#    test_with_RF
#    test_with_LVQ
#    search_for_ORFs
#    make_vectors
#    standardise_vectors
#    standardise_vectors_SVM
#    LVQ_split
#    forward_vs_reverse
#    repeat_type
#    build_vectors
#    write_lib_and_stat_files
#    write_html_file
#    decide
#    decide_three
#    cleanup
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Purpose   : "Walks" through the decision making tree; first decides whether
#             the repeat is in forward or reverse order, next whether DNA
#             or Retroelement, etc.
# Usage     : make_classification($subroutine-reference, $string);
# Arguments : a reference to the subroutine that will do the classification, and 
#             a string defining the classification type
# Returns   : nothing, prints
# Globals   : @files
#**********************
sub make_classification {
    my ($decide_whether, $classification_type) = @_;
    foreach my $i (0..$#files) {
        # The file exists
        if (-s "$path_to_outdir/$files[$i].tmp" 
         && $classification_type eq 'forward_vs_reverse') {
            $decide_whether->( 'forward', 'reverse', $files[$i] );
        }
        elsif (-s "$path_to_outdir/$files[$i].tmp") {
            # Makes the binary classifications
            $decide_whether->( 'DNA', 'Retro',  $files[$i] );
            $decide_whether->( 'LTR', 'nonLTR', $files[$i] );
            # There are no long or xlong SINEs
            if ($files[$i] =~ m/short|medium/) {
                $decide_whether->( 'LINE', 'SINE', $files[$i] );
            }
        }
    }
}


#-------------------------------------------------------------------------------
# Purpose   : linearizes a fasta file
# Usage     : linearize($filename);
# Arguments : string
# Returns   : prints to a new file
# Globals   : path variables
#******************
sub linearize {
    open( IN,  "<", "$path_to_outdir/$_[0]" )     or die "cant open IN $!";
    open( LIN, ">", "$path_to_outdir/$_[0].lin" ) or die "cant open LIN $!";
    my $I      = 0;
    my $String = '';
    while (<IN>) {
        if ( $_ =~ m/^>/ ) {
            s/\r//;
            if ( $I == 1 ) {
                print LIN $String, "\n";
            }
            print LIN $_;
            $String = '';
        }
        else {
            $I = 1;
            chomp;
            s/\r//;
            $String .= $_;
        }
    }
    print LIN $String, "\n";

    close IN;
    close LIN;

    unlink "$path_to_outdir/$_[0]";
    rename "$path_to_outdir/$_[0].lin", "$path_to_outdir/$_[0]";
}


#-------------------------------------------------------------------------------
# Purpose   : Classifies the repeat (using libsvm) into one of two categories.
#             Makes several system calls.
# Usage     : &test_with_SVM($string1, $string2, $string3);
# Arguments : category1=$string1, category2=string2, sequence-length(e.g.
#             short, medium, etc.)=string3, normally the $_ of the @files array
# Returns   : nothing, prints
# Globals   : path variables
#******************
sub test_with_SVM {
    my ($repeat_type1, $repeat_type2, $file) = @_;
    my $analysis_type = $repeat_type1 . '_vs_' . $repeat_type2;            
    
    print STDERR "  $file $analysis_type ... ";
    
    # Decides about the variables used (tetra vs pentamers)
    my $path_to_oligomers;
    if ($TEclass::mer == 5) {
        $path_to_oligomers = $path_to_pentamers;
        $SVM_ref = \%SVM5;
    }
    elsif ($TEclass::mer == 4) {
        $path_to_oligomers = $path_to_tetramers;
        $SVM_ref = \%SVM4;
    }
    
    standardise_vectors_SVM($file, $analysis_type, $path_to_oligomers);
    
    # Makes the prediction
    my $cmd = "$path_to_svm/svm-predict " 
         . "$path_to_outdir/$file.$analysis_type.scaled " 
         . "$path_to_oligomers/$file.$analysis_type.vect.model "
         . "$path_to_outdir/$file.$analysis_type.predict " 
         . ">/dev/null";
    system( $cmd );
    
    # Writes the result of the classification to the classification hash
    open( IN1, "<", "$path_to_outdir/$file.tmp") or die "cannot open IN1 $!";
    open( IN2, "<", "$path_to_outdir/$file.$analysis_type.predict") 
                                                    or die "cannot open IN2 $!";
    # Makes the classification and writes the result to the SVM hashes
    while (1) {
        my $repeat = <IN1>;
        my $seq    = <IN1>;
        my $pred   = <IN2>;
        chomp $repeat;
        chomp $seq; # this variable is not used anywhere
        chomp $pred;
        # Prevents classifying known, DNA and LTR transposon later in the cycle
        if ( defined(${$SVM_ref}{$repeat}) 
        &&  ${$SVM_ref}{$repeat}[-1] =~ m/DNA|\bLTR/ ) {        
            # Do nothing
        }
        elsif ( $pred eq '-1' ) {
            push @{${$SVM_ref}{$repeat}}, $repeat_type1;
        }
        elsif ( $pred eq '1' ) {
            push @{${$SVM_ref}{$repeat}}, $repeat_type2;
            # There are no SINEs in the long and xlong categories
            if ( $file =~ m/long|xlong/ && $analysis_type eq 'LTR_vs_nonLTR' ) {
                push @{${$SVM_ref}{$repeat}}, 'LINE';
            }            
        }
        last if eof(IN1);
    }
    close IN1;
    close IN2;
    print STDERR "done.\n";    
}


#-------------------------------------------------------------------------------
# Purpose   : Tests the repeats using the RF models
# Usage     : &test_with_RF($string1, $string2, $string3);
# Arguments : category1=$string1, category2=string2, sequence-lenght(e.g.
#             short, medium, etc.)=string3, normally the $_ of the @files array
# Returns   : nothing, prints
# Globals   : path variables, %RFs
#******************
sub test_with_RF {
    my ($repeat_type1, $repeat_type2, $file) = @_;
    my $analysis_type = $repeat_type1 . '_vs_' . $repeat_type2;
            
    print STDERR "  $file $analysis_type ... ";    
    standardise_vectors($file, $analysis_type, $path_to_RFmodels);
    translate_from_SVM_to_CSV($file, $analysis_type, $path_to_outdir);
    
    # Makes the classification
    my $cmd =   "$path_to_librf/rfpredict "
           . "-m $path_to_RFmodels/$file.$analysis_type.rf.model "
           . "-d $path_to_outdir/$file.$analysis_type.csv "
           . "-f 256 "
           . "-o $path_to_outdir/$file.$analysis_type.rf.out "
           . "-l $path_to_outdir/$file.$analysis_type.lab "
           . " 1>$path_to_outdir/$file.$analysis_type.rf.log 2>/dev/null";
    system($cmd);

    # Writes the result of the classification to the classification hash
    open( IN1, "<", "$path_to_outdir/$file.tmp" ) or die "cannot open IN1 $!";
    open( IN2, "<", "$path_to_outdir/$file.$analysis_type.rf.out") or die "cannot open IN2 $!";
    while (1) {
        my $repeat = <IN1>;
        my $seq    = <IN1>;
        my $pred   = <IN2>;
        chomp $repeat;
        chomp $seq; # this variable is not used anywhere
        chomp $pred;
        # Prevents classifying known, DNA and LTR transposon later in the cycle
        if ( defined($RFs{$repeat}) 
        &&  @{$RFs{$repeat}}[-1] =~ m/DNA|\bLTR/ ) {
            # Do nothing
        }
        elsif ( $pred <  0.5 ) {
            push @{$RFs{$repeat}}, $repeat_type1;
        }
        elsif ( $pred >= 0.5 ) {
            push @{$RFs{$repeat}}, $repeat_type2;
            # There are no SINEs in the long and xlong categories
            if ( $file =~ m/long|xlong/ && $analysis_type eq 'LTR_vs_nonLTR' ) {
                push @{$RFs{$repeat}}, 'LINE';
            }            
        }
        last if eof(IN1);
    }
    close IN1;
    close IN2;
    print STDERR "done.\n";
}


#-------------------------------------------------------------------------------
# Purpose   : Tests the repeats using the LVQ codebooks
# Usage     : &test_with_LVQ($string1, $string2, $string3);
# Arguments : category1=$string1, category2=string2, sequence-lenght(e.g.
#             short, medium, etc.)=string3, normally the $_ of the @files array
# Returns   : nothing, prints
# Globals   : path variables
#******************
sub test_with_LVQ {
    my ($repeat_type1, $repeat_type2, $file) = @_;
    
    # Skips the unnecessary rounds
    if ($repeat_type1 eq 'forward' || $repeat_type1 eq 'DNA') {     
        my $analysis_type = &LVQ_split($repeat_type1);
        
        print STDERR "  $file $analysis_type ... ";
        
        standardise_vectors($file, $analysis_type, $path_to_LVQ_codebooks);
        translate_from_SVM_to_LVQ($file, $analysis_type, 256, $path_to_outdir);        
        
        rename "$path_to_outdir/$file.$analysis_type.codebook", 
               "$path_to_outdir/$file.$analysis_type.data";
        
        # Makes the classification
        my $cmd = "$path_to_LVQ/classify " 
                 . "-din  $path_to_outdir/$file.$analysis_type.data " 
                 . "-cin  $path_to_LVQ_codebooks/$file.$analysis_type.codebook " 
                 . "-dout $path_to_outdir/$file.$analysis_type.lvq "
                 . "1>/dev/null 2>/dev/null";
        system($cmd);
    
        # Writes the result of the classification to the classification hash
        open( IN1, "<", "$path_to_outdir/$file.tmp" ) or die "cannot open IN1 $!";
        open( IN2, "<", "$path_to_outdir/$file.$analysis_type.lvq") or die "cannot open IN2 $!";
        my $pred   = <IN2>; # This is to skip the first line containing '256'
        while (1) {
            my $repeat = <IN1>;
            my $seq    = <IN1>;
            $pred      = <IN2>;
            my @result = split/\s+/, $pred; 
            chomp $repeat;
            chomp $seq; # this variable is not used anywhere
            chomp $pred;
            # Prevents classifying known, DNA and LTR transposon later in the cycle
            if (defined($LVQ{$repeat}) 
            && @{$LVQ{$repeat}}[-1] =~ m/DNA|LTR|LINE|SINE/) {
                # Do nothing
            }
            else {
                my @output;
                if ($result[-1] =~ m/DNA/) {
                    @output = ( 'DNA' );
                }
                elsif ($result[-1] =~ m/LTR/) {
                    @output = ( 'Retro', 'LTR' );
                }
                elsif ($result[-1] =~ m/LINE/) {
                    @output = ( 'Retro', 'nonLTR', 'LINE' );
                }
                elsif ($result[-1] =~ m/SINE/) {
                    @output = ( 'Retro', 'nonLTR', 'SINE' );
                }
                else {
                    @output = ( "$result[-1]", );
                }
                push @{$LVQ{$repeat}}, @output;
            }        
            last if eof(IN1);
        }
        close IN1;
        close IN2;
        print STDERR "done.\n";
    }
}


#-------------------------------------------------------------------------------
# Purpose   : Determine the position of ORFs if there are any 
# Usage     : &test_for_ORFs();
# Arguments : no arguments
# Returns   : The %ORFs hash, and prints
# Globals   : path and hash variables, 
#******************
sub search_for_ORFs {
    # The processed file
    my $cmd = "$path_to_glimmer/glimmer3 -A atg -l -X "
            . "$path_to_outdir/$inputfile.rev " 
            . "$path_to_ICM_model/RepBase.icm " 
            . "$path_to_outdir/$inputfile.rev "
            . "1>/dev/null 2>/dev/null";        
    system($cmd);
    
    my $repeat; 
    open (IN, "<", "$path_to_outdir/$inputfile.rev.predict") or die "cant open IN $!";
    # Reads in the content of the ORFs into a hash
    while (<IN>) {
        if ($_ =~ m/^>/) {
            chomp;              
            $repeat = $_;
        }     
        # Adds one ORF to the hash 
        else {
            my @coords = split;
            if ($coords[-1] == 3) {
                my $one_ORF = $coords[1] . '..' . $coords[2] . ':' . $coords[3];
                push @{$ORFs{$repeat}}, $one_ORF;
            }
        }
    }     
    close IN;
}


#-------------------------------------------------------------------------------
# Purpose   : Makes oligomer frequency vectors from linear DNA sequences   
# Usage     : make_vectors(file tag, analysis type)
# Arguments : $string, $string
# Returns   : nothing, prints
# Globals   : none
#******************
sub make_vectors {
    my ($file) = @_;

    open(IN,  "<", "$path_to_outdir/$file.tmp")  or die "cannot open IN  $!";
    open(OUT, ">", "$path_to_outdir/$file.vect") or die "cannot open OUT $!";
    my $outfile = *OUT;
 
    # Processes the linearised .lib file
    while (1) {
        my $repeat = <IN>;
        $repeat =~ s/>//;
        chomp($repeat); 
        my $sequence = <IN>;
        process_sequence($sequence, 0, $outfile);                
        last if eof(IN);
    }
    close OUT; 
    close IN;
}


#-------------------------------------------------------------------------------
# Purpose   : Standardises the vector file, so that the different oligomer 
#             frequencies will get te same weight              
# Usage     : &standardise_vectors(file, analysis type, rangefile path);
# Arguments : $string, $string, $string 
# Returns   : nothing, prints   
# Globals   : none
#**********************
sub standardise_vectors {
    my ($file, $analysis_type, $path_to_rangefile) = @_; 
    my $cmd = "$path_to_svm/svm-scale -r" 
          . "  $path_to_rangefile/$file.$analysis_type.range" 
          . "  $path_to_outdir/$file.vect"
          . " >$path_to_outdir/$file.$analysis_type.scaled";
    system($cmd);   
}


#-------------------------------------------------------------------------------
# Purpose   : Standardises the vector file, so that the different oligomer 
#             frequencies will get te same weight              
# Usage     : &standardise_vectors(file, analysis type, rangefile path);
# Arguments : $string, $string, $string 
# Returns   : nothing, prints   
# Globals   : none
#******************
sub standardise_vectors_SVM {
    my ($file, $analysis_type, $path_to_rangefile) = @_;     
    my $cmd = "$path_to_svm/svm-scale -r" 
          . "  $path_to_rangefile/$file.$analysis_type.vect.range" 
          . "  $path_to_outdir/$file.vect"
          . " >$path_to_outdir/$file.$analysis_type.scaled";
    system($cmd);   
}


#-------------------------------------------------------------------------------
# Purpose   : Decide which codebooks to use. 
# Usage     : &LVQ_split(repeat_type)
# Arguments : $string,
# Returns   : a string, determining the analysis type
# Globals   : none
#******************
sub LVQ_split {
    my ($arg1) = @_;
    my $analysis_type;
    if ($arg1 eq 'forward') {
        $analysis_type = 'forward_vs_reverse';
    }
    else {
        $analysis_type = 'repeats';
    }
    return $analysis_type;
}


#-------------------------------------------------------------------------------
# Purpose   : Decides whehter the sequence is in the forward or reverse direction 
# Usage     : &forward_vs_reverse()
# Arguments : no arguments
# Returns   : nothing, prints
# Globals   : @files, path variables, and the main hash variables
#******************
sub forward_vs_reverse {
    open (GLM, ">", "$path_to_outdir/$inputfile.rev") or die "can't open GLM $!";
    foreach (@files) {  
        open (IN,  "<", "$path_to_outdir/$_.tmp")  or die "can't open IN  $!";
        open (OUT, ">", "$path_to_outdir/$_.tmp2") or die "can't open OUT $!";         
        # Rads in the repeats from the tmp files
        while (1) {
            last if eof(IN); 
            my $decision = 0;
            my $reverse_complemented;
            my $repeat = <IN>;
            my $seq = <IN>;
            chomp $repeat;
            chomp $seq;
            # Decides, based on majority.
            if ($SVM4{$repeat}[0] eq 'forward') { $decision++; }
            if ($SVM5{$repeat}[0] eq 'forward') { $decision++; }
            if ($LVQ{$repeat}[0]  eq 'forward') { $decision++; }
            if ($RFs{$repeat}[0]  eq 'forward') { $decision++; }
            
            # Reverse complement the sequence if necessary
            if ($decision <= 1) {
                $reverse_complemented = reverse_complement($seq);
                print OUT "$repeat\n$reverse_complemented\n";

                if ($opt_r) {
                    $OUT{$repeat} = $reverse_complemented; 
                    print GLM "$repeat\n$reverse_complemented\n";
                }
                else {
                    $OUT{$repeat} = $seq;
                    print GLM "$repeat\n$seq\n";                
                }
                
                $rev_comp{$repeat} = 1;
            }
            else {
                $OUT{$repeat} = $seq;
                print OUT "$repeat\n$seq\n";
                print GLM "$repeat\n$seq\n";                
            }            
        }
        close IN;
        close OUT;
        rename "$path_to_outdir/$_.tmp2", "$path_to_outdir/$_.tmp"; 
    }
    close GLM;
}


#-------------------------------------------------------------------------------
# Purpose   : to make the final classification of the repeats 
# Usage     : &repeat_type()
# Arguments : no arguments
# Returns   : nothing, prints
# Globals   : @files, path and hash variables,
#******************
sub repeat_type {
    foreach (@files) {
        open (IN,  "<", "$path_to_outdir/$_.tmp") or die "can't open IN  $!";         
        # Reads in the repeats from the tmp files
        while (1) {
            last if eof(IN); 
            my $result = '';
            my $repeat = <IN>;
            chomp $repeat;
            my $seq = <IN>;
            chomp $seq;  
            # Chooses about the final decision algorithm
            if ($opt_t) {
                $RESULT{$repeat} = decide_three($repeat);
            }
            else {
                $RESULT{$repeat} = decide($repeat);
            }          
        }
        close IN; 
    }
}


#-------------------------------------------------------------------------------
# Purpose   : To build oligomer vectors 
# Usage     : &build_vectors($int)
# Arguments : integer (4 or 5)
# Returns   : the global $xmer variable and prints the vectors
# Globals   : @files, path variables,
#******************
sub build_vectors {
    my ($oligomer) = @_;
    print STDERR "  - Making vectors ... ";
    
    if ($oligomer == 4) {
        build_tetramer_hash();
    }
    elsif ($oligomer == 5) {
        build_pentamer_hash();
    }
    
    foreach my $i (0..$#files) {
        if (-s "$path_to_outdir/$files[$i].tmp" ) {
            make_vectors($files[$i]);
        }
    }
    print STDERR "done.\n";
}



#-------------------------------------------------------------------------------
# Purpose   : Creates the final lib file with the result of classification 
# Usage     : &write_lib_file();
# Arguments : none
# Returns   : nothing, prints
# Globals   : path and hash variables
#******************
sub write_lib_and_stat_files {
    open (OUT, ">", "$path_to_outdir/$inputfile.lib") or die "Can't open OUT $!";
    open (STA, ">", "$path_to_outdir/$inputfile.stat") or die "Can't open STA $!";
    
    # The number of different repeat types 
    my $DNAs    = 0;
    my $LTRs    = 0;
    my $LINEs   = 0;
    my $SINEs   = 0;
    my $unknown = 0;

    # Prints out the header
    foreach my $key (sort keys %OUT) {
        print OUT "$key|TEclass result: $RESULT{$key}";
        # Prints info on reverse complementation
        if (defined $rev_comp{$key} && defined $opt_r) {
            print OUT "|reverse complemented";
        }   
        elsif (defined $opt_r) {
            print OUT "|forward";
        }             
        elsif (defined $rev_comp{$key}) {
            print OUT "|reverse";
        }         
        else {
            print OUT "|forward";
        }         
                
        # Prints out the ORFs if there are any 
        if ( defined $ORFs{$key} && @{$ORFs{$key}} ) {
            print OUT "|ORFs: ";
            my $orf = join " ", @{$ORFs{$key}};            
            print OUT "$orf";
        }        
        print OUT "\n";
        
        # Prints out the sequence
        print OUT "$OUT{$key}\n"; 
        
        # Makes repeat statistics
        if    ( $RESULT{$key} eq 'DNA'  ) { $DNAs++;  }
        elsif ( $RESULT{$key} eq 'LTR'  ) { $LTRs++;  }
        elsif ( $RESULT{$key} eq 'LINE' ) { $LINEs++; }
        elsif ( $RESULT{$key} eq 'SINE' ) { $SINEs++; }
        else                              { $unknown++; }        
    }

    # Prints repeat statistics
    print STA "\nRepeat statistics:"
            . "\nDNA transposons: $DNAs"
            . "\nLTRs:  $LTRs"
            . "\nLINEs: $LINEs"
            . "\nSINEs: $SINEs"
            . "\nUnclear:  $unknown"
            . "\n\nTotal: ", ($DNAs+$LTRs+$LINEs+$SINEs+$unknown), "\n" ;
    
    # Prints repeat statistics
    print STDERR "\nRepeat statistics:"
               . "\nDNA transposons: $DNAs"
               . "\nLTRs:  $LTRs"
               . "\nLINEs: $LINEs"
               . "\nSINEs: $SINEs"
               . "\nUnclear:  $unknown"
               . "\n\nTotal: ", ($DNAs+$LTRs+$LINEs+$SINEs+$unknown), "\n" ;               
    close OUT;
    close STA;
}


#-------------------------------------------------------------------------------
# Purpose   : Creates the final html file with the details of classification 
# Usage     : &write_html_file();
# Arguments : none
# Returns   : nothing, prints
# Globals   : path and hash variables
#******************
sub write_html_file {
    # The outputfile with the summary of the classification procedure
    open( HTML, ">", "$path_to_outdir/$inputfile.html") or die "can't open HTML $!";
    
    # Prints header
    print HTML "<html>\n<body>\n<table>\n";
    print HTML "<tr bgcolor='#EEEEEE'> " 
             . "<td width=50px><b>Nr.</b></td> " 
             . "<td width=240px><b>ID</b></td> "   
             . "<td width=120px><b>Result</b></td> "
             . "<td width=200px><b>ORFs  </b></td> "                
             . "<td width=120px><b>Forw./Rev.</b></td> "
             . "<td width=300px><b>SVM classification (4mer)</b></td> "
             . "<td width=300px><b>SVM classification (5mer)</b></td> "
             . "<td width=300px><b>LVQ classification</b></td> "                                
             . "<td width=300px><b>RF classification</b></td> " 
             . "</tr>\n";
    
    # Prints a table row
    my $TE_nr = 0;
    foreach my $repeat (sort keys %OUT) {
        my $repeat_name = $repeat;
        $repeat_name =~ s/>//;
        $TE_nr++;
        # The background color is gray
        if ( $TE_nr / 2 == int( $TE_nr / 2 ) ) {
            print HTML "<tr bgcolor='#EEEEEE'> ";
        }
        # The background color is white
        else {
            print HTML "<tr                  > ";
        }
        # Prints the repeat name and final classification result
        print HTML "<td width=50px>$TE_nr</td> " 
                 . "<td width=240px>$repeat_name</td> "    
                 . "<td width=120px>$RESULT{$repeat}</td> ";
                 
        # ORFs, if detected         
        if ( defined $ORFs{$repeat} && @{$ORFs{$repeat}} ) {         
            my $orfs = join ' ', @{$ORFs{$repeat}};         
            print HTML "<td width=200px>$orfs</td> ";
        }
        else {
            print HTML "<td width=200px> - </td> ";
        }
        
        # Reverse complemented   
        if ( defined $rev_comp{$repeat} ) {                  
            if ($opt_r) {
                print HTML "<td width=120px>rev. compl.</td> ";                
            }
            else {              
                print HTML "<td width=120px>reverse</td> ";
            }
        }
        else {
            print HTML "<td width=200px>forward</td> ";
        }
        
        # Prints clasification details
        my $svm4_string = join " ", @{$SVM4{$repeat}};
        my $svm5_string = join " ", @{$SVM5{$repeat}};         
        my $lvq_string = join " ", @{$LVQ{$repeat}};        
        my $rf_string  = join " ", @{$RFs{$repeat}};                        
        print HTML "<td width=300px>$svm4_string</td> "
                 . "<td width=300px>$svm5_string</td> "
                 . "<td width=300px>$lvq_string</td> "         
                 . "<td width=300px>$rf_string </td> "                                
                 . "</tr>\n";
    }
    print HTML "</table>\n</body>\n</html>\n";
    close HTML;
}


#-------------------------------------------------------------------------------
# Purpose   : It makes the decision about the repeat type. 
# Usage     : decide($string)
# Arguments : the name of the repeat
# Returns   : string - the result of the classification
# Globals   : global hash variables
#***********************************
sub decide {
    my ($repeat) = @_;
    my $result = '';
    # Decides about the result.
    # Pentamer SVM and tetramer SVM result is the same  
    if ($SVM5{$repeat}[-1] eq $SVM4{$repeat}[-1]) {
        $result = $SVM5{$repeat}[-1]; 
    }
    # Pentamer SVM and any other method is the same
    elsif ($SVM5{$repeat}[-1] eq $LVQ{$repeat}[-1] 
        && $SVM5{$repeat}[-1] eq $RFs{$repeat}[-1]) {
        $result = $SVM5{$repeat}[-1]; 
    }
    # Any three are the same
    else {
        my @rounds = qw/ -1 2 1 /;
        # If any three are similar
        foreach my $i (@rounds) {
            if (   defined $SVM5{$repeat}[$i]
                && defined $SVM4{$repeat}[$i]
                && defined $LVQ{$repeat}[$i]
                && $SVM5{$repeat}[$i] eq $SVM4{$repeat}[$i] 
                && $SVM5{$repeat}[$i] eq $LVQ{$repeat}[$i]) {
                $result = $SVM5{$repeat}[$i];
                return $result;    
            }    
            elsif (defined $SVM5{$repeat}[$i]
                && defined $SVM4{$repeat}[$i]
                && defined $RFs{$repeat}[$i]                        
                && $SVM5{$repeat}[$i] eq $SVM4{$repeat}[$i] 
                && $SVM5{$repeat}[$i] eq $RFs{$repeat}[$i]) {
                $result = $SVM5{$repeat}[$i];
                return $result;    
            }    
            elsif (defined $SVM5{$repeat}[$i]
                && defined $LVQ{$repeat}[$i]
                && defined $RFs{$repeat}[$i]        
                && $SVM5{$repeat}[$i] eq $LVQ{$repeat}[$i] 
                && $SVM5{$repeat}[$i] eq $RFs{$repeat}[$i]) {
                $result = $SVM5{$repeat}[$i];
                return $result;    
            }    
            elsif (defined $SVM4{$repeat}[$i]
                && defined $LVQ{$repeat}[$i]
                && defined $RFs{$repeat}[$i]                        
                && $SVM4{$repeat}[$i] eq $LVQ{$repeat}[$i] 
                && $SVM4{$repeat}[$i] eq $RFs{$repeat}[$i]) {
                $result = $SVM4{$repeat}[$i];
                return $result;    
            }   
        }
        $result = 'unclear';
    }
    return $result;
}

#-------------------------------------------------------------------------------
# Purpose   : An alternative subroutine to make the decision about the repeat 
#             type. It takes the one that is supported by at least three methods. 
# Usage     : decide_three($string)
# Arguments : the name of the repeat
# Returns   : string - the result of the classification
# Globals   : global hash variables
#***********************************
sub decide_three {
    my ($repeat) = @_;
    my $result = '';
    my @rounds = qw/ -1 2 1 /;
    # If any three are similar
    foreach my $i (@rounds) {
        if (   defined $SVM5{$repeat}[$i]
            && defined $SVM4{$repeat}[$i]
            && defined $LVQ{$repeat}[$i]
            && $SVM5{$repeat}[$i] eq $SVM4{$repeat}[$i] 
            && $SVM5{$repeat}[$i] eq $LVQ{$repeat}[$i]) {
            $result = $SVM5{$repeat}[$i];
            return $result;    
        }    
        elsif (defined $SVM5{$repeat}[$i]
            && defined $SVM4{$repeat}[$i]
            && defined $RFs{$repeat}[$i]                        
            && $SVM5{$repeat}[$i] eq $SVM4{$repeat}[$i] 
            && $SVM5{$repeat}[$i] eq $RFs{$repeat}[$i]) {
            $result = $SVM5{$repeat}[$i];
            return $result;    
        }    
        elsif (defined $SVM5{$repeat}[$i]
            && defined $LVQ{$repeat}[$i]
            && defined $RFs{$repeat}[$i]        
            && $SVM5{$repeat}[$i] eq $LVQ{$repeat}[$i] 
            && $SVM5{$repeat}[$i] eq $RFs{$repeat}[$i]) {
            $result = $SVM5{$repeat}[$i];
            return $result;    
        }    
        elsif (defined $SVM4{$repeat}[$i]
            && defined $LVQ{$repeat}[$i]
            && defined $RFs{$repeat}[$i]                        
            && $SVM4{$repeat}[$i] eq $LVQ{$repeat}[$i] 
            && $SVM4{$repeat}[$i] eq $RFs{$repeat}[$i]) {
            $result = $SVM4{$repeat}[$i];
            return $result;    
        }   
    }
    $result = 'unclear';
    return $result;
}

#-------------------------------------------------------------------------------
# Purpose   : Deletes the unnecesary files 
# Usage     : cleanup();
# Arguments : none
# Returns   : nothing
# Globals   : path variable
#******************
sub cleanup { 
    unlink <$path_to_outdir/*.tmp>;
    unlink <$path_to_outdir/*.lab>;
    unlink <$path_to_outdir/*.log>;    
    unlink <$path_to_outdir/*.data>;    
    unlink <$path_to_outdir/*.detail>;
    unlink <$path_to_outdir/*.predict>;    
    unlink <$path_to_outdir/*.lvq>;    
    unlink <$path_to_outdir/*.csv>;    
    unlink <$path_to_outdir/*.out>;       
    unlink <$path_to_outdir/*.scaled>;    
    unlink <$path_to_outdir/*.tmp>;
    unlink <$path_to_outdir/*.vect>;
    unlink <$path_to_outdir/*.rev>;    
}


#------------------------------------- POD -------------------------------------
#
#-------------------------------------------------------------------------------

=head1 NAME

   TEclassTest.pl: assigns a transposon type (DNA, LTR, LINE, SINE)
   to a submitted sequence. 
 
=head1 VERSIONS 

   v1.0 - 14. Nov. 2008.
   v2.0 - 18. Jul. 2009. 
   v2.1 - 27. Jan. 2011.  

=head1 ARGUMENTS AND USAGE

    ./TEclassTest.pl [flags] file.fa
    file.fa is a file with the sequences to test, in fasta format.  

=head1 FLAGS (optional)
  
   -h help
   -c '/some/directory' : uses a user provided path to the classifiers. 
      Without the flag the defalut is the TEclass-2.1/classifiers 
      directory.
   -o '/some/directory' : writes the output to a user provided directory. 
      Without the flag a directory is created for the output in the 
      current folder using the input filename and a timestamp.
   -t Uses a decision algorithm that performs better when the SVM 
      models were built with a realtively small subset of RepBase 
      (~ -s 300-500).
   -r Reverse complements the input sequence if it is predicted to 
      be on the reverse strand
   
=head1 OUTPUT

   It creates an output folder, which contains four files: the original 
   processed file, a *.lib file with the sequences and the result of the 
   classification, an *.html file with information on the classification 
   procedure and a *.stat file with the basic classification statistics.
  
=head1 DEPENDENCIES 

   To run it, you need the following programs installed (besides perl):

   libsvm 
   librf 
   lvq_pak
   glimmer3.

   It also needs several pre-built classifiers in the following 
   directories:
   
   tetramer_models 
   pentamer_models 
   LVQ_codebooks
   RandomForests
   ICM_models
   
   These models have to be built with the TEclass_build.pl script, or 
   can be downloaded from the site http://www.compgen.uni-muenster.de. 

=head1 CONFIGURATION 

   See README for the details of building the models and configuring 
   the script. Should run with Perl 5.8 or higher    
   
=head1 AUTHOR, ETC
   
   Gyorgy Abrusan, (gyorgy.abrusan@gmail.com)

=cut
