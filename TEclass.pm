#!/usr/bin/env perl

package TEclass;
require Exporter;
use warnings;

our @ISA       = qw / Exporter /;
our @EXPORT    = qw / choose_oligomer
                      classify_length
                      classify_type
                      process_sequence
                      reverse_complement
                      build_trimer_hash
                      build_tetramer_hash
                      build_pentamer_hash  
                      standardise_vectors_train                      
                      translate_from_SVM_to_LVQ
                      translate_from_SVM_to_CSV
                      build_repeat_vectors 
                      makedir /;
               
#-------------------------------------------------------------------------------        
#  Subroutines used by the following scripts:
#
#                              Build_Random_forests.pl
#                              Build_LVQ_sets.pl
#                              Build_SVM_models.pl
#                              Build_ICM_models 
#                              TEclassTest.pl
#                              TEclassBuild.pl  
#                              Preprocessing.pl
#                              Configure.pl                    
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------    
#                               - Subroutines -   
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Purpose   : decides whether tri-, tetra-, or pentamers are used, and chooses 
#             one of them.
# Usage     : choose_oligomer($opt_o) 
# Arguments : the oligomer variable of the Getopts::Std.
# Returns   : the subroutines it calls (e.g. build_trimer_hash, see below) 
#             return values, but not this subroutine.
# Globals   : none
#******************
sub choose_oligomer {
    my $oligo = $_[0];
    if (defined ($oligo)) {  
        if ($oligo == 3) {
            build_trimer_hash();
        }
        elsif ($oligo == 5) {
            build_pentamer_hash();
        }
        elsif ($oligo == 4) {
            build_tetramer_hash();
        }
        else {
            print "You have to use 3-, 4- or 5-mers (e.g. -o 4)\n";
            exit;        
        }
    }
    else {         
        print STDERR "You need to specify the oligomer (3-, 4- or 5-mer)\n";
        exit;     
    }
}

#-------------------------------------------------------------------------------
# Purpose   : Writes sequences to different files (.lib format), depending on 
#             their lenght  
# Usage     : classify_length($sequence-name, $sequence, @filehandle-list)
#             sequence name should contain '>'
# Arguments : two strings and a list, the list can have maximum four elements
# Returns   : nothing, prints
# Globals   : none
#******************
sub classify_length {
    my $description = $_[0];
    my $seq = $_[1];
    
    if (length($seq) < 50) {
        return;
    } 
    elsif ( length($seq) < 620 ) { 
        my $FILEHANDLE = $_[2]; 
        print $FILEHANDLE "$description\n", "$seq\n";            
    }
    elsif ( length $seq < 1800 ) {
        my $FILEHANDLE = $_[3];
        print $FILEHANDLE "$description\n", "$seq\n";                     
    }
    elsif ( length $seq < 4000 ) {
        my $FILEHANDLE = $_[4];
        print $FILEHANDLE "$description\n", "$seq\n";             
    }
    else {
        my $FILEHANDLE = $_[5];
        print $FILEHANDLE "$description\n", "$seq\n";            
    }                   
}

#-------------------------------------------------------------------------------
# Purpose   : Writes sequences to different files (.lib format), depending on 
#             their extension (int, LTR, 3end, orf etc)
# Usage     : classify_length($sequence-name, $sequence, @filehandle-list)
#             sequence name should contain '>'
# Arguments : two strings and a list, 
# Returns   : nothing, prints
# Globals   : none
#******************
sub classify_type {
    my $TE = $_[0];
    my $description = $_[1];
    my $seq = $_[2];
    my $int = '-int';

    if ($TE =~ m/_int\b|$int|_I\b|-I_/) { 
        my $FILEHANDLE = $_[3]; 
        print $FILEHANDLE ">$description\n", "$seq\n";               
    }
    elsif ($TE =~ m/_LTR\b|-LTR_/) {
        my $FILEHANDLE = $_[4]; 
        print $FILEHANDLE ">$description\n", "$seq\n";            
    }                    
    elsif ($TE =~ m/_3end\b/) {
        my $FILEHANDLE = $_[5]; 
        print $FILEHANDLE ">$description\n", "$seq\n";            
    }                    
    elsif ($TE =~ m/_5end\b/) {
        my $FILEHANDLE = $_[6]; 
        print $FILEHANDLE ">$description\n", "$seq\n";            
    }                    
    elsif ($TE =~ m/_orf2\b/) {
        my $FILEHANDLE = $_[7]; 
        print $FILEHANDLE ">$description\n", "$seq\n";            
    }                    
    else {
        my $FILEHANDLE = $_[8]; 
        print $FILEHANDLE "$description\n", "$seq\n";                
    }                     
}


#-------------------------------------------------------------------------------
# Purpose   : Builds a feature-vector for each sequence. The features are  
#             oligomer frequencies.
# Usage     : process_sequence($string, $number, *filehandle) 
# Arguments : sequence (string), prefix (number), output (filehandle) 
# Returns   : Nothing, prints.
# Globals   : Uses the global variables of the build_tetramer_hash (or pentamer) 
#             subroutines
#*******************
sub process_sequence {
    my $sequence = uc($_[0]);
    my $prefix = $_[1]; # The the SVM "classifier"
    my $OUT = $_[2];    # Filehandle

    # The vector of SVM features. Initially every element is zero
    my @features; 
    foreach my $x (0..$oligos) {
        $features[$x] = 0;
    }
    my $iterations = length($sequence)-$mer-1;
    # Counts the number of occurences of each x-mer in the sequence
    foreach my $x (0..$iterations) {
        if ( exists( $x_mers{substr($sequence, $x, $mer)} ) ) {
            $features[$x_mers{substr($sequence, $x, $mer)}]++;
            $features[$oligos+1]++;
        }
    }
    # Calculates and prints the frequencies for each oligomer
    print $OUT $prefix, "\t";
    foreach my $x (0..$oligos) {
        $features[$x] /= $features[$oligos+1]; 
        print $OUT "", ($x+1), ':', $features[$x], "\t";
    }
    print $OUT "\n";    
}


#-------------------------------------------------------------------------------
# Purpose   : Reverse complements a linear nucleotide sequence
# Usage     : reverse complement($string) 
# Arguments : nucleotide sequence 
# Returns   : reverse complemented sequence(string)
# Globals   : none
#*********************
sub reverse_complement {
    my $seq = uc($_[0]);
    $seq = reverse($seq);   
    $seq =~ tr/ATCG/TAGC/;       
    $seq =~ tr/WSYRKMBDHVN/N/;
    return $seq;
}


#-------------------------------------------------------------------------------
# Purpose   : Builds a hash of trimers
# Usage     : build_trimer_hash() 
# Arguments : none 
# Returns   : three globals
# Globals   : $oligos, $mer, %x_mers
#********************
sub build_trimer_hash {
    $oligos = 63;  # 64-1
    $mer = 3;
    %x_mers = ();
    my @nucleotides = qw / A C G T /;
    my $index = 0;
    foreach my $a (@nucleotides) {
        foreach my $b (@nucleotides) {
            foreach my $c (@nucleotides) {
                my $trimer = $a . $b . $c;
                $x_mers{$trimer} = $index++;
            }
        }
    }
    return %x_mers; 
}


#-------------------------------------------------------------------------------
# Purpose   : Builds a hash of tetramers
# Usage     : build_tetramer_hash() 
# Arguments : none 
# Returns   : three globals
# Globals   : $oligos, $mer, %x_mers
#**********************
sub build_tetramer_hash {
    $oligos = 255;  # 256-1
    $mer = 4;
    %x_mers = ();
    my @nucleotides = qw / A C G T /;
    my $index = 0;
    foreach my $a (@nucleotides) {
        foreach my $b (@nucleotides) {
            foreach my $c (@nucleotides) {
                foreach my $d (@nucleotides) {
                    my $tetramer = $a . $b . $c . $d;
                    $x_mers{$tetramer} = $index++;
                }
            }
        }
    }
    return %x_mers;
}


#-------------------------------------------------------------------------------
# Purpose   : Builds a hash of pentamers
# Usage     : build_pentamer_hash() 
# Arguments : none 
# Returns   : three globals 
# Globals   : $oligos, $mer, %x_mers
#**********************
sub build_pentamer_hash {
    $oligos = 1023;  # 1024-1
    $mer = 5;
    %x_mers = ();
    my @nucleotides = qw / A C G T /;
    my $index = 0;
    foreach my $a (@nucleotides) {
        foreach my $b (@nucleotides) {
            foreach my $c (@nucleotides) {
                foreach my $d (@nucleotides) {
                    foreach my $e (@nucleotides) {
                        my $pentamer = $a . $b . $c . $d . $e;
                        $x_mers{$pentamer} = $index++;
                    }
                }
            }
        }
    }
    return %x_mers;
}


#-------------------------------------------------------------------------------
# Purpose   : Standardises the vector file, so that the different oligomer 
#             frequencies will get te same weight              
# Usage     : &standardise_repeat_vectors($file, $analysis type);
# Arguments : $string, $string
# Returns   : nothing, prints   
# Globals   : none
#******************
sub standardise_vectors_train {
    my ($file, $analysis_type, $path_to_svm, $lower_boundary) = @_; 
    print STDERR "  - Scaling training vectors ... ";
    my $cmd = "$path_to_svm/svm-scale -l $lower_boundary -u 1 -s "
            . " $file.$analysis_type.range " 
            . " $file.$analysis_type.vect "
            . ">$file.$analysis_type.scaled";
    system($cmd);
    print STDERR "done.\n";
    unlink "$file.$analysis_type.vect";
}


#-------------------------------------------------------------------------------
# Purpose   : Translate the output format of SVM to the output format of LVQ     
# Usage     : &translate_from_SVM_to_LVQ(file, analysis_type dimensions,)
# Arguments : string, string, integer, 
# Returns   : nothing, prints
# Globals   : none
#******************
sub translate_from_SVM_to_LVQ {
    my ($file, $analysis_type, $dimensions, $path) = @_;
    
    if (defined($path)) {
        open (IN, "<", "$path/$file.$analysis_type.scaled") or die "can't open IN $!";
        open (OUT, ">", "$path/$file.$analysis_type.codebook") or die "can't open OUT $!";
    }
    else {
        open (IN, "<", "$file.$analysis_type.scaled") or die "can't open IN $!";
        open (OUT, ">", "$file.$analysis_type.codebook") or die "can't open OUT $!";
    }
    print OUT "$dimensions\tlvq\n";
    while (<IN>) {
        chomp;
        my @array = split/\s+1\:/, $_;
        $array[1] =~ s/\d+://g;
        
        my $label = 0;
        if ($analysis_type eq 'repeats') {
            if    ($array[0] == 1) { $label = 'DNA' ; } 
            elsif ($array[0] == 2) { $label = 'LTR' ; }
            elsif ($array[0] == 3) { $label = 'LINE'; }
            elsif ($array[0] == 4) { $label = 'SINE'; }
        }
        else {
            if    ($array[0] == 1) { $label = 'forward' ; } 
            elsif ($array[0] == 0) { $label = 'reverse' ; }        
        }
        print OUT $array[1], "\t", $label, "\n";
    } 
    close IN;
    close OUT;
    
    unlink "$file.$analysis_type.scaled";
}  


#-------------------------------------------------------------------------------
# Purpose   : Translate the output format of SVM to the output format of LVQ     
# Usage     : &translate_from_SVM_to_LVQ(file, analysis_type dimensions,)
# Arguments : string, string, integer, 
# Returns   : nothing, prints
# Globals   : none
#******************
sub translate_from_SVM_to_CSV {
    my ($file, $analysis_type, $path) = @_;
    if ( defined($path) ) {
        open (IN, "<", "$path/$file.$analysis_type.scaled") or die "can't open IN $!";
        open (CSV, ">", "$path/$file.$analysis_type.csv") or die "can't open CSV $!";
        open (LAB, ">", "$path/$file.$analysis_type.lab") or die "can't open LAB $!";
    }
    else {
        open (IN, "<", "$file.$analysis_type.scaled") or die "can't open IN  $!";
        open (CSV, ">", "$file.$analysis_type.csv")   or die "can't open CSV $!";
        open (LAB, ">", "$file.$analysis_type.lab")   or die "can't open LAB $!";
    }
    
    while (<IN>) {
        chomp;
        my @array = split/\s+1\:/, $_;
        $array[1] =~ s/\d+://g;
        $array[1] =~ s/\s+/,/g;
        print LAB $array[0], "\n";
        print CSV $array[1], "\n";
    } 
    close IN;
    close CSV;
    close LAB;
    
    unlink "$file.$analysis_type.scaled";
}  

#-------------------------------------------------------------------------------
# Purpose   : To build a vector of oligomer frequencies for every repeat           
# Usage     : build_repeat_vectors($string, $string, $lower, $upper);
# Arguments : $file, $analysis type, $int, $string  
# Returns   : nothing, prints  
# Globals   : none 
#******************
sub build_repeat_vectors {
    my ($file, $analysis_type, $label_1, $label_2) = @_;
    open(IN, "<", "$file.lib") or die "cannot open $file.lib $!";
    open(OUT, ">", "$file.$analysis_type.vect") 
                             or die "cannot open $file.$analysis_type.vect $!";
    my $outfile = *OUT;
    print STDERR "$file repeats:\n"; 
    print STDERR "  - Writing training set ... ";
    # Processes the linearised .lib file
    while (1) {
        last if eof(IN);
        my $repeat = <IN>;
        $repeat =~ s/>//;
        chomp($repeat); 
        my $sequence = <IN>;
        my @tags = split/\s+/, $repeat;
        my @prefix = split/_/, $tags[0];  
        $prefix[-1] =~ s/\?//;
        # Assigns a numeric value for the different repeat types
        if ($analysis_type eq 'forward_vs_reverse') {
            # Sequence in forward direction
            process_sequence($sequence, $label_1, $outfile);
            # Reverse complemented sequence
            $sequence = reverse_complement($sequence);
            process_sequence($sequence, $label_2, $outfile);
        }
        elsif ($analysis_type eq 'DNA_vs_Retro')  { 
            my $label = '';
            if    ($prefix[-1] eq DNA) { $label = $label_1; } 
            else                       { $label = $label_2; }
            process_sequence($sequence, $label, $outfile);
        }
        elsif ($analysis_type eq 'LTR_vs_nonLTR')  {
            my $label = '';
            if    ($prefix[-1] eq DNA) { next; } 
            elsif ($prefix[-1] eq LTR) { $label = $label_1; }
            else                       { $label = $label_2; }
            process_sequence($sequence, $label, $outfile);
        }
        elsif ($analysis_type eq 'LINE_vs_SINE')  {
            my $label = '';
            if    ($prefix[-1] eq DNA)  { next; } 
            elsif ($prefix[-1] eq LTR)  { next; }
            elsif ($prefix[-1] eq LINE) { $label = $label_1; }
            elsif ($prefix[-1] eq SINE) { $label = $label_2; }
            process_sequence($sequence, $label, $outfile);
        }
    }
    print STDERR "done.\n";
    close OUT; 
    close IN;
}


#-------------------------------------------------------------------------------
# Purpose   : Makes a directory or dies          
# Usage     : makedir($string);
# Arguments : $dirname  
# Returns   : the directory name 
# Globals   : none
#******************
sub makedir {
    my $dir = $_[0];
    if (! -d $dir) { 
        mkdir $dir or die "Cannot create $dir directory! \n";
    }
    elsif (defined $_[1])  {
        print STDERR "Directory $dir exists.\n";
        my $t = time;        
        $dir .= '_' . $t;
        mkdir $dir or die "Cannot create $dir directory! \n";
        print STDERR "Directory $dir created instead.\n";                 
    }
    return $dir;
}


1;
