#!/usr/bin/env perl

# A script to extract linear sequences from RepeatBase, to remove redundant 
# (very similar) repeats and to create separate files for repeats with 
# different length. 


use warnings;
use Getopt::Std;
use List::Util qw/ max /;
use Dumpvalue;
use lib X 
use Configuration;
use TEclass;

getopts('hW:'); 
our($opt_h, $opt_W);


# Pints help
if ($opt_h) {
    print STDERR "\n"
       . "Usage:\n" 
       . " ./Preprocessing.pl [-h] -W 'workdir' \n\n"
       . "flags:\n"       
       . "-h : this help\n";
    exit;
}

# Sets the working directory
my $workdir;
if ($opt_W) {
    $workdir = makedir($opt_W); 
}
else {
    $workdir = makedir("workdir");
} 


#--------------------- Extracts sequences from RepBase -------------------------
#
#-------------------------------------------------------------------------------

# Builds an embl file from *.ref 
system ("cat $path_to_RepBase/*.ref >$workdir/RepBase.embl");
system ("cat $path_to_RepBase/appendix/*.ref >>$workdir/RepBase.embl"); 

our @files      = qw/  short  medium  long  xlong /;
our @filehandle = ();

# Opening the seq files 
foreach my $i (0..$#files) {
    open ($i, ">", "$workdir/$files[$i].seq") or die "Can't open $files[$i] $!";
    push @filehandle, $i;
}


# Processing the GIRInst eition of Repbase 
open (IN,  "<", "$workdir/RepBase.embl")  or die "can't open RepBase.embl $!";
print STDERR "Processing RepBase ... ";

our %processed = ();
# This is really ugly, but uses very little memory 
while (<IN>) {
    if ($_ =~ m/^ID/) {
        identify_new_repeat($_);  
    }
    elsif ($_ =~ m/^KW/) {
        RB_select_repeat_type($_); 
    }
    elsif ($_ =~ m/^SQ/)  {
        read_in_sequence();  
    }         
}
close IN; 
print STDERR "done.\n";


# Processing the Repeatmasker edition of RepBase
open (IN, "<", "$path_to_RepeatMaskerLib/RepeatMaskerLib.embl") 
                                    or die "can't open RepeatMaskerLib.embl $!";
print STDERR "Processing RepeatMaskerLib.embl ... ";
while (<IN>) {
    if ($_ =~ m/^ID/) {
        identify_new_repeat($_);
    }
    elsif ($_ =~ m/^DE/) {    
        translate_to_RepBaseID($_);
    }    
    elsif ($_ =~ m/^CC/) {
        RM_select_repeat_type($_);
    }
    elsif ($_ =~ m/^SQ/)  {        
        read_in_sequence();
    }
}
close IN; 
print STDERR "done.\n";

# closing the *.seq files
foreach my $i (0..$#files) {
    close $filehandle[$i];
}

#--------------------------- Removes very similar sequences -------------------- 
#
#-------------------------------------------------------------------------------

# The repeat sequences that were excluded by blastclust
open (EXCL, ">", "$workdir/Excluded.lib") or die "Cannot open EXCL $!";

# Runs blastcust for each .seq file
foreach my $i (0..$#files) {
    # Reads in the sequence file into a hash      
    open (SEQ, "<", "$workdir/$files[$i].seq" ) or die "cannot open SEQ $!"; 
    my %seq_hash = ();
    while (! eof(SEQ) ) {   
        my $header = <SEQ>;
        $header =~ s/>//;
        $header =~ s/\(//g;
        $header =~ s/\)//g;
        chomp($header);       
        my $sequence = <SEQ>;                         
        chomp($sequence);
        $seq_hash{$header} = $sequence;
    }
    close SEQ; 

    # Runs blastclust 
    my $cmd = "$path_to_blastclust/blastclust "
            . "-i $workdir/$files[$i].seq "
            . "-o $workdir/blastclust.$i.txt " 
            . "-p F -b F -L 0.95 -S 90 2>/dev/null";
    system ( $cmd ); 
   `rm error.log`; 
    
    # Processes the clusters
    open (BLC, "<", "$workdir/blastclust.$i.txt") or die "cannot open BLC $!";     
    while (<BLC>) {
        my @clusters = split;
        if (@clusters == 1) {
           next;
        }    
        
        # Keeps only the longest
        my @lengths;        
        foreach my $i (0..$#clusters) {
            push @lengths, length($clusters[$i]);
        }
        
        my $longest = max(@lengths); 
        my $count = 0;            
        foreach my $i (0..$#clusters) {     
            if ($lengths[$i] < $longest) {
                print EXCL ">$clusters[$i]\n$seq_hash{$clusters[$i]}\n";            
                delete $seq_hash{$clusters[$i]};
            }
            # In case there are multiple entries with the $longest length
            elsif ($count > 0) {
                print EXCL ">$clusters[$i]\n$seq_hash{$clusters[$i]}\n";            
                delete $seq_hash{$clusters[$i]};            
            }
            else {
                $count++;
            }            
        }     
    }
    close BLC;
    
    open (LIB, ">", "$workdir/$files[$i].lib" ) or die "Can't open LIB $!"; 
    foreach my $key (sort keys %seq_hash) {                                                                                                                       
        my $sequence = $seq_hash{$key};
        print LIB ">$key\n", $sequence, "\n"; 
    }
    close LIB;
}
close EXCL;

cleanup(@files);

print STDERR "Preprocessing ready.\n";


#--------------------------------- SUBROUTINES ---------------------------------
#    identify_new_repeat
#    translate_to_RepBaseID
#    RM_select_repeat_type
#    RB_select_repeat_type
#    read_in_sequence
#    add_truncated_LINEs
#    write_sequence                           
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
# Purpose   : To identify a new repeat, and set the basic parameters              
# Usage     : identify_new_repeat($line)
# Arguments : string
# Returns   : $TE, $TE_ID,  
# Globals   : $TE, $TE_ID,  
#******************
sub identify_new_repeat {  
    my $line = $_[0]; 
    $TE = 0;
    $line =~ s/\(//g;
    $line =~ s/\)//g; 
    $TE_ID = (split /\s+/, $line)[1]; 
}    


#-------------------------------------------------------------------------------
# Purpose   : In case of different RepBase and RepeatMasker IDs it uses the one 
#             from RepBase
# Usage     : translate_to_RepBaseID($_);
# Arguments : $string_
# Returns   : $TE_ID  
# Globals   : $TE_ID
#******************
sub translate_to_RepBaseID {
    my $line = $_[0]; 
    if ($line =~ m/\bRepbaseID:/) {
        $TE_ID = (split /\s+/, $line)[2];
    }        
}

#-------------------------------------------------------------------------------
# Purpose   : To identife whehter the repeat is DNA, LTR, LINE or SINE when the 
#             input is the RepeatMasker embl file           
# Usage     : RM_select_repeat_type($_);
# Arguments : $string
# Returns   : $type  
# Globals   : $TE
#******************
sub RM_select_repeat_type {
    my $line = $_[0];
    if ($line =~ m/\bType:/) {  
        $type = (split/\s+/, $line)[2];       
        if ($type =~ m/DNA|LTR|LINE|SINE/) {
            $TE = 1;
        } 
    }
}


#-------------------------------------------------------------------------------
# Purpose   : To identife whehter the repeat is DNA, LTR, LINE or SINE when the 
#             input is the RepBase embl file             
# Usage     : RM_select_repeat_type($_);
# Arguments : $string
# Returns   :  
# Globals   : $TE, $type
#******************
sub RB_select_repeat_type {
    my $line = $_[0];
    if ($_ =~ m/DNA transposon/) {     
        $type = 'DNA';  
        $TE = 1;     
    }
    elsif ($_ =~ m/\sLTR Retrotransposon/ 
        || $_ =~ m/Endogenous retrovirus/
        || $_ =~ m/Endogenous Retrovirus/) {
        $type = 'LTR';
        $TE = 1;
    } 
    elsif ($_ =~ m/SINE\;/) {
        $type = 'SINE';
        $TE = 1;
    }
    elsif ($_ =~ m/Non-LTR Retrotransposon/ && $type ne 'SINE') {
        $type = 'LINE';
        $TE = 1 
    }
}


#-------------------------------------------------------------------------------
# Purpose   : Reads in the sequence from  embl format            
# Usage     : read_in_sequence()
# Arguments : none 
# Returns   : nothing, prints  
# Globals   : none
#*******************
sub read_in_sequence {
    if ($TE == 1) {
        my $sequence = '';
        # Reads in the sequence
        my $line ='x'; # prevents complains
        until ($line =~ m/^\/\/$/) {
            $line = <IN>;             
            $line =~ s/[0-9]| //g;
            chomp($line);
            $sequence .= $line;                
        } 
        # Transforms to uppercase + removes the last two slashes;
        $sequence = uc( substr($sequence,  0, length($sequence)-2) );
        write_sequence($sequence);
    }
}  


#-------------------------------------------------------------------------------
# Purpose   : Writes the sequence into the right *.lib file            
# Usage     : write_sequence();
# Arguments : no arguments
# Returns   : nothing, prints  
# Globals   : $TE_ID, $type, $sequence
#******************
sub write_sequence {
    my $sequence = $_[0];
    
    # Makes modified IDs for the repeats
    my ($LINE_450, $LINE_900);
    if ($type eq 'LINE') {
        $LINE_450 = $TE_ID . '_450_LINE'; 
        $LINE_900 = $TE_ID . '_900_LINE';     
        $TE_ID .= '_LINE';
    }
    else { 
        $TE_ID .= "_$type";
    }             
    
    # If not yet processed, prints the sequence to the sequence file 
    if (! exists $processed{$TE_ID}) {            
        classify_length(">$TE_ID", $sequence, @filehandle);
        # In case the repeat is a LINE it makes also 'fake' truncated sequences
        if ($type eq 'LINE') {
            my $length = length $sequence;
            if ($length > 3000) {
                my $sequence_450 = substr($sequence, $length-450);
                my $sequence_900 = substr($sequence, $length-900);
                print {$filehandle[0]} ">$LINE_450\n$sequence_450\n";
                print {$filehandle[1]} ">$LINE_900\n$sequence_900\n";                
            }
        }
    }
    
    $processed{$TE_ID} = 1;
}


#-------------------------------------------------------------------------------
# Purpose   : deletes the unnecessary files            
# Usage     : cleanup(@files);
# Arguments : @files
# Returns   : nothing  
# Globals   : $workdir
#******************
sub cleanup {
    my @files = @_;
    foreach my $i (0..$#files) {
        unlink "$workdir/blastclust.$i.txt";
        unlink "$workdir/$files[$i].seq";
    }
    unlink "$workdir/RepBase.embl";
}
