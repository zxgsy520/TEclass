#!/usr/bin/env perl

# A script to create ICM models for the identification of ORFs of 
# unknown repetitive elements. 

#First version finished 27th June 2009.

use warnings;
use Getopt::Std;
use List::Util qw/ max /;
use Dumpvalue;
use lib X 
use Configuration;
use TEclass;

getopts('hO:W:'); 
our($opt_h, $opt_O, $opt_W);


#-------------------------------------------------------------------------------
#                                     MAIN
#-------------------------------------------------------------------------------

# Pints help
if ($opt_h) {
    print STDERR "\n"
       . "Usage:\n" 
       . " ./Build_ICM_models.pl [-h] -O 'outdir' -W 'workdir' \n\n"
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
my $path_to_model;
if ($opt_O) {
    $path_to_model = makedir("$opt_O/ICM_model", 1); # TEclass_build checks the path for correctness  
}
else {
    $path_to_model = makedir("$path_to_TEclass/classifiers/ICM_model", 1);
} 


#------------------------- Builds *.embl files from ref ------------------------

system ("cat $path_to_RepBase/*.ref >$workdir/RepBase.embl");
system ("cat $path_to_RepBase/appendix/*.ref >>$workdir/RepBase.embl"); 

# Extracts the ORFs from the repbase files
open (IN, "<", "$workdir/RepBase.embl") or die "cant open IN $!";
our %seq_hash = ();
while (<IN>) {
    extract_orfs();
}
close IN; 
unlink "$workdir/RepBase.embl";

# Prints out the sequences for blastclust. Blastclust doesn't like very long 
# names so they get substituted with numbers. The translate hash is used to 
# convert them back from integers.
%translate = ();
open (RB, ">", "$workdir/RepBase.tmp") or die "Cannot open RB $!";
$j = 0;
foreach my $key (sort keys %seq_hash) {
    $j++;
    print RB ">$j\n$seq_hash{$key}\n";
    $translate{$j} = $key; 
}
close RB;

# Runs Blastclust on the extraced ORFs          
system ("$path_to_blastclust/blastclust -i $workdir/RepBase.tmp -o $workdir/blastclust.txt -p F -b F -L 0.95 -S 90");  
`rm error.log`;
unlink "$workdir/RepBase.tmp";

          
# The repeat sequences that were excluded by blastclust
open (EXCL, ">", "$workdir/Excluded.lib")   or die "Cannot open EXCL $!";
open (BLC,  "<", "$workdir/blastclust.txt") or die "Cannot open BLC  $!";    
while (<BLC>) {
    my @clusters = split;
    # Stops when there are no more clusters (i.e. a cluster has one element)
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
            print EXCL ">$clusters[$i]\n$seq_hash{$translate{$clusters[$i]}}\n";            
            delete $seq_hash{$clusters[$i]};
        }
        # In case there are multiple entries with the $longest length
        elsif ($count > 0) {
            print EXCL ">$clusters[$i]\n$seq_hash{$translate{$clusters[$i]}}\n";            
            delete $seq_hash{$clusters[$i]};            
        }
        else {
            $count++;
        }            
    }     

}
close BLC;
close EXCL;
unlink "$workdir/blastclust.txt";

# Prints out the rest of the hash    
open (TMP, ">", "$workdir/RepBase.orf") or die "can't open ORF $!";        
foreach my $key (sort keys %seq_hash) {                                                                                                                       
    print TMP ">$key\n$seq_hash{$key}\n";        
}
close TMP;

print STDERR "ORF extraction done.\n";

# Building ICM models
print STDERR "Building ICM models with glimmer3 ... ";
system("$path_to_glimmer/build-icm $path_to_model/RepBase.icm < $workdir/RepBase.orf");
system("mv $workdir/RepBase.orf $path_to_model");
print STDERR "done.\n";


#-------------------------------------------------------------------------------    
#                               - Subroutines -   
#    extract_orfs     
#
#-------------------------------------------------------------------------------

###  should be rewritten completely ...
#-------------------------------------------------------------------------------
# Purpose   : To extract the sequence of the open reading frames of the TEs              
# Usage     : &extract_orfs(), 
# Arguments : no arguments
# Returns   : nothing, prints  
# Globals   : 
#******************
sub extract_orfs {
    # Identifies a new repeat (and resets the variables) 
    if ($_ =~ m/^ID/) {
        s/\(//;
        s/\)//;
        $CDS = 0; # A switch   
        @description = ();  # Name, coordinates, coordinates ...   
        push @description, (split)[1];                         
        return;
    }

    # Identifies the ORF coordinates   
    if ($_ =~ m/^FT\s/ && $_ =~ m/\sCDS\s/ && $_ =~ /\.\./) {
        chomp;
        $CDS = 1; 
        # All fragments are in one line
        if ($_ =~ m/join/ && $_ =~ m/\)/) {
            s/\)//; 
            my $fragments = (split/\(/, $_)[1];                              
            push @description, $fragments;
        }
        # Fragments are in multiple lines
        elsif ($_ =~ m/join/ && $_ !~ m/\)/){ 
            my $fragments .= (split/\(/, $_)[1]; 
            my $line = '';        
            until ($line =~ m/\)/) {
                $line = <IN>;
                chomp($line);
                $line =~ s/FT//;
                $line =~ s/\s//g;  
                $fragments .= $line;
            }                        
            $fragments =~ s/\)//; 
            push @description, $fragments; 
                      
        } 
        # There is only one exon (in one line)
        else {
            my @exon = split;
            # Some TEs aren't correctly annotated in RepBase - excludes those
            if ($exon[2] eq '0..0') {
                $CDS = 0;
                return;
            } 
            push @description, $exon[2];                
        }          
        
        return;
    }

    # Reads in the sequence 
    if ($_ =~ m/^SQ/ && $CDS == 1) {
        # Reads in the sequence
        my $sequence;
        my $line = '';
        until ($line =~ m/^\/\/$/) {
            $line = <IN>;             
            $line =~ s/[0-9]| //g;
            chomp($line);
            $sequence .= $line;                
        } 
        # Transforms to uppercase + removes the last two slashes;
        $sequence = uc( substr($sequence, 0, length($sequence)-2) );
        
        # writes the orfs to the $seq_hash
        foreach my $i (1..$#description) {
            my $orf_seq = '';
            # One orf has several exons/fragments
            if ($description[$i] =~ m/,/) {
                my @domains = split/,/, $description[$i];
                foreach $j (0..$#domains) {
                    my @coords = split/\.\./, $domains[$j];
                    #next if ($coords[0] =~ m/\D/ || $coords[1] =~ m/\D/);
                    next if ($coords[0] > length $sequence || $coords[1] > length $sequence);
                    $orf_seq .= get_orf_seq(\@coords, $sequence);                                                                              
                }
                add_orf($orf_seq, $description[0], $description[$i]);
            }             
            # One ORF has only one exon 
            else {
                my @coords = split/\.\./, $description[$i];
                #next if ($coords[0] =~ m/\D/ || $coords[1] =~ m/\D/);
                next if ($coords[0] > length $sequence || $coords[1] > length $sequence);
                my $orf_seq = get_orf_seq(\@coords, $sequence);                
                add_orf($orf_seq, $description[0], $description[$i]);
            }          
        }
        
        return;
    }    
}

#-------------------------------------------------------------------------------
# Purpose   : To extract the sequence of the open reading frames              
# Usage     : get_orf_seq($array_ref, $sequence) 
# Arguments : $array_ref, $sequence 
# Returns   : orf_sequence  
# Globals   : none
#******************
sub get_orf_seq {    
    my ($coords, $sequence) = @_; 
    my $orf_seq;               
    # The ORF is in the forward strand
    if ($coords->[1] > $coords->[0]) {
        my $start = $coords->[0] - 1;
        my $offset = $coords->[1] - $coords->[0] + 1;
        $orf_seq = substr($sequence, $start, $offset);
    }
    # The ORF is in the reverse strand (e.g. Polintons)
    elsif ($coords->[0] > $coords->[1]) {
        my $start = $coords->[1] - 1;
        my $offset = $coords->[0]-$coords->[1] + 1;
        $orf_seq = substr($sequence, $start, $offset);
        $orf_seq = &reverse_complement($orf_seq);
    }
                
    return $orf_seq;                
} 


#-------------------------------------------------------------------------------
# Purpose   : Ads the orf to the seq_hash           
# Usage     : add_orf($sequence, $name, $coords) 
# Arguments : $sequence, $name, $coords 
# Returns   :   
# Globals   : %seq_hash
#******************
sub add_orf {
    my ($orf, $TEname, $coords) = @_;
    if (length $orf > 90) { 
        my $id = $TEname . '_' . $coords;
        $id =~ s/\(//g;
        $id =~ s/\)//g;
        $seq_hash{$id} = $orf;
    }
}    
      
      
