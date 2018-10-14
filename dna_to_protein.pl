#!/usr/bin/env perl

#author   : Gou Xiangjian (862137261@qq.com)
#date     : 2018/09/24

#usage    : perl dna_to_protein.pl <dna_fasta_file>

#function : Batch translation of DNA sequences into protein sequences by using the six-frame translation method.
#explain  : X indicates any amino acid, and * indicates stop codon (you can modify it on line 16).
#notice   : Except for the beginning, character > can't appear in the sequence name !!!
#output   : A DNA sequence produces six protein sequences, + and - indicates positive and negative strand respectively.

#load the modules
use strict;
use warnings;
use constant { STOP => '*', ANY => 'X' };
use File::Basename qw/basename/;

#get all codons and corresponding code
my %codon_code = &codon_list;

#start the main program
open I_DNA, '<', $ARGV[0], or die "can't open file $ARGV[0]:$!";
my $ba_name = basename $ARGV[0].'.pro';
open O_PRO, '>', $ba_name, or die "can't  get file $ba_name:$!";
$| = 1;
$/ = '>';
while(<I_DNA>){
    #step1 : get the id and seq of dna
    chomp;
    next unless $_;
    my ($id, $seq) = /(.*)\n([\d\D]*)\Z/;
    $seq =~ s/\s//g;
    #step2 : get the dna sequence of the different frame
    my $six_frame = &produce_six_frame($id, $seq);
    #step3 : get the protein sequence of the corresponding dna sequence, and output to file
    foreach my $pro_id (sort keys %$six_frame){
        my $pro_seq = &trans_to_protein($six_frame->{$pro_id});
        print O_PRO ">$pro_id\n";
        while($pro_seq){
            my $out_seq = substr $pro_seq, 0, 60;
            print O_PRO "$out_seq\n";
            substr($pro_seq, 0, 60) = '';
        }
    }
}
close I_DNA;
close O_PRO;

#Some subroutines as follows :
#----------------------------------------------------------------------------------------------------

#get the dna sequence of the different frame
sub produce_six_frame{
    my ($id, $seq) = @_;
    (my $an_seq = reverse $seq) =~ tr/AGTCagtc/TCAGtcag/;
    my %six_frame;
    foreach my $frame (qw/+1 +2 +3 -1 -2 -3/){
        my $new_id  = "${id}_$frame";
        my $seq_tmp = $frame =~ /\+/ ? $seq : $an_seq;
        my ($start) = $frame =~ /(\d)/;
        my $new_seq = substr $seq_tmp, $start-1;
        $six_frame{$new_id} = $new_seq;
    }
    return \%six_frame;
}

#get the protein sequence of the corresponding dna sequence
sub trans_to_protein{
    my $seq = uc shift;
    my $num = length($seq) % 3 ? int( length($seq)/3 )+1 : length($seq)/3;
    my ($sub, $protein_seq) = (0, '');
    foreach my $i (1 .. $num){
        my $codon = $i != $num ? substr $seq, $sub, 3 : substr $seq, $sub;
        $sub += 3;
        if(length $codon == 3){
            if($codon !~ /N/){
                $protein_seq .= $codon_code{$codon};
            }
            else{
                my $n_num = $codon =~ s/N/N/g;
                if($n_num >= 2){
                    $protein_seq .= ANY;
                }
                else{
                    my @codon = ();
                    foreach my $base (qw/A T G C/){
                        my $replace = $codon =~ s/N/$base/r;
                        push @codon,$replace;
                    }
                    $protein_seq .= ($codon_code{$codon[0]} eq $codon_code{$codon[1]} and $codon_code{$codon[1]} eq $codon_code{$codon[2]} and  $codon_code{$codon[2]} eq $codon_code{$codon[3]}) ? $codon_code{$codon[0]} : ANY;
                }
            }
        }
        elsif(length $codon == 2){
            if($codon =~ /N/){
                $protein_seq .= ANY;
            }
            else{
                my @codon = ();
                push @codon, $codon.$_ foreach qw/A T G C/;
                $protein_seq .= ($codon_code{$codon[0]} eq $codon_code{$codon[1]} and $codon_code{$codon[1]} eq $codon_code{$codon[2]} and  $codon_code{$codon[2]} eq $codon_code{$codon[3]}) ? $codon_code{$codon[0]} : ANY;
            }
        }
        else{
            $protein_seq .= ANY;
        }
    }
    return $protein_seq;
}

#get all codons and corresponding code
sub codon_list{
    my %temp = (
        GCA => 'A', GCC => 'A', GCG => 'A', GCT => 'A',                             # Ala   4
        TGC => 'C', TGT => 'C',                                                     # Cys   2
        GAC => 'D', GAT => 'D',                                                     # Asp   2
        GAA => 'E', GAG => 'E',                                                     # Glu   2
        TTT => 'F', TTC => 'F',                                                     # Phe   2
        GGA => 'G', GGC => 'G', GGG => 'G', GGT => 'G',                             # Gly   4
        CAC => 'H', CAT => 'H',                                                     # His   2
        ATA => 'I', ATC => 'I', ATT => 'I',                                         # Ile   3
        AAA => 'K', AAG => 'K',                                                     # Lys   2
        TTA => 'L', TTG => 'L', CTA => 'L', CTC => 'L', CTG => 'L', CTT => 'L',     # Leu   6
        ATG => 'M',                                                                 # Met   1
        AAC => 'N', AAT => 'N',                                                     # Asn   2
        CCA => 'P', CCC => 'P', CCG => 'P', CCT => 'P',                             # Pro   4
        CAA => 'Q', CAG => 'Q',                                                     # Gln   2
        AGA => 'R', AGG => 'R', CGA => 'R', CGC => 'R', CGG => 'R', CGT => 'R',     # Arg   6
        AGC => 'S', AGT => 'S', TCA => 'S', TCC => 'S', TCG => 'S', TCT => 'S',     # Ser   6
        ACA => 'T', ACC => 'T', ACG => 'T', ACT => 'T',                             # Thr   4
        GTA => 'V', GTC => 'V', GTG => 'V', GTT => 'V',                             # Val   4
        TGG => 'W',                                                                 # Trp   1
        TAC => 'Y', TAT => 'Y',                                                     # Tyr   2
        TAA => STOP,TAG => STOP,TGA => STOP,                                        # Stop  3
    );
    return %temp;
}
