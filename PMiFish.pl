#!/usr/bin/perl
use strict;
use warnings;

system "perl ./Scripts/1_Preprocessing.pl";
system "perl ./Scripts/1_4_Rarefaction.pl";
system "perl ./Scripts/2_1_Find_unique.pl";
system "perl ./Scripts/2_2_Denoise.pl";
system "perl ./Scripts/2_3_Separate_chimera.pl";
system "perl ./Scripts/2_4_Rarefaction.pl";
system "perl ./Scripts/3_1_Usearch_global.pl";
system "perl ./Scripts/4_1_Annotation.pl";
system "perl ./Scripts/5_1_Fasta_for_Phylogenetic_Analysis.pl";
system "perl ./Scripts/5_2_Summary_table.pl";
system "perl ./Scripts/5_3_Fasta_classified_by_family.pl";
system "perl ./Scripts/6_1_Portal.pl";
system "perl ./Scripts/6_2_Phylogenetic_trees.pl";
