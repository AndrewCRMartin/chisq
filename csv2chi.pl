#!/usr/bin/perl
#*************************************************************************
#
#   Program:    csv2chi
#   File:       csv2chi.pl
#   
#   Version:    V1.0
#   Date:       01.03.11
#   Function:   Convert a CSV file to the required format for the chisq
#               program
#   
#   Copyright:  (c) Dr. Andrew C. R. Martin 2011
#   Author:     Dr. Andrew C. R. Martin
#   Address:    Biomolecular Structure & Modelling Unit,
#               Department of Biochemistry & Molecular Biology,
#               University College,
#               Gower Street,
#               London.
#               WC1E 6BT.
#   EMail:      andrew@bioinf.org.uk
#               
#*************************************************************************
#
#   This program is not in the public domain, but it may be copied
#   according to the conditions laid out in the accompanying file
#   COPYING.DOC
#
#   The code may be modified as required, but any modifications must be
#   documented so that the person responsible can be identified. If 
#   someone else breaks this code, I don't want to be blamed for code 
#   that does not work! 
#
#   The code may not be sold commercially or included as part of a 
#   commercial product except as described in the file COPYING.DOC.
#
#*************************************************************************
#
#   Description:
#   ============
#
#*************************************************************************
#
#   Usage:
#   ======
#   Requires the CSV data to contain a title row and column with the
#   top left cell empty - e.g.
#       A   B   C   
#   a   1   5   7
#   b   2   7   2
#   c   3   5   2
#   would be represented as:
#   "","A","B","C"
#   "a","1","5","7"
#   "b","2","7","2"
#   "c","3","5","2"
#
#*************************************************************************
#
#   Revision History:
#   =================
#   V1.0  01.03.11 Original  By: ACRM
#
#*************************************************************************
use strict;
my @rowNames = ();
my @matrix;

# Read the column names
$_ = <>;
chomp;
s/\"//g;
my @colNames = split(/\,/, $_);
shift @colNames;
my $nCols = @colNames;

# Read the rows
my $nRows = 0;
while(<>)
{
   chomp; 
   s/\"//g;
   s/^\s+//;
   if(length)
   {
       my @values = split(/\,/, $_);
       push @rowNames, shift(@values);
       
       my $count = 0;
       foreach my $value (@values)
       {
           $matrix[$nRows][$count++] = $value;
       }
       $nRows++;
   }
}

for(my $r=0; $r<$nRows; $r++)
{
    for(my $c=0; $c<$nCols; $c++)
    {
        print "$rowNames[$r] $colNames[$c] $matrix[$r][$c]\n";
    }
}
