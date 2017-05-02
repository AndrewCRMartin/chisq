#!/usr/bin/perl
#*************************************************************************
#
#   Program:    csvh2chi
#   File:       csvh2chi.pl
#   
#   Version:    V1.0
#   Date:       02.05.17
#   Function:   Convert a CSV file to the required format for the chisq
#               program
#   
#   Copyright:  (c) Dr. Andrew C. R. Martin 2017
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
#   Another program for converting CSV files to chisq format.
#
#*************************************************************************
#
#   Usage:
#   ======
#   Requires the CSV data to contain a title row only and labels in the
#   first columns which are combined with the labels in the title row.
#
#   This program can deal with data files for both chisq and chisq3
#
#   For example for a standard chisq test
#      Folded,Curved,Count
#      1,1,2
#      1,2,20
#      2,1,31
#      2,2,13
#   would be converted to
#      Folded-1 Curved-1 2
#      Folded-1 Curved-2 20
#      Folded-2 Curved-1 31
#      Folded-2 Curved-2 13
#
#   For a 3D chisq test, for example:
#      Folded,Curved,Extended,Count
#      0,1,3,2
#      0,2,1,20
#      1,1,2,31
#      1,2,3,13
#   would be converted to
#      Folded-0 Curved-1 Extended-3 2
#      Folded-0 Curved-2 Extended-1 20
#      Folded-1 Curved-1 Extended-2 31
#      Folded-1 Curved-2 Extended-3 13
#
#*************************************************************************
#
#   Revision History:
#   =================
#   V1.0  02.05.17 Original  By: ACRM
#
#*************************************************************************

use strict;

$_ = <>;
chomp;
my @headers = split(/,/);

while(<>)
{
    chomp;
    my @data = split(/,/);
    my $nData = scalar(@headers);
    for(my $i=0; $i<$nData-1; $i++)
    {
        print "$headers[$i]-$data[$i] ";
    }
    print "$data[$nData-1]\n";
}
