#!/usr/bin/perl -s
#*************************************************************************
#
#   Program:    
#   File:       
#   
#   Version:    
#   Date:       
#   Function:   
#   
#   Copyright:  (c) Dr. Andrew C. R. Martin 2011
#   Author:     Dr. Andrew C. R. Martin
#   Address:    Biomolecular Structure & Modelling Unit,
#               Department of Biochemistry & Molecular Biology,
#               University College,
#               Gower Street,
#               London.
#               WC1E 6BT.
#   EMail:      INTERNET: martin@biochem.ucl.ac.uk
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
#   cellsignificance [-nobf] [-low] file.chi 
#   -nobf Don't do Bonferroni correction
#   -low  Allow significance to be indicated for cells with low expected
#         values
#
#*************************************************************************
#
#   Revision History:
#   =================
#   V1.0  01.03.11   By: ACRM
#
#*************************************************************************
use strict;
my %data;
my %rows;
my %cols;
my %rowTotal;
my %colTotal;
my $lowOK;

my $tfile = "/tmp/cellsig_$$.chi";
my $threshold = 0.05;

# Read  the data file
my $totalObs = 0;
while(<>)
{
    chomp;
    s/^\s+//;
    if(length)
    {
        my @line = split;
        $data{$line[0]}{$line[1]} = $line[2];
        $totalObs += $line[2];

        $rows{$line[0]} = 1;
        $cols{$line[1]} = 1;
    }
}

# Count totals for the rows
foreach my $row (keys %rows)
{
    $rowTotal{$row} = 0;
    foreach my $col (keys %cols)
    {
        $rowTotal{$row} += $data{$row}{$col};
    }
}

# Count totals for the columns
foreach my $col (keys %cols)
{
    $colTotal{$col} = 0;
    foreach my $row (keys %rows)
    {
        $colTotal{$col} += $data{$row}{$col};
    }
}

# Bonferonni correction - set to 1 if -nobf used on command line
my $Bonferroni = (defined($::nobf))?1:((keys %rows) * (keys %cols));

# Run the analysis for each cell
foreach my $row (sort keys %rows)
{
    foreach my $col (sort keys %cols)
    {
        # Create a file for running the chisq program and run it
        my $cell = $data{$row}{$col};
        open(TMP, ">$tfile") || die "Can't write $tfile";
        print TMP "$row $col $cell\n";
        printf TMP "$row X$col %d\n", $rowTotal{$row} - $cell;
        printf TMP "X$row $col %d\n", $colTotal{$col} - $cell; 
        printf TMP "X$row X$col %d\n", $totalObs - $cell;
        close TMP;
        my $result = `chisq -d -y $tfile`;

        # Extract the significance results
        $result =~ s/\n/ /g;
        $result =~ /=\s+(.*)\s+with\s+(\d+)\s+degrees/;
        my $chi=$1;
        print "$row $col Chi=$chi ";
        $result =~ /$row\,\s+$col\:\s+Obs\s+(.*?)\s+Exp\s+(.*?)\s+/;
        my $obs = $1;
        my $exp = $2;
        print " [$obs $exp] ";

        printf "%s ", (($obs>$exp)?"over":"under");

        # Extract expected values and check they are OK
        if(!defined($::low))
        {
            $result =~ /.*Exp\s+(\d+\.\d).*Exp\s+(\d+\.\d).*Exp\s+(\d+\.\d).*Exp\s+(\d+\.\d)/;
            if(($1 >= 5) && ($2 >= 5) && ($3 >= 5) && ($4 >= 5))
            {
                $lowOK = 1;
            }
            else
            {
                $lowOK = 0;
            }
        }


        # Run the chisig program to get the p-value
        $result = `chisig $chi 1`;
        $result =~ /the\s+(.*)\s+level/;
        my $pvalue = $1;

        # Apply Bonferroni correction
        $pvalue *= $Bonferroni;

        # Print significance
        if($pvalue<$threshold)
        {
            if(defined($::low) || $lowOK)
            {
                printf "SIGNIFICANT (p=%.5g)", $pvalue, (($Bonferroni==1)?"":" - BF corrected");
            }
        }
        print "\n";

        unlink $tfile;
    }
}
