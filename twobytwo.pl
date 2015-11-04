#!/usr/bin/perl -s
#*************************************************************************
#
#   Program:    twobytwo
#   File:       twobytwo.pl
#   
#   Version:    V1.0
#   Date:       02.03.12
#   Function:   Create a 2x2 contingency table from a larger table where
#               the cells are x,y  | ~x,y
#                             x,~y | ~x,~y
#   
#   Copyright:  (c) UCL / Dr. Andrew C. R. Martin 2012
#   Author:     Dr. Andrew C. R. Martin
#   Address:    Biomolecular Structure & Modelling Unit,
#               Department of Biochemistry & Molecular Biology,
#               University College,
#               Gower Street,
#               London.
#               WC1E 6BT.
#   EMail:      andrew@bioinf.org.uk
#               andrew.martin@ucl.ac.uk
#   Web:        http://www.bioinf.org.uk/
#               
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
#   Used the file format of the chisq and fisher programs
#
#*************************************************************************
#
#   Usage:
#   ======
#   twobytwo key1 key2 file.chisq
#
#*************************************************************************
#
#   Revision History:
#   =================
#
#*************************************************************************
if(defined($::h))
{
    UsageDie();
}

# Get keys from command line
my $req1 = shift(@ARGV);
my $req2 = shift(@ARGV);

if(($req1 eq "") || ($req2 eq ""))
{
    UsageDie();
}

# Initialize values
my %values;
my %keys1 = ();
my %keys2 = ();

# Read data file
while(<>)
{
    chomp;
    my ($key1, $key2, $value) = split();
    $keys1{$key1} = 1;
    $keys2{$key2} = 1;
    $values{$key1}{$key2} = $value;
}

# Check specified keys exist
if(!defined($keys1{$req1}))
{
    die "$req1 did not appear in the data file";
}
if(!defined($keys2{$req2}))
{
    die "$req2 did not appear in the data file";
}

# Obtain the 4 cells
my $tl = $values{$req1}{$req2};

my $tr = 0;
foreach my $key2 (keys %keys2)
{
    if($key2 ne $req2)
    {
        $tr += $values{$req1}{$key2};
    }
}

my $bl = 0;
foreach my $key1 (keys %keys1)
{
    if($key1 ne $req1)
    {
        $bl += $values{$key1}{$req2};
    }
}

my $br = 0;
foreach my $key1 (keys %keys1)
{
    if($key1 ne $req1)
    {
        foreach my $key2 (keys %keys2)
        {
            if($key2 ne $req2)
            {
                $br += $values{$key1}{$key2};
            }
        }
    }
}

# Print results
print "$req1 $req2 $tl\n";
print "$req1 !$req2 $tr\n";
print "!$req1 $req2 $bl\n";
print "!$req1 !$req2 $br\n";


sub UsageDie
{
    print <<__EOF;

twobytwo V1.0 (c) 2012, Dr. Andrew C.R. Martin, UCL

Usage: twobytwo key1 key2 [datafile]

Creates a 2x2 contingency table from a larger table where the cells are
     x,y  | ~x,y
     x,~y | ~x,~y
This is used with the datafiles for the chisq and fisher programs to
create a 2x2 matrix for testing individual cell significance.

__EOF

   exit 1;
}
