chisquared
==========

(c) 1994-2017 UCL, Andrew C.R. Martin
-------------------------------------

Simple programs to perform chi-squared tests:

- chisq - chi-squared calculation
- chisig - calculate the significance for a given chi-squared value 
and degrees of freedom 
- chitab - calculate critical chi-squared value for a given
significance and degrees of freedom 
- chisq3 - 3-way chi-squared calculation
- cellsignificance.pl - calculate significance for a single cell
- csv2chi.pl - rewrite a CSV file with table and column headers in
the required format
- csvh2chi.pl - rewrite a CSV file with column headers only in the
required format
- tabulate.pl - writes a table of chi-squared values for different
levels of significance and degrees of freedom
- twobytwo.pl - Creates a 2x2 contingency table from a larger table

Compile this with:

```
   make
```

Note it doesn't really need g++ except that the library routines use
C++ style comments and therefore need to be compiled with g++. To
get the linker to work properly, this then has to be compiled with
g++ as well!
