/*************************************************************************

   Program:    chitab
   File:       chitab.c
   
   Version:    V1.0
   Date:       02.03.00
   Function:   Calculate critical Chi-squared value for a given 
               significance value and number of degrees of freedom
   
   Copyright:  (c) University of Reading / Dr. Andrew C. R. Martin 2000
   Author:     Dr. Andrew C. R. Martin
   Address:    School of Animal and Microbial Sciences,
               The University of Reading,
               Whiteknights,
               P.O. Box 228,
               Reading RG6 6AJ.
               England.
   Phone:      +44 (0)118 987 5123 Extn. 7022
   Fax:        +44 (0)118 931 0180
   EMail:      a.c.r.martin@reading.ac.uk
               andrew@stagleys.demon.co.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============
   Simply prompts for a significance value and a number of degrees of 
   freedom and returns the critical Chi squared.

   This program uses the C numerics library v1.1 from WordenWare.
   Available on the web at:
      http://www.Brent.Worden.org/products/numericsc.html

   Compile as follows:
   g++ -o chitab chitab.c -I$HOME/libsrc/numerics -L$HOME/lib -lnumerics
   

**************************************************************************

   Usage:
   ======
   Just run the program - it will prompt you

**************************************************************************

   Revision History:
   =================

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include "numerics/normdist.h"

/************************************************************************/
/* Defines and macros
*/

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
void Usage(void);

/************************************************************************/
int main(int argc, char **argv)
{
   long dof;
   double p, chis;

   if(argc > 1)
   {
      if(argc == 3)
      {
         sscanf(argv[1], "%lf", &p);
         sscanf(argv[2], "%ld", &dof);
      }
      else
      {
         Usage();
         return(0);
      }
   }
   else
   {
      printf ("Enter the significance value      :  ");
      scanf("%lf", &p);
      
      printf ("Enter number of degrees of freedom: ");
      scanf("%ld", &dof);
   }
   p = (double)1.0 - p;

   chis = chisqv(p, dof);
   if(argc>1)
   {
      printf("%f\n", chis);
   }
   else
   {
      printf("Critical Chi-squared:  %f\n", chis);
   }
   
   return(0);
}


/************************************************************************/
void Usage(void)
{
   fprintf(stderr,"\nchitab V1.0 (c) 2000, Dr. Andrew C.R. Martin, \
University of Reading\n");

   fprintf(stderr,"\nUsage: chitab [significance dof]\n");

   fprintf(stderr,"\nchitab calculates the critical Chi-squared value \
for a specified \n");
   fprintf(stderr,"significance and number of degrees of freedom. These \
may be given on\n");
   fprintf(stderr,"the command line. If not, then the program will \
prompt for them.\n");
   fprintf(stderr,"\nThe significance is specified as a value < 1.0 with \
smaller values\n");
   fprintf(stderr,"indicating higher significance.\n");
   fprintf(stderr,"\nThis program uses the WordenWare C numerics library \
v1.1\n\n");
}

