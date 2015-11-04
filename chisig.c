/*************************************************************************

   Program:    chisig
   File:       chisig.c
   
   Version:    V1.1
   Date:       04.03.08
   Function:   Calculate significance for a Chi-squared value
   
   Copyright:  (c) University of Reading / Dr. Andrew C. R. Martin 2000
   Author:     Dr. Andrew C. R. Martin
   Address:    School of Animal and Microbial Sciences,
   EMail:      andrew@bioinf.org.uk
               
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
   Simply prompts for a Chi-squared value and a number of degrees of 
   freedom and returns the significance level.

   This program uses the C numerics library v1.1 from WordenWare.
   Available on the web at:
      http://www.Brent.Worden.org/products/numericsc.html
      http://www.codebeach.com/index.asp?authorName=Brent%20Worden
      http://sourceforge.net/projects/numeric/
      http://developers.soft112.com/www-brent-worden-org.html
      http://www.free-downloads-center.com/author/wordenware.html
      
      (See ~/acrm/libsrc/numerics/)

   Compile as follows:
   g++ -o chisig chisig.c -I$HOME/include -L$HOME/lib -lnumerics
   

**************************************************************************

   Usage:
   ======
   Just run the program - it will prompt you

**************************************************************************

   Revision History:
   =================
   V1.0  02.03.00  Original   By: ACRM
   V1.1  04.03.08  Allow values on the command line

*************************************************************************/
/* Includes
*/

/************************************************************************/
/* Defines and macros
*/

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/

#include <stdio.h>
#include "numerics/normdist.h"

int main(int argc, char **argv)
{
   long dof;
   double p, chis, chisq;

   if(argc == 3)
   {
      sscanf(argv[1], "%lf", &p);
      sscanf(argv[2], "%ld", &dof);
   }
   else
   {
      printf ("Enter Chi-squared value           :  ");
      scanf("%lf", &p);
      printf ("Enter number of degrees of freedom: ");
      scanf("%ld", &dof);
   }

   chisq = chisqp(p, dof);
   chis = (double)1.0 - chisq;
   printf("Significant at the %.13g level (1-%.20g)\n", chis, chisq);
   return(0);
}
