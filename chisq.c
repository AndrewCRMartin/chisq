/*************************************************************************

   Program:    chisq
   File:       chisq.c
   
   Version:    V1.3
   Date:       16.12.94
   Function:   Do general chi squared analysis
   
   Copyright:  (c) Dr. Andrew C. R. Martin, 1994
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure and Modelling,
               University College London,
   Phone:      HOME: +44 (0)372 275775
   EMail:      JANET:    martin@uk.ac.ucl.bioc.bsm
               INTERNET: martin@bsm.bioc.ucl.ac.uk
                         amartin@scitec.adsp.sub.org
               UUCP:     ...cbmehq!cbmuk!scitec!amartin
               
**************************************************************************

   This program is not in the public domain, but it may be freely copied
   and distributed for no charge providing this header is included.
   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! The code may not be sold commercially without prior permission 
   from the author, although it may be given away free with commercial 
   products, providing it is made clear that this program is free and that 
   the source code is provided with the program.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Notes:
   ======
   
**************************************************************************

   Revision History:
   =================
   V1.0  09.02.94 Original specialised version
   V1.1  21.06.94 Original
   V1.2  07.07.94 Documented
   V1.3  16.12.94 Changed matrix printing field width. Fixed minor bug
                  in calc of expected values

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"

/************************************************************************/
/* Defines
*/
#define MAXITEM 40
#define MAXBUFF 160
#define SMALL   (0.1e-20)

/************************************************************************/
/* Globals
*/
BOOL gDisplay = FALSE;
char gItemList1[MAXITEM][MAXBUFF],
     gItemList2[MAXITEM][MAXBUFF];
int  gNItem1 = 0, gNItem2 = 0;

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void ZeroMatrix(int matrix[MAXITEM][MAXITEM]);
BOOL ReadData(FILE *in, int matrix[MAXITEM][MAXITEM]);
REAL CalcChiSq(int matrix[MAXITEM][MAXITEM], int *NDoF);
int CalcNDoF(int Tot1[MAXITEM], int Tot2[MAXITEM]);
void Usage(void);
void PrintMatrix(int matrix[MAXITEM][MAXITEM]);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for chi squared calculation

   21.06.94 Original    By: ACRM
*/
int main(int argc, char **argv)
{
   FILE *in;
   REAL chisq;
   int  dof;
   int  matrix[MAXITEM][MAXITEM];

   if(argc < 2 || argc > 3)
   {
      Usage();
      return(1);
   }

   if(argc == 3)
   {
      argv++;
      argc--;
      if(!strcmp(argv[0],"-d"))
      {
         gDisplay = TRUE;
      }
      else
      {
         Usage();
         return(1);
      }
   }
   

   if((in = fopen(argv[1],"r"))==NULL)
   {
      fprintf(stderr,"Unable to open input file: %s\n",argv[1]);
      return(1);
   }
   
   ZeroMatrix(matrix);
   
   if(ReadData(in, matrix))
   {
      if(gDisplay)
         PrintMatrix(matrix);
      
      chisq = CalcChiSq(matrix, &dof);
      printf("ChiSq = %f with %d degrees of freedom\n", chisq, dof);
   }
}

/************************************************************************/
/*>void ZeroMatrix(int matrix[MAXITEM][MAXITEM])
   ---------------------------------------------
   Zeroes the data matrix

   09.02.94 Original    By: ACRM
*/
void ZeroMatrix(int matrix[MAXITEM][MAXITEM])
{
   int i,j;

   /* Zero the matrix                                                   */
   for(i=0; i<MAXITEM; i++)
   {
      for(j=0; j<MAXITEM; j++)
      {
         matrix[i][j] = 0;
      }
   }
}

/************************************************************************/
/*>BOOL ReadData(FILE *in, int matrix[MAXITEM][MAXITEM])
   -----------------------------------------------------
   Read data into the matrix

   21.06.94 Original    By: ACRM
*/
BOOL ReadData(FILE *in, int matrix[MAXITEM][MAXITEM])
{
   int  count;
   char buffer[MAXBUFF];
   char item1[MAXBUFF], item2[MAXBUFF];
   int  MatPos1,    MatPos2,
        i, j;

   for(i=0; i<MAXITEM; i++)
      for(j=0; j<MAXITEM; j++)
         matrix[i][j] = 0;

   while(fgets(buffer,MAXBUFF,in))
   {
      sscanf(buffer,"%s %s %d",item1,item2,&count);

      /* Find the matrix position for the first item                    */
      MatPos1 = (-1);
      for(i=0; i<gNItem1; i++)
      {
         if(!strcmp(item1,gItemList1[i]))
         {
            MatPos1 = i;
            break;
         }
      }
      if(MatPos1 == (-1))
      {
         strcpy(gItemList1[gNItem1],item1);
         MatPos1 = gNItem1++;
         if(MatPos1 >= MAXITEM)
         {
            fprintf(stderr,"Too many items in first column\n");
            return(FALSE);
         }
      }
      
      /* Find the matrix position for the second item                   */
      MatPos2 = (-1);
      for(i=0; i<gNItem2; i++)
      {
         if(!strcmp(item2,gItemList2[i]))
         {
            MatPos2 = i;
            break;
         }
      }
      if(MatPos2 == (-1))
      {
         strcpy(gItemList2[gNItem2],item2);
         MatPos2 = gNItem2++;
         if(MatPos2 >= MAXITEM)
         {
            fprintf(stderr,"Too many items in second column\n");
            return(FALSE);
         }
      }

      /* Fill in the value in the matrix                                */
      matrix[MatPos1][MatPos2] = count;
   }

   return(TRUE);
}

/************************************************************************/
/*>REAL CalcChiSq(int matrix[MAXITEM][MAXITEM], int *NDoF)
   -------------------------------------------------------
   Actually calculate the Chi squared value

   09.02.94 Original    By: ACRM
   16.12.94 Cast values in calculation of expected (was being done as int)
*/
REAL CalcChiSq(int matrix[MAXITEM][MAXITEM], int *NDoF)
{
   REAL chisq = (REAL)0.0,
        observed,
        expected;
   int  i, j, NObs = 0,
        Tot1[MAXITEM], Tot2[MAXITEM];
   
   /* Find total number of observations                                 */
   for(i=0; i<MAXITEM; i++)
   {
      Tot1[i] = Tot2[i] = 0;
      for(j=0; j<MAXITEM; j++)
      {
         NObs += matrix[i][j];
      }
   }

   if(gDisplay)
      printf("\nTotal observations: %d\n\n",NObs);

   /* Find the matrix totals                                            */
   for(i=0; i<MAXITEM; i++)
      for(j=0; j<MAXITEM; j++)
         Tot1[i] += matrix[i][j];
   for(j=0; j<MAXITEM; j++)
      for(i=0; i<MAXITEM; i++)
         Tot2[j] += matrix[i][j];

   /* Step through positions in matrix                                  */
   for(i=0; i<MAXITEM; i++)
   {
      for(j=0; j<MAXITEM; j++)
      {
         if(Tot1[i] && Tot2[j])
         {
            /* Calculate expected value at this cell                    */
            expected = (REAL)Tot1[i] * (REAL)Tot2[j] / (REAL)NObs;
            observed = (REAL)matrix[i][j];

            if(gDisplay)
               printf("%s, %s: Obs %5.1f Exp %5.1f\n",
                      gItemList1[i],gItemList2[j],observed,expected);
            
            /* Add to chisq value                                       */
            if(expected > SMALL)
               chisq += (observed - expected)*(observed - expected) / 
                  expected;
         }
      }
   }

   /* Calc DoFs                                                         */
   *NDoF = CalcNDoF(Tot1,Tot2);
   
   return(chisq);
}

/************************************************************************/
/*>int CalcNDoF(int Tot1[MAXITEM], int Tot2[MAXITEM])
   --------------------------------------------------
   Calculate number of degress of freedom

   09.02.94 Original    By: ACRM
*/
int CalcNDoF(int Tot1[MAXITEM], int Tot2[MAXITEM])
{
   int i, rows, cols;
   
   for(i=0,rows=0; i<MAXITEM; i++)
      if(Tot1[i]) rows++;

   for(i=0,cols=0; i<MAXITEM; i++)
      if(Tot2[i]) cols++;

   return((rows-1) * (cols-1));
}

/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message

   21.06.94 Original    By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"ChiSq V1.3 (c) 1994 Andrew C.R. Martin, UCL\n");
   fprintf(stderr,"Usage: chisq [-d] <in>\n");
   fprintf(stderr,"       -d Display observed and expected values\n");
   fprintf(stderr,"Input file has format: item1 item2 NObs\n");
   fprintf(stderr,"Max of each item = %d\n",MAXITEM);
}

/************************************************************************/
/*>void PrintMatrix(int matrix[MAXITEM][MAXITEM])
   ----------------------------------------------
   Print out the matrix

   09.02.94 Original    By: ACRM
   15.12.94 Changed print field from 3 to 5
*/
void PrintMatrix(int matrix[MAXITEM][MAXITEM])
{
   int i, j, jtot;
   

   for(i=0; i<gNItem1; i++)
   {
      jtot = 0;
      for(j=0; j<gNItem2; j++)
      {
         printf("%5d ",matrix[i][j]);
         jtot += matrix[i][j];
      }
      printf(" : %d\n",jtot);
   }
}


