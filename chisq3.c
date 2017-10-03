/*************************************************************************

   Program:    chisq3
   File:       chisq3.c
   
   Version:    V1.0
   Date:       28.05.17
   Function:   Do general chi squared analysis
   
   Copyright:  (c) Dr. Andrew C. R. Martin, 2017
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure and Modelling,
               University College London,
   EMail:      andrew@bioinf.org.uk
               
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
   V1.0  28.05.17 Original based on ChiSq

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Defines
*/
#define MAXITEM 100
#define MAXBUFF 160
#define SMALL   (0.1e-20)

/************************************************************************/
/* Globals
*/
BOOL gDisplay      = FALSE,
     gGotExpecteds = FALSE;
char gItemList1[MAXITEM][MAXBUFF],
     gItemList2[MAXITEM][MAXBUFF],
     gItemList3[MAXITEM][MAXBUFF];
int  gNItem1 = 0, 
     gNItem2 = 0, 
     gNItem3 = 0;
REAL gExpecteds[MAXITEM][MAXITEM][MAXITEM];

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void ZeroMatrix(int matrix[MAXITEM][MAXITEM][MAXITEM]);
BOOL ReadData(FILE *in, int matrix[MAXITEM][MAXITEM][MAXITEM]);
REAL CalcChiSq(int matrix[MAXITEM][MAXITEM][MAXITEM], int *NDoF);
int CalcNDoF(void);
void Usage(void);
void PrintMatrix(int matrix[MAXITEM][MAXITEM][MAXITEM]);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile);
int CountObservations(int matrix[MAXITEM][MAXITEM][MAXITEM]);
void CalculateExpecteds(int matrix[MAXITEM][MAXITEM][MAXITEM], int NObs);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for chi squared calculation

   21.06.94 Original    By: ACRM
*/
int main(int argc, char **argv)
{
   FILE *in = stdin,
        *out = stdout;
   REAL chisq;
   int  dof;
   static int  matrix[MAXITEM][MAXITEM][MAXITEM];
   char InFile[160], OutFile[160];

   if(!ParseCmdLine(argc, argv, InFile, OutFile))
   {
      Usage();
   }
   else
   {
      if(blOpenStdFiles(InFile, OutFile, &in, &out))
      {
         ZeroMatrix(matrix);
   
         if(ReadData(in, matrix))
         {
            if(gDisplay)
               PrintMatrix(matrix);
            
            chisq = CalcChiSq(matrix, &dof);
            printf("ChiSq = %f with %d degrees of freedom\n", chisq, dof);
         }
      }
      else
      {
         fprintf(stderr,"Unable to open I/O files\n");
         return(1);
      }
   }
   return(0);
}

/************************************************************************/
/*>void ZeroMatrix(int matrix[MAXITEM][MAXITEM][MAXITEM])
   ------------------------------------------------------
   Zeroes the data matrix

   09.02.94 Original    By: ACRM
*/
void ZeroMatrix(int matrix[MAXITEM][MAXITEM][MAXITEM])
{
   int i,j,k;

   /* Zero the matrix                                                   */
   for(i=0; i<MAXITEM; i++)
   {
      for(j=0; j<MAXITEM; j++)
      {
         for(k=0; k<MAXITEM; k++)
         {
            matrix[i][j][k] = 0;
         }
      }
   }
}

/************************************************************************/
/*>BOOL ReadData(FILE *in, int matrix[MAXITEM][MAXITEM][MAXITEM])
   -----------------------------------------------------
   Read data into the matrix

   21.06.94 Original    By: ACRM
   04.03.08 Added reading of expecteds
*/
BOOL ReadData(FILE *in, int matrix[MAXITEM][MAXITEM][MAXITEM])
{
   int  count;
   char buffer[MAXBUFF];
   char item1[MAXBUFF], item2[MAXBUFF], item3[MAXBUFF];
   int  MatPos1,    MatPos2, MatPos3,
      i, j, k;
   REAL expect = (REAL)0.0;

   for(i=0; i<MAXITEM; i++)
      for(j=0; j<MAXITEM; j++)
         for(k=0; k<MAXITEM; k++)
            matrix[i][j][k] = 0;

   while(fgets(buffer,MAXBUFF,in))
   {
      TERMINATE(buffer);
      
      if(gGotExpecteds)
      {
         REAL countDouble;
         sscanf(buffer,"%s%s%s%lf%lf",item1,item2,item3,&countDouble,&expect);
         count = (int)countDouble;
      }
      else
      {
         REAL countDouble;
         sscanf(buffer,"%s %s %s %lf",item1,item2,item3,&countDouble);
         count = (int)countDouble;
      }

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

      /* Find the matrix position for the third item                   */
      MatPos3 = (-1);
      for(i=0; i<gNItem3; i++)
      {
         if(!strcmp(item3,gItemList3[i]))
         {
            MatPos3 = i;
            break;
         }
      }
      if(MatPos3 == (-1))
      {
         strcpy(gItemList3[gNItem3],item3);
         MatPos3 = gNItem3++;
         if(MatPos3 >= MAXITEM)
         {
            fprintf(stderr,"Too many items in third column\n");
            return(FALSE);
         }
      }

      /* Fill in the value in the matrix                                */
      matrix[MatPos1][MatPos2][MatPos3] = count;
      if(gGotExpecteds)
      {
         gExpecteds[MatPos1][MatPos2][MatPos3] = expect;
      }
   }

   return(TRUE);
}

/************************************************************************/
int CountObservations(int matrix[MAXITEM][MAXITEM][MAXITEM])
{
   int row, col, plane, 
       NObs = 0;

   /* Find total number of observations                                 */
   for(row=0; row<MAXITEM; row++)
   {
      for(col=0; col<MAXITEM; col++)
      {
         for(plane=0; plane<MAXITEM; plane++)
         {
            NObs += matrix[row][col][plane];
         }
      }
   }

   return(NObs);
}

/************************************************************************/
void CalculateExpecteds(int matrix[MAXITEM][MAXITEM][MAXITEM], int NObs)
{
   int RowTotals[MAXITEM],
       ColTotals[MAXITEM],
       PlaneTotals[MAXITEM],
      row, col, plane, i;

   /* Zero the arrays */
   for(i=0; i<MAXITEM; i++)
   {
      RowTotals[i] = ColTotals[i] = PlaneTotals[i] = 0;
   }

   /* Sum the totals */
   for(row=0; row<MAXITEM; row++)
   {
      for(col=0; col<MAXITEM; col++)
      {
         for(plane=0; plane<MAXITEM; plane++)
         {
            RowTotals[row]     += matrix[row][col][plane];
            ColTotals[col]     += matrix[row][col][plane];
            PlaneTotals[plane] += matrix[row][col][plane];
         }
      }
   }

   /* Now calculate the expecteds */
   for(row=0; row<MAXITEM; row++)
   {
      for(col=0; col<MAXITEM; col++)
      {
         for(plane=0; plane<MAXITEM; plane++)
         {
            if(RowTotals[row] &&
               ColTotals[col] &&
               PlaneTotals[plane])
            {
               gExpecteds[row][col][plane] = ((REAL)RowTotals[row] *
                                              (REAL)ColTotals[col] *
                                              (REAL)PlaneTotals[plane]) /
                  ((REAL)NObs * (REAL)NObs); 
            }
         }
      }
   }
}

/************************************************************************/
/*>REAL CalcChiSq(int matrix[MAXITEM][MAXITEM][MAXITEM], int *NDoF)
   -------------------------------------------------------
   Actually calculate the Chi squared value

   09.02.94 Original    By: ACRM
   16.12.94 Cast values in calculation of expected (was being done as int)
   06.08.03 Added Yates correction
   03.04.08 Added obtaining expecteds from file
*/
REAL CalcChiSq(int matrix[MAXITEM][MAXITEM][MAXITEM], int *NDoF)
{
   REAL chisq = (REAL)0.0,
        observed,
        expected;
   int  row, col, plane, NObs = 0, NCells = 0, NSmall = 0, NZero = 0;

   NObs = CountObservations(matrix);
   if(!gGotExpecteds)
   {
      CalculateExpecteds(matrix, NObs);
   }

   if(gDisplay)
      printf("\nTotal observations: %d\n\n",NObs);

   /* Calc DoFs                                                         */
   *NDoF = CalcNDoF();

   NCells = gNItem1 * gNItem2 * gNItem3;

   /* Step through positions in matrix                                  */
   for(row=0; row<gNItem1; row++)
   {
      for(col=0; col<gNItem2; col++)
      {
         for(plane=0; plane<gNItem3; plane++)
         {
            expected = (REAL)gExpecteds[row][col][plane];
            observed = (REAL)matrix[row][col][plane];
            
            if(gDisplay)
               printf("%s %s %s: Obs %5.1f Exp %5.1f\n",
                      gItemList1[row],gItemList2[col],gItemList3[plane],
                      observed,expected);
            
            /* Add to chisq value                                       */
            if(expected > SMALL)
            {
               if(expected < (REAL)5)
               {
                  NSmall++;
               }
               
               chisq += (observed - expected)*(observed - expected) / 
                  expected;
            }
            else
            {
               NZero++;
            }
         }
      }
   }

   if(NZero)
   {
      fprintf(stderr,"Warning: %d expecteds were < %g and not included\n",
              NZero, SMALL);
   }
               
   if ((NSmall / (REAL)NCells) > 0.25)
   {
      fprintf(stderr,"Warning: More than 25%% of expecteds were < 5\n");
   }

   return(chisq);
}

/************************************************************************/
/*>int CalcNDoF(void)
   --------------------------------------------------
   Calculate number of degress of freedom

   09.02.94 Original    By: ACRM
*/
int CalcNDoF(void)
{
   return((gNItem1-1) * (gNItem2-1) * (gNItem3-1));
}

/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message

   28.05.17 Original    By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"ChiSq3 V1.0 (c) 2017 Andrew C.R. Martin, UCL\n");
   fprintf(stderr,"Usage: chisq3 [-d] [-f] [-e] [in [out]]\n");
   fprintf(stderr,"       -d Display observed and expected values\n");
   fprintf(stderr,"       -f Use first dataset observeds as expecteds\n");
   fprintf(stderr,"       -e Expected values appear in 5th column\n");
   fprintf(stderr,"\nInput file has format: item1 item2 item3 NObs [Exp]\n");
   fprintf(stderr,"Max dimensions of contingency table: %d x %d x %d\n\n",
           MAXITEM, MAXITEM, MAXITEM);
}

/************************************************************************/
/*>void PrintMatrix(int matrix[MAXITEM][MAXITEM][MAXITEM])
   ----------------------------------------------
   Print out the matrix

   09.02.94 Original    By: ACRM
   15.12.94 Changed print field from 3 to 5
*/
void PrintMatrix(int matrix[MAXITEM][MAXITEM][MAXITEM])
{
   int i, j, k, jtot;
   
   for(k=0; k<gNItem3; k++)
   {
      printf("Plane %d\n", k+1);
      
      for(i=0; i<gNItem1; i++)
      {
         jtot = 0;
         for(j=0; j<gNItem2; j++)
         {
            printf("%5d ",matrix[i][j][k]);
            jtot += matrix[i][j][k];
         }
         printf(" : %d\n",jtot);
      }
      printf("\n");
   }
   
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
   Globals: int    gDisplay
   Returns: BOOL                Success?

   Parse the command line
   
   06.08.03 Original    By: ACRM
   16.06.09 Added -f
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile)
{
   int minArgs = 0,
       maxArgs = 2;

   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   
   if(argc < minArgs)
      return(FALSE);
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'd':
            gDisplay = TRUE;
            break;
         case 'e':
            gGotExpecteds = TRUE;
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are correct number of arguments left       */
         if((argc < minArgs) || (argc > maxArgs))
            return(FALSE);
         
         /* Copy the first to infile                                    */
         if(argc)
         {
            strcpy(infile, argv[0]);
            argc--;
            argv++;

            /* If there's another, copy it to outfile                   */
            if(argc)
            {
               strcpy(outfile, argv[0]);
               argc--;
               argv++;
            }
         }
         
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}

