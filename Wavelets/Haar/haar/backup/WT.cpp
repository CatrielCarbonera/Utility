/*
Here is simple code for direct and inverse wavelet transforms.
If NORM = Sqrt(2), then you get orthonormal Haar transform.
For NORM = 1 there is a factor of 1/2 already included in the transformation,
that takes care of character overflow.
The program is set up to work with charachter data types, but
character types are commented out for test. They will produce rounding errors
if uncommented and int floats are removed.
I tested the program with DATALENGTH =  32, 48, 64, 128 and with appropriate 
changes to LOG_DATALENGTH and WT_steps definitions.
To test the program, array with elements equal their indexes is transformed
back and forth, and the resulting arrays are printed on the screen.
*/

#include <stdio.h>
#include <conio.h>
#include <math.h>

#define DATALENGTH 64    // Must be power of 2 for the inverse transform to work.
						 // Direct transform can be run with DATALENGTH=48.

#define LOG_DATALENGTH 6 // Must be equal to log(2,DATALENGTH),  logarithm of datalength
#define WT_steps 6       // the largest scale level. Must be smaller or equal to LOG_DATALENGTH.

#define WTfloat double

WTfloat NORM = 1;   //  /(sqrt((WTfloat)214)); // Relative to 1/Sqrt(2) normalization of Haar functions.

/* Wavelet transform */

WTfloat * FWT(WTfloat * c, int data_length, int wt_steps){

int tmpi,i,j,offset;

WTfloat b[DATALENGTH/2];

int j64 = DATALENGTH, j32 = DATALENGTH/2;

 for(j=0; j<wt_steps; j++, j64=j64/2, j32=j32/2){
 		offset = data_length-j64;    /* offset = data_length*(1-1/2^j) */
     	for(i=0; i < j32; i++){
      	tmpi = offset + 2*i;
         b[i] = (WTfloat) NORM * (c[tmpi] + c[tmpi+1])/((WTfloat)2);
       	c[offset+i] = (WTfloat) NORM * (c[tmpi] - c[tmpi+1])/((WTfloat)2);
     	}
     	for(i=0; i<j32; i++){
      	c[data_length-j32+i] = b[i]; /* copy sums back into the array */
      }
  	}
	return(c);
}  /* Wavelet transform */


/* Inverse wavelet transform */

//char
WTfloat * FIWT( WTfloat * c, int data_length, int wt_steps){

int i,j,j2;

WTfloat b[DATALENGTH];

for(j= LOG_DATALENGTH-wt_steps+1, j2=pow((double) 2,(LOG_DATALENGTH-wt_steps)); j<=LOG_DATALENGTH; j++, j2*=2){
      for(i=0; i < j2; i++){
        b[data_length-2*i-1] = (c[data_length-i-1] - c[data_length-i-j2-1])/((WTfloat)NORM);
        b[data_length-2*i-2] = (c[data_length-i-1] + c[data_length-i-j2-1])/((WTfloat)NORM);
      }
      for(i=data_length-2*j2; i<data_length; i++){
			c[i] = b[i];
      }
};
 return(c);
}  /* Inverse wavelet transform */



/* Simple test follows. */

void main(){

/* Create array for the test. */

WTfloat aaa[DATALENGTH]={ 
0x9a, 0xa5, 0x1b, 0xab, 0x9e, 0xaa, 0x1c, 0x2d, 0xa7, 0x26,
0x9e, 0x9f, 0x9a, 0x1d, 0xaa, 0xa5, 0x2d, 0x1c, 0xa7, 0x26,
0x9b, 0x2b, 0x1c, 0xaa, 0x1d, 0x9f, 0xa5, 0x2d, 0xa7, 0x26
, 0x1b, 0x92, 0x9e, 0x9a, 0xa5, 0x1b, 0xab, 0x9e, 0xaa, 0x1c,
0x2d, 0xa7, 0x26, 0x9e, 0x9f, 0x9a, 0x1d, 0xaa, 0xa5, 0x2d, 
0x1c, 0xa7, 0x26,0x9b, 0x2b, 0x1c, 0xaa, 0x1d, 0x9f, 0xa5, 
0x2d, 0xa7, 0x26, 0x1b} ;

int i;



printf("We are goiing to transform the following data:\n");
for(i=0;i<DATALENGTH;i++){
	aaa[i]=sin((double)3*i);
   printf("%4d %f	", i, ((WTfloat)aaa[i]));
   if(((i+1)%8) == 0){
   	printf("\n");
   }
}
printf("\n");
getch();



FWT(aaa,DATALENGTH,WT_steps);

printf("haar transform follow\n");
for(i=0;i<DATALENGTH;i++){
   printf("%d %f	", i, ((WTfloat)aaa[i]));
   if(((i+1)%8) == 0){
   	printf("\n");
   }
}
printf("\n");
getch();


FIWT(aaa,DATALENGTH,WT_steps);
printf("inverse haar transform follow\n");
for(i=0;i<DATALENGTH;i++){
   printf("%d %f	", i, (aaa[i]));
   if(((i+1)%8) == 0){
   	printf("\n");
   }
}
printf("\n");
getch();


} /* main */

