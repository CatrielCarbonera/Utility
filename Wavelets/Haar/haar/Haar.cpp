/* Wavelet coefficients calculation. */
/* Haar transform follows. */

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>

#define INTEL
#define NEWCODE

#define DATALENGTH 64

void main(){

int i, j, j64, j32, tmpi, offset, tmp2;


j64 = DATALENGTH;
j32 = DATALENGTH/2;

int b[DATALENGTH];

int a[DATALENGTH] = {
0xfd, 0xfd, 0xfe, 0xfe, 0xfe, 0xfd, 0xfe, 0xfd,
0xfe, 0x1d, 0x27, 0x29, 0x52, 0x67, 0x6a, 0x51,
0x3f, 0x38, 0xfb, 0xc6, 0xb9, 0x95, 0x82, 0x82,
0x82, 0x8a, 0x94, 0xb1, 0xed, 0xff, 0x7, 0xb,
0xe, 0x10, 0x11, 0x11, 0x10, 0x10, 0xe, 0xc,
0xc, 0xc, 0xc, 0x10, 0xf, 0x10, 0x16, 0x17};



for(j=0; j<=5; j++, j32>>=1, j64>>=1) {
  	offset = (signed char)(DATALENGTH-j64);    /* offset = DATALENGTH*(1-1/2^j) */
	for(i=0; i < j32; i++) {
      tmpi = offset+(i<<1);    /* offset + 2*i   */
/*
   	b[i]        =	(a[tmpi]>>1) + (a[tmpi+1]>>1);
      a[offset+i] = 	(a[tmpi]>>1) - (a[tmpi+1]>>1);
*/
    	tmp2 = a[tmpi+1]/2;   /* compiler must invoke ASR for >>1 to perform */
    	tmpi = a[tmpi]/2;     /* division by 2 of negative numbers correctly */
      b[i] = tmpi + tmp2;    /* sums are accumulated in b */
		tmpi -=tmp2;           /* this is the wavelet coefficient */

      a[offset+i] = tmpi; /* wavelet coeffs are stored in the original array */
	}

#ifdef INTELl
/****test****/
printf("\ncopying %d numbers starting at %d\n", j32,DATALENGTH-j32);
for(i=0;i<j32;i++){
   printf("%d %d	", i, ((int)b[i]));
   if(((i+1)%8) == 0){
   	printf("\n");
   }
}
/****end test****/
#endif

   for(i=DATALENGTH-j32; i<DATALENGTH; i++) {
   	a[i] = b[i-DATALENGTH+j32];  /* copy sums back into the original array */
   }
}


#ifdef INTEL
/****test****/
printf("Haar transform follows:\n");
for(i=0;i<DATALENGTH;i++){
   printf("%d %x	", i, ((unsigned char)a[i]));
   if(((i+1)%8) == 0){
   	printf("\n");
   }
}
printf("\n");
getch();
/****end test****/
#endif


}