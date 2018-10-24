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

#include <stdlib.h>
#include <stdio.h>
#include <conio.h>
#include <math.h>


#define DATALENGTH 64   // Must be power of 2 for the inverse transform to work.
// Direct transform can be run with DATALENGTH=48.

#define LOG_DATALENGTH 6  // Must be equal to log(2,DATALENGTH),  logarithm of datalength
#define WT_steps 5        // the largest scale level. Must be smaller or equal to LOG_DATALENGTH.

#define WTfloat double

WTfloat NORM = 1;   //  /(sqrt((WTfloat)214)); // Relative to 1/Sqrt(2) normalization of Haar functions.

/* Wavelet transform */


// CDF(1,1) or HaarWaveletTransform

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
WTfloat * IFWT( WTfloat * c, int data_length, int wt_steps){

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



// CDF(5,3) or LinearWaveletTransform

// direct transform
void LWT_stage(const WTfloat * y, int leny2, WTfloat * d, WTfloat *  s)
{
	// s - smooth coefficients
	// d - detail coefficients

	d[0] = y[0]-y[1];	// front boundary treatment equivalent to inserting a zero

	// dual lifting
	for(int i=1; i<leny2; i++)
	{
		d[i] = y[2*i] - (y[2*i+1] + y[2*i-1])/2;
	}

	// primal lifting
	for(int i=0; i<(leny2-1); i++)
	{
		s[i] = y[2*i+1] + (d[i] + d[i+1])/4;
	}

	s[leny2-1] = y[2*leny2-1] + d[leny2-1]/2;	// end boundary treatment equivalent to inserting a zero

	return;
}

void LWT(WTfloat * x,  int wt_steps) 
{
	// LWT is the linear wavelet transform also known as CDF(2,2) or DCF(5,3)
	// This code is based on the lifting scheme
	// s - smooth samples at various scales, also storage for the accumulated detail coefficients
	// d - detail at various scale

	// limit the maximum number of transform stages 

	WTfloat s[DATALENGTH/2];	// longest needed wavelet smooth array
	WTfloat d[DATALENGTH/2];	// longest needed wavelet detail array

	int offset=0;			// current offset from the BEGINNING of the data array

	int lens=DATALENGTH/2;	// current length of smooth array

	for (int j=0; j<wt_steps; j++)		//computer the LWT stages	
	{
		LWT_stage(x+offset, lens, d, s);   // compute LWT_stage
		for(int i=offset, j2=0; i<offset+lens; i++, j2++)
		{
			x[i] = d[j2];		// pack coefficients into the output stream
			x[i+lens]=s[j2];
		}

		offset=offset + lens;	// next start of copy back operation
		lens=lens>>1;			// half the length of the smooth array 
	}
}

// inverse transform
void ILWT_stage(WTfloat * wt, int lenwt2, WTfloat * x)
{
	// yr - resulting inverse transform
	// s - smooth coefficients
	// d - detail coefficients


	WTfloat * d = wt;			// reflects how we packed the wavelet coefficients: first - d, then - s
	WTfloat * s = wt + lenwt2;

	for (int i=0; i<lenwt2; i++)
	{
		x[2*i+1] = s[i] - ((d[i]+d[i+1])/4);
	}

	// front boundary treatment: no division corresponds to symmetrical
	// reflection over the first element, division by 2 correspond to prepending with zero

	x[0]= d[0] + x[1];

	// back boundary treatment: more symmetrical treatment is achieved by the factor of 2

	x[2*lenwt2-1]= s[lenwt2-1]-d[lenwt2-1]/ 2;

	for (int i=1; i<lenwt2; i++)
	{
		x[2*i] = d[i] + (x[2*i-1]+x[2*i+1])/2;
	}

	return; 
}

void ILWT(WTfloat * x, int wt_steps)
{
	// s - smooth coefficients
	// d - detail coefficients
	// n is the data length, but it is not used in this implementation yet

	WTfloat s[DATALENGTH]; 
	WTfloat d[DATALENGTH/2];

	int lens=DATALENGTH/2;
	int offset=DATALENGTH>>(wt_steps-1);   //  offset from the END of the data array

	for (int j=0; j<wt_steps; j++)
	{
		ILWT_stage(x+DATALENGTH-offset, offset>>1, s);  // compute LWT_stage

		for(int i=DATALENGTH-offset, j2=0; i<DATALENGTH; i++, j2++)
		{
			x[i]=s[j2];   // pack coefficients into the output stream
		}
		offset=offset<<1;
	}
	return;
}





/* Simple test follows. */

void main(){

	/* Create array for the test. */

	WTfloat bbb[DATALENGTH];
	WTfloat aaa[DATALENGTH]
//	; WTfloat	xxx[DATALENGTH]
	={ 	0x9a, 0xa5, 0x1b, 0xab, 0x9e, 0xaa, 0x1c, 0x2d, 0xa7, 0x26,
		0x9e, 0x9f, 0x9a, 0x1d, 0xaa, 0xa5, 0x2d, 0x1c, 0xa7, 0x26,
		0x9b, 0x2b, 0x1c, 0xaa, 0x1d, 0x9f, 0xa5, 0x2d, 0xa7, 0x26
		, 0x1b, 0x92, 0x9e, 0x9a, 0xa5, 0x1b, 0xab, 0x9e, 0xaa, 0x1c,
		0x2d, 0xa7, 0x26, 0x9e, 0x9f, 0x9a, 0x1d, 0xaa, 0xa5, 0x2d, 
		0x1c, 0xa7, 0x26,0x9b, 0x2b, 0x1c, 0xaa, 0x1d, 0x9f, 0xa5, 
		0x2d, 0xa7, 0x26, 0x1b} ;

		int i;



		printf("We are going to transform the following data:\n");
		for(i=0;i<DATALENGTH;i++){
			aaa[i]=sin((double)i*i);  //rand();
			bbb[i]=aaa[i];
			printf("%2d %4.4f	", i, ((WTfloat)aaa[i]));
			if(((i+1)%4) == 0){
				printf("\n");
			}
		}
		printf("\n");
		getch();



		LWT(aaa, WT_steps);

		printf("LWT transform follows\n");
		for(i=0;i<DATALENGTH;i++){
			printf("%2d %4.4f	", i, ((WTfloat)aaa[i]));
			if(((i+1)%4) == 0){
				printf("\n");
			}
		}
		printf("\n");
		getch();


		ILWT(aaa, WT_steps);

		printf("inverse LWT transform follows\n");
		for(i=0;i<DATALENGTH;i++){
			printf("%2d %4.4f	", i, ((WTfloat)aaa[i]));
			if(((i+1)%4) == 0){
				printf("\n");
			}
		}
		printf("\n");


		printf("difference follows\n");
		double sdif=0;
		for(i=0;i<DATALENGTH;i++){
			printf("%2d %4.4f	", i, (WTfloat)(bbb[i]-aaa[i]));
			if(((i+1)%4) == 0){
				printf("\n");
			}
			sdif=sdif+abs((bbb[i]-aaa[i]));
		}
		printf("\ntotal sum of abs of differences = %f",sdif);
		getch();


} /* main */

