#include <stdlib.h>
#include <stdio.h>
#include <conio.h>
#include <math.h>


const int LOG_DATALENGTH(6);  // Must be equal to log(2,DATALENGTH),  logarithm of datalength
const int DATALENGTH(1<<LOG_DATALENGTH); // Must be power of 2 for the inverse transform to work.
const int WT_steps(LOG_DATALENGTH-1);        // the largest scale level. Must be smaller or equal to LOG_DATALENGTH.

#define WTfloat double

/* Wavelet transform */
// CDF(5,3) or LinearWaveletTransform

// direct transform

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
		WTfloat * y = x+offset;
		{
			d[0] = y[0]-y[1];	// front boundary treatment equivalent to inserting a zero

			// dual lifting
			for(int i=1; i<lens; i++)
			{
				d[i] = y[2*i] - (y[2*i+1] + y[2*i-1])/2;
			}

			// primal lifting
			for(int i=0; i<(lens-1); i++)
			{
				s[i] = y[2*i+1] + (d[i] + d[i+1])/4;
			}

			s[lens-1] = y[2*lens-1] + d[lens-1]/2;	// end boundary treatment equivalent to inserting a zero

		}

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
		int lenwt2=offset>>1;
		WTfloat * d1 = x+DATALENGTH-offset;			// reflects how we packed the wavelet coefficients: first - d, then - s
		WTfloat * s1 = d1 + lenwt2;

		for (int i=0; i<lenwt2; i++)
		{
			s[2*i+1] = s1[i] - ((d1[i]+d1[i+1])/4);
		}

		// front boundary treatment: no division corresponds to symmetrical
		// reflection over the first element, division by 2 correspond to prepending with zero

		s[0]= d1[0] + s[1];

		// back boundary treatment: more symmetrical treatment is achieved by the factor of 2

		s[2*lenwt2-1]= s1[lenwt2-1]-d1[lenwt2-1]/ 2;

		for (int i=1; i<lenwt2; i++)
		{
			s[2*i] = d1[i] + (s[2*i-1]+s[2*i+1])/2;
		}

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
	// getch();



	LWT(aaa, WT_steps);

	printf("LWT transform follows\n");
	for(i=0;i<DATALENGTH;i++){
		printf("%2d %4.4f	", i, ((WTfloat)aaa[i]));
		if(((i+1)%4) == 0){
			printf("\n");
		}
	}
	printf("\n");
	// getch();


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
	printf("\ntotal sum of abs of differences = %f\n",sdif);
	// getch();


} /* main */

