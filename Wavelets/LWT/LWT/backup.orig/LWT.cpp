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

void LWT_stage(const WTfloat * y, int leny, WTfloat * d, WTfloat *  s)
{
	// s - smooth coefficients
	// d - detail coefficients

	d[0] = y[0]-( 0 + y[1]);	// front boundary treatment equivalent to inserting a zero
	
	for(int i=1; i<(leny/2); i++)
	{
		d[i] = y[2*i] - (y[2*i+1] + y[2*i-1])/2;
	}
	for(int i=0; i<(leny/2-1);i++)
	{
		// primal lifting
		s[i] = y[2*i+1] + (d[i] + d[i+1])/4;
	}
	
	s[leny/2-1] = y[2*(leny/2-1)+1] + (d[leny/2-1] + 0 )/2;	// end boundary treatment equivalent to inserting a zero

	return;
}

void LWT(WTfloat * x, int loglen, int wt_steps) 
{
	// LWT is the linear wavelet transform also known as CDF(2,2) or DCF(5,3)
	// This code is based on the lifting scheme
	// s - smooth samples at various scales, also storage for the accumulated detail coefficients
	// d - detail at various scale

	// limit the maximum number of transform stages 

	int steps=wt_steps;
	if (loglen<wt_steps)
		steps=loglen;

	WTfloat ts[DATALENGTH], s[DATALENGTH/2];	// longest needed wavelet smooth array
	WTfloat d[DATALENGTH/2];	// longest needed wavelet detail array
	int j2=0, j3=DATALENGTH/2;

	int lenx=DATALENGTH;
	int lens=DATALENGTH;
	int lenj=DATALENGTH/2;

	for(int jj=0; jj< lenx; jj++)	//Prepare transform data
	{
		ts[jj]=x[jj];
	}

	for (int j=0; j<steps; j++)		//computer the LWT stages	
	{
		LWT_stage(ts, lens, d, s);   // compute LWT_stage
		for(int i=j2; i<j3; i++)
		{
				x[i] = d[i-j2]; // pack coefficients into the output stream
		}
		lens=lens/2;
		for(int jj=0; jj< lens; jj++)	//Prepare transform data for the next stage
		{
			ts[jj]=s[jj];
		}
		j2=j2 + lenj; // DATALENGTH/pow(2,(j+1));			// next start of copy back operation
		lenj=lenj/2;
		j3=j2 + lenj; //DATALENGTH/pow(2,(j+2));		// next end of copy back operation
	}

	// also pack the last stage smooth coeffs
	for(int i=j2; i<lenx; i++)
	{
		x[i] = s[i-j2]; // pack coefficients into the output stream
	}

}



void ILWT_stage(WTfloat * wt, int lenwt, WTfloat * yr)
{
// yr - resulting inverse transform
// s - smooth coefficients
// d - detail coefficients

WTfloat * d = wt;			// we packed the wavelet coefficients first - d, then - s
WTfloat * s = wt + lenwt/2;

//// inverse transform
for (int i=0; i< lenwt/2; i++)
{
    yr[2*i+1] = s[i] - ((d[i]+d[i+1])/4);
}

// front boundary treatment: no division corresponds to symmetrical
// reflection over the first element, division by 2 correspond to prepending with zero

yr[0]= d[0] + yr[1];

// back boundary treatment: more symmetrical treatment is achieved by the factor of 2

yr[lenwt-1]= s[lenwt/2-1]-d[lenwt/2-1]/ 2;

for (int i=1; i<(lenwt/2); i++)
{
   yr[2*i] = d[i] + (yr[2*i-1]+yr[2*i+1])/2;
}



return; 
}



void ILWT(WTfloat * x, int loglen,  int wt_steps)
{
	// s - smooth coefficients
	// d - detail coefficients
	// n is the data length, but it is not used in this implementation yet

	int steps=wt_steps;
	if (loglen<wt_steps)
		steps=loglen;

	WTfloat s[DATALENGTH]; 
	WTfloat d[DATALENGTH/2];

	int lenx=DATALENGTH;
	int lens=DATALENGTH/2;

	int offset=lenx/(1<<steps); // (double)pow((double)2,(steps-1)); // offset from the END of the data array

	for (int j=0; j<steps; j++)
	{
//		wt=wt((lenx-offset+1):lenx);

		ILWT_stage(x+(lenx-(offset<<1)), offset<<1, s);  // compute LWT_stage

		for(int i=lenx-(offset<<1); i<DATALENGTH; i++ )
		{
			x[i]=s[i-lenx+(offset<<1)];
//		wt = [x(1:(lenx-offset)) s];             // pack coefficients into the output stream
		}
		offset=offset<<1;
	}
	return;
}


















/* Simple test follows. */

void main(){

/* Create array for the test. */

WTfloat  bbb[DATALENGTH];
WTfloat  aaa[DATALENGTH]={ 
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
//	aaa[i]=sin((double)3*i);
	bbb[i]=aaa[i];
   printf("%2d %4.4f	", i, ((WTfloat)aaa[i]));
   if(((i+1)%4) == 0){
   	printf("\n");
   }
}
printf("\n");
getch();



LWT(aaa,DATALENGTH,WT_steps);

printf("LWT transform follows\n");
for(i=0;i<DATALENGTH;i++){
   printf("%2d %4.4f	", i, ((WTfloat)aaa[i]));
   if(((i+1)%4) == 0){
   	printf("\n");
   }
}
printf("\n");
getch();


ILWT(aaa,DATALENGTH,WT_steps);

printf("inverse LWT transform follows\n");
for(i=0;i<DATALENGTH;i++){
   printf("%2d %4.4f	", i, ((WTfloat)aaa[i]));
   if(((i+1)%4) == 0){
   	printf("\n");
   }
}
printf("\n");


printf("difference follows\n");
for(i=0;i<DATALENGTH;i++){
   printf("%2d %4.4f	", i, (WTfloat)(bbb[i]-aaa[i]));
   if(((i+1)%4) == 0){
   	printf("\n");
   }
}
printf("\n");
getch();


} /* main */

