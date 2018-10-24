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


void LWT( WTfloat * x, int loglen , int WT_steps, WTfloat* wt) 
{
// LWT is the linear wavelet transform also known as CDF(2,2)
// This code is based on the lifting scheme
// s - smooth samples at various scales, also storage for the accumulated detail coefficients
// d - detail at various scale

// limit the maximum number of transform stages 

int steps=WT_steps;
if (loglen<WT_steps)
    steps=loglen;

s=x;    
wt=[];
for (int j=1; j<steps; j++)
{
    LWT_stage(s, d, new_s);   // compute LWT_stage
    wt = [wt d];          // pack coefficients into the output stream
}

// also pack the last stage smooth coeffs
wt=[wt new_s];
 
}


void LWT_stage(WTfloat * y, int leny, WTfloat * d, WTfloat *  s)
{
	// s - smooth coefficients
	// d - detail coefficients

	for(int i=0; i<(leny/2); i++)
	{
		if (i==0)
			// front boundary treatment equivalent to inserting a zero
			d[0] = y[0]-( 0 + y[1]);
		else
			d[i] = y[2*i-1] - (y[2*i] + y[2*i-2])/2;
	}

	for(int i=0; i<(leny/2);i++)
	{
		// primal lifting
		if (i<(leny/2-1))
			s[i] = y[2*i] + (d[i] + d[i+1])/4;
		else
			// end boundary treatment equivalent to inserting a zero
			s[i] = y[2*i] + (d[i] + 0 )/2;
	}

	return;
}







void ILWT(x, loglen, wt, WT_steps)
{
	// s - smooth coefficients
	// d - detail coefficients
	// n is the data length, but it is not used in this implementation yet

	lenx=2^loglen;

	// limit the maximum number of transform stages 
	steps=WT_steps;
	if (loglen<WT_steps)
		steps=loglen;

	offset=lenx/2^(steps-1); // offset from the END of the data array
	wt=x;
	for (int j=0; j<steps;j++)
	{
		wt=wt((lenx-offset+1):lenx);
		ILWT_stage(wt, s);  // compute LWT_stage
		wt = [x(1:(lenx-offset)) s];             // pack coefficients into the output stream
		offset=2*offset;
	}
	return;
}

void ILWT_stage(wt, yr)
{
// s - smooth coefficients
// d - detail coefficients

len=length(wt);
d=wt(1:(len/2));
s=wt((len/2+1):len);

//// inverse transform
yr=zeros(1,length(wt));

for i=1:(len/2-1)
    yr(2*i) = s(i) - rounding((d(i)+d(i+1))/4);
end

// front boundary treatment: no division corresponds to symmetrical
// reflection over the first element, division by 2 correspond to prepending
// with zero
yr[0]= d[0] + yr[1];
// back boundary treatment: more symmetrical treatment is achieved by the
// factor of 2
yr[len-1]= s[len/2-1]-d[len/2-1]/ 2;

for (int i=0; i<(len/2-1); i++)
    yr[2*i+1] = d[i+1] + (yr[2*i]+yr[2*i+2])/2;

return; 
}




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


IFWT(aaa,DATALENGTH,WT_steps);
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

