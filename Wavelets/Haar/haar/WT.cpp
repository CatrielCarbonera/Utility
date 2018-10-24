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
#include <cmath>
#include <iostream>
#include <string>
#include "SampleInput.h"
#include "FWT.h"

#define DATALENGTH 64    // Must be power of 2 for the inverse transform to work.
// Direct transform can be run with DATALENGTH=48.

#define LOG_DATALENGTH 6 // Must be equal to log(2,DATALENGTH),  logarithm of datalength
#define WT_steps 6       // the largest scale level. Must be smaller or equal to LOG_DATALENGTH.

#define WTfloat double

WTfloat NORM = 1;   //  /(sqrt((WTfloat)214)); // Relative to 1/Sqrt(2) normalization of Haar functions.

/* Wavelet transform */
WTfloat * FWT(WTfloat * c, int data_length, int wt_steps)
{

    int tmpi, i, j, offset;

    WTfloat b[DATALENGTH / 2];

    int j64 = DATALENGTH, j32 = DATALENGTH / 2;

    for (j = 0; j < wt_steps; j++, j64 = j64 / 2, j32 = j32 / 2)
    {
        offset = data_length - j64;    /* offset = data_length*(1-1/2^j) */
        for (i = 0; i < j32; i++)
        {
            tmpi = offset + 2 * i;
            b[i] = (WTfloat)NORM * (c[tmpi] + c[tmpi + 1]) / ((WTfloat)2);
            c[offset + i] = (WTfloat)NORM * (c[tmpi] - c[tmpi + 1]) / ((WTfloat)2);
        }
        for (i = 0; i < j32; i++)
        {
            c[data_length - j32 + i] = b[i]; /* copy sums back into the array */
        }
    }
    return(c);
}  /* Wavelet transform */


/* Inverse wavelet transform */

//char
WTfloat * FIWT(WTfloat * c, int data_length, int wt_steps)
{

    int i, j, j2;

    WTfloat b[DATALENGTH];
    auto l_max = (int) std::pow(2 , (LOG_DATALENGTH - wt_steps));
    for (j = LOG_DATALENGTH - wt_steps + 1, j2 = l_max; j <= LOG_DATALENGTH; j++, j2 *= 2)
    {
        for (i = 0; i < j2; i++)
        {
            b[data_length - 2 * i - 1] = (c[data_length - i - 1] - c[data_length - i - j2 - 1]) / ((WTfloat)NORM);
            b[data_length - 2 * i - 2] = (c[data_length - i - 1] + c[data_length - i - j2 - 1]) / ((WTfloat)NORM);
        }
        for (i = data_length - 2 * j2; i < data_length; i++)
        {
            c[i] = b[i];
        }
    };
    return(c);
}  /* Inverse wavelet transform */

/* Simple test follows. */
void test_default()
{

    /* Create array for the test. */

    WTfloat aaa[DATALENGTH] = {
    0x9a, 0xa5, 0x1b, 0xab, 0x9e, 0xaa, 0x1c, 0x2d, 0xa7, 0x26,
    0x9e, 0x9f, 0x9a, 0x1d, 0xaa, 0xa5, 0x2d, 0x1c, 0xa7, 0x26,
    0x9b, 0x2b, 0x1c, 0xaa, 0x1d, 0x9f, 0xa5, 0x2d, 0xa7, 0x26
    , 0x1b, 0x92, 0x9e, 0x9a, 0xa5, 0x1b, 0xab, 0x9e, 0xaa, 0x1c,
    0x2d, 0xa7, 0x26, 0x9e, 0x9f, 0x9a, 0x1d, 0xaa, 0xa5, 0x2d,
    0x1c, 0xa7, 0x26,0x9b, 0x2b, 0x1c, 0xaa, 0x1d, 0x9f, 0xa5,
    0x2d, 0xa7, 0x26, 0x1b };

    int i;


    printf("We are goiing to transform the following data:\n");
    for (i = 0; i < DATALENGTH; i++) {
        aaa[i] = sin((double)3 * i);
        printf("%4d %f	", i, ((WTfloat)aaa[i]));
        if (((i + 1) % 8) == 0) {
            printf("\n");
        }
    }
    printf("\n");
    //getch();

    FWT(aaa, DATALENGTH, WT_steps);

    printf("haar transform follow\n");
    for (i = 0; i < DATALENGTH; i++) {
        printf("%d %f	", i, ((WTfloat)aaa[i]));
        if (((i + 1) % 8) == 0) {
            printf("\n");
        }
    }
    printf("\n");
    //getch();


    FIWT(aaa, DATALENGTH, WT_steps);
    printf("inverse haar transform follow\n");
    for (i = 0; i < DATALENGTH; i++) {
        printf("%d %f	", i, (aaa[i]));
        if (((i + 1) % 8) == 0) {
            printf("\n");
        }
    }
    printf("\n");
    //getch();


} /* main */

void test_mine() {

    /* Create array for the test. */

    WTfloat aaa[DATALENGTH] = {
         1, 1,  1, 1,  1, 1,  2, 1,  1, 1,  1, 1,  1, 1,  10, 10,
         1, 1,  1, 1,  1, 1,  1, 1,  1, 1,  1, 1,  1, 1,  1, 1,
         1, 1,  1, 1,  1, 1,  1, 1,  1, 1,  1, 1,  1, 1,  1, 1,
        -1, -1,  -1, 1,  1, 1, -1, 1,  1, 1,  1, -1,  1, 1, -1,  1,
    };
	WTfloat bbb[DATALENGTH];
	std::array<double, DATALENGTH> ccc;
    int i;
    for (i = 0; i < DATALENGTH; ++i)
        ccc[i] = bbb[i] = aaa[i];


    printf("We are goiing to transform the following data:\n");
    for (i = 0; i < DATALENGTH; i++) {
        printf("%4d %f	", i, ((WTfloat)aaa[i]));
        if (((i + 1) % 8) == 0)
        {
            printf("\n");
        }
    }
    printf("\n");
    _getch();

    FWT(aaa, DATALENGTH, WT_steps);
	FastWaveletTransform<>::FWT(ccc);
	for (i = 0; i < DATALENGTH; ++i)
	{
		std::cout << ccc[i] << " - " << aaa[i] << " = " << ccc[i] - aaa[i] << std::endl;
	}
	printf("\n");
	_getch();
	printf("haar transform follow\n");
    for (i = 0; i < DATALENGTH; i++) {
        printf("%d %f	", i, ((WTfloat)aaa[i]));
        if (((i + 1) % 8) == 0) {
            printf("\n");
        }
    }
    //printf("Remove small scales transform follow\n");
    //for (i = DATALENGTH / 4; i < DATALENGTH; ++i) {
    //    if (std::abs(aaa[i]) < 2 /*&& aaa[i] < 0*/)
    //        aaa[i] = 0;
    //    printf("%d %f	", i, ((WTfloat)aaa[i]));
    //    if (((i + 1) % 8) == 0) {
    //        printf("\n");
    //    }
    //}

    printf("\n");
    _getch();


    FIWT(aaa, DATALENGTH, WT_steps);
    printf("inverse haar transform follow\n");
    for (i = 0; i < DATALENGTH; i++) {
        printf("%d %f	", i, (aaa[i]));
        if (((i + 1) % 8) == 0) {
            printf("\n");
        }
    }
    printf("\n");
    _getch();

    printf("inverse haar transform minus original follow\n");
    for (i = 0; i < DATALENGTH; i++) {
        printf("%d %f	", i, (aaa[i] - bbb[i]));
        if (((i + 1) % 8) == 0) {
            printf("\n");
        }
    }
    printf("\n");
    _getch();
} /* main */

template <typename CONTAINER>
void print_container(std::string const & message, CONTAINER const & a)
{
    std::cout << message << std::endl;
    for (size_t i = 0; i < a.size(); i++) {
        std::cout << "["<<i <<"] " << a[i] << " ";
        if (((i + 1) % 8) == 0) {
            std::cout << std::endl;
        }
    }
}

template <typename CONTAINER>
double norm_inf(CONTAINER const & a, CONTAINER const & b)
{
    double max_diff = 0;
    for (size_t i = 0; i < a.size(); ++i)
    {
        auto diff = std::abs(b[i] - a[i]);
        if (diff > max_diff)
            max_diff = diff;
    }
    return max_diff;
}

template <typename CONTAINER>
void count_small_variations(CONTAINER const & signal, double epsilon = 2)
{
    int num_negative = 0;
    for (size_t i = 0; i < signal.size() - 1; ++i)
        num_negative += (signal[i] > signal[i + 1]) && (signal[i] - signal[i + 1]) < 2 ? 0 : 1;
    std::cout << "Number of small variations = " << num_negative << std::endl;
}

void test_points()
{
    double arc_length = 0;
    auto prev_point = points[0];
    for (auto point : points)
    {
        arc_length += length(point - prev_point);
        prev_point = point;
    }
    std::cout << "arc length = " << arc_length << std::endl;
    size_t const log_datalength = 9;
    typedef FastWaveletTransform<log_datalength, 3> fwt_type;

    prev_point = points[0];
    double prev_arc_length = 0;
    auto prev_point_arc_length = prev_arc_length;
    double delta = arc_length / fwt_type::data_length;
    fwt_type::array_type signal;
    auto convert = [](double z) {return std::round(10000 * z); };
    signal.front() = convert(points[0].z);
    signal.back() = convert(points[num_points - 1].z);
    size_t i = 0;
    arc_length = 0;
    auto delta_next = std::nextafter(delta, 2.0);
    for (auto point : points)
    {
        arc_length += length(point - prev_point);
        while (i + 1 < signal.size() && arc_length - prev_arc_length > delta_next)
        {
            auto t = (prev_arc_length + delta - prev_point_arc_length) / (arc_length - prev_point_arc_length);
            
            auto pt = ((1. - t) * prev_point + t * point);
            auto dt = length(pt - prev_point);
            prev_arc_length += delta;

            signal[++i] = convert(pt.z);
        }
        prev_point = point;
        prev_point_arc_length = arc_length;
    }
    auto original_signal = signal;
    count_small_variations(signal, 2);
    std::cout << "hit any key to continue";
    _getch();
    std::cout << std::endl;
    fwt_type::FWT(signal);
	std::cout << "hit any key to continue";
	print_container("Haar transform:", signal);
    fwt_type::FilterSmallerThan(signal, 2);
    fwt_type::FIWT(signal);
    print_container("Inverse Haar transform:", signal);

    std::cout << "Max diff between original and new signal = " << norm_inf(signal, original_signal) << std::endl;

    std::cout << "hit any key to continue";
    _getch();
    std::cout << std::endl;

    count_small_variations(signal, 2);
}

int main(int, char *[])
{
	test_mine();

    std::cout << "Hit any key to finish ..." << std::endl;;
    _getch();
    return 0;
} 