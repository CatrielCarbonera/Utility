#pragma once
#include <array>

template <size_t LOG_DATALENGTH = 6 , size_t wt_steps = LOG_DATALENGTH /* the largest scale level. Must be smaller or equal to LOG_DATALENGTH.*/>
struct FastWaveletTransform
{

    typedef double scalar_type;
    static size_t const data_length = 1<<LOG_DATALENGTH;    // Must be power of 2 for the inverse transform to work.

    // Direct transform can be run with DATALENGTH=48.
    typedef std::array<scalar_type, data_length> array_type;
    static_assert(wt_steps <= LOG_DATALENGTH, "Error - wt_steps > LOG_DATALENGTH");

    //static size_t const wt_steps;      // the largest scale level. Must be smaller or equal to LOG_DATALENGTH.
    static scalar_type NORM() { return 1; }   //  /(sqrt((WTfloat)214)); // Relative to 1/Sqrt(2) normalization of Haar functions.

    static array_type & FWT(array_type & c)
    {

        size_t tmpi, i, j, offset;

        std::array<scalar_type, data_length / 2> b;

        size_t j_data_length = data_length, j_half_data_length = data_length / 2;

        for (j = 0; j < wt_steps; j++, j_data_length = j_data_length / 2, j_half_data_length = j_half_data_length / 2)
        {
            offset = data_length - j_data_length;    /* offset = data_length*(1-1/2^j) */
            for (i = 0; i < j_half_data_length; i++)
            {
                tmpi = offset + 2 * i;
                b[i] = NORM() * (c[tmpi] + c[tmpi + 1]) / ((scalar_type)2);
                c[offset + i] = NORM() * (c[tmpi] - c[tmpi + 1]) / ((scalar_type)2);
            }
            for (i = 0; i < j_half_data_length; i++)
            {
                c[data_length - j_half_data_length + i] = b[i]; /* copy sums back into the array */
            }
        }
        return c;
    }  /* Wavelet transform */

    /* Inverse wavelet transform */

    static array_type & FIWT(array_type & c)
    {
        size_t i, j, j2;

        std::array<scalar_type, data_length> b;
        size_t const ln_max = 1<<(LOG_DATALENGTH - wt_steps);
        for (j = LOG_DATALENGTH - wt_steps + 1, j2 = ln_max; j <= LOG_DATALENGTH; j++, j2 *= 2)
        {
            for (i = 0; i < j2; i++)
            {
                b[data_length - 2 * i - 1] = (c[data_length - i - 1] - c[data_length - i - j2 - 1]) / (NORM());
                b[data_length - 2 * i - 2] = (c[data_length - i - 1] + c[data_length - i - j2 - 1]) / (NORM());
            }
            for (i = data_length - 2 * j2; i < data_length; i++)
            {
                c[i] = b[i];
            }
        };
        return(c);
    }  /* Inverse wavelet transform */

    static array_type & FilterSmallerThan(array_type & c, scalar_type epsilon)
    {
        for (size_t i = 0; i < c.size(); ++i)
            if ( std::abs(c[i]) < epsilon && c[i] < 0)
                c[i] = 0;
        return c;
    }
};
