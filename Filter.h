
/*
 * Filter.h
 *
 *
 */

#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "eq.h"

/* FILTER INFORMATION STRUCTURE FOR FILTER ROUTINES */

typedef struct FILTER{
    unsigned int length;       /* size of filter */
    float *history;            /* pointer to history in filter */
    float *coef;               /* pointer to coefficients of filter */
} FILTER;


#define FILTER_SECTIONS   2   /* 2 filter sections for 24 db/oct filter */

typedef struct {
    double a0, a1, a2;       /* numerator coefficients */
    double b0, b1, b2;       /* denominator coefficients */
} BIQUAD;

BIQUAD ProtoCoef[FILTER_SECTIONS];      /* Filter prototype coefficients,
                                         1 for each filter section
                                         */

void filter_constructor(FILTER * iir, filter_type_enum ftype, Eq_type_enum eq, uint8_t band_num, float freq, float gain, float bw);

//S to Z transform

void szxform(
             double *a0, double *a1, double *a2,     /* numerator coefficients */
             double *b0, double *b1, double *b2,   /* denominator coefficients */
             float fc,           /* Filter cutoff frequency */
             float fs,           /* sampling rate */
             float *gain,           /* overall gain factor */
             float *coef);         /* pointer to 4 iir coefficients */





void prewarp(double *a0, double *a1, double *a2, float fc, float fs);
void bilinear(
              double a0, double a1, double a2,    /* numerator coefficients */
              double b0, double b1, double b2,    /* denominator coefficients */
              float *gain,                                   /* overall gain factor */
              float fs,                                   /* sampling rate */
              float *coef);                         /* pointer to 4 iir coefficients */


/*
 * ----------------------------------------------------------
 *      Pre-warp the coefficients of a numerator or denominator.
 *      Note that a0 is assumed to be 1, so there is no wrapping
 *      of it.
 * ----------------------------------------------------------
 */

void prewarp(
             double *a0, double *a1, double *a2,
             float fc, float fs)
{
    double wp, pi;

    pi = 4.0 * atan(1.0);
    wp = 2.0 * fs * tan(pi * fc / fs);

    *a2 = (*a2) / (wp * wp);
    *a1 = (*a1) / wp;
}


/*
 * ----------------------------------------------------------
 * bilinear()
 *
 * Transform the numerator and denominator coefficients
 * of s-domain biquad section into corresponding
 * z-domain coefficients.
 *
 *      Store the 4 IIR coefficients in array pointed by coef
 *      in following order:
 *             beta1, beta2    (denominator)
 *             alpha1, alpha2  (numerator)
 *
 * Arguments:
 *             a0-a2   - s-domain numerator coefficients
 *             b0-b2   - s-domain denominator coefficients
 *             k               - filter gain factor. initially set to 1
 *                                and modified by each biquad section in such
 *                                a way, as to make it the coefficient by
 *                                which to multiply the overall filter gain
 *                                in order to achieve a desired overall filter
 *                                gain, specified in initial value of k.
 *             fs             - sampling rate (Hz)
 *             coef    - array of z-domain coefficients to be filled in.
 *
 * Return:
 *             On return, set coef z-domain coefficients
 * ----------------------------------------------------------
 */
void bilinear(
              double a0, double a1, double a2,    /* numerator coefficients */
              double b0, double b1, double b2,    /* denominator coefficients */
              float *gain,           /* overall gain factor */
              float fs,           /* sampling rate */
              float *coef         /* pointer to 4 iir coefficients */
              )
{
    double ad, bd;

    /* alpha (Numerator in s-domain) */
    ad = 4. * a2 * fs * fs + 2. * a1 * fs + a0;
    /* beta (Denominator in s-domain) */
    bd = 4. * b2 * fs * fs + 2. * b1* fs + b0;

    /* update gain constant for this section */
    *gain *= ad/bd;

    /* Denominator */
    *coef++ = (2. * b0 - 8. * b2 * fs * fs)
    / bd;         /* beta1 */
    *coef++ = (4. * b2 * fs * fs - 2. * b1 * fs + b0)
    / bd; /* beta2 */

    /* Nominator */
    *coef++ = (2. * a0 - 8. * a2 * fs * fs)
    / ad;         /* alpha1 */
    *coef = (4. * a2 * fs * fs - 2. * a1 * fs + a0)
    / ad;   /* alpha2 */
}


void szxform(
             double *a0, double *a1, double *a2, /* numerator coefficients */
             double *b0, double *b1, double *b2, /* denominator coefficients */
             float fc,         /* Filter cutoff frequency */
             float fs,         /* sampling rate */
             float *gain,         /* overall gain factor */
             float *coef)         /* pointer to 4 iir coefficients */
{
    /* Calculate a1 and a2 and overwrite the original values */
    prewarp(a0, a1, a2, fc, fs);
    prewarp(b0, b1, b2, fc, fs);
    bilinear(*a0, *a1, *a2, *b0, *b1, *b2, gain, fs, coef);
}


















