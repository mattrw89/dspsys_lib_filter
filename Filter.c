/*
 * Filter.c
 *
 *  Created on: Feb 13, 2012
 *      Author: Daren
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "filter.h"

float iir_filter(input,iir)
float input;        /* new input sample */
FILTER *iir;        /* pointer to FILTER structure */
{
    unsigned int i;
    float *hist1_ptr,*hist2_ptr,*coef_ptr;
    float output,new_hist,history1,history2;
    
    /* allocate history array if different size than last call */
    
    if(!iir->history) {
        iir->history = (float *) calloc(2*iir->length,sizeof(float));
        if(!iir->history) {
            printf("\nUnable to allocate history array in iir_filter\n");
            exit(1);
        }
    }
    
    coef_ptr = iir->coef;                /* coefficient pointer */
    
    hist1_ptr = iir->history;            /* first history */
    hist2_ptr = hist1_ptr + 1;           /* next history */
    
    /* 1st number of coefficients array is overall input scale factor,
     * or filter gain */
    output = input * (*coef_ptr++);
    
    for (i = 0 ; i < iir->length; i++)
    {
        history1 = *hist1_ptr;           /* history values */
        history2 = *hist2_ptr;
        
        output = output - history1 * (*coef_ptr++);
        new_hist = output - history2 * (*coef_ptr++);    /* poles */
        
        output = new_hist + history1 * (*coef_ptr++);
        output = output + history2 * (*coef_ptr++);      /* zeros */
        
        *hist2_ptr++ = *hist1_ptr;
        *hist1_ptr++ = new_hist;
        hist1_ptr++;
        hist2_ptr++;
    }
    
    return(output);
}

int main (void) {
 
 FILTER iir;
 enum filter_type_enum ftype =BWORTH;
 enum Eq_type_enum eq = LPF;
 uint8_t band_num = 4;         //Filter Order (1-4)
 float freq = 15000;           //FILTER CUT OFF FREQ
 float gain = -15;               //en8ter desired dB gain
       gain = pow(10,gain/20); //convert gain to return dB magnitude.

// bw must be greater than zero!!
 float bw = 1;
 
 
 filter_constructor(&iir, ftype, eq, band_num, freq, gain, bw);
 
    return 0;
    system("PAUSE");
    
}

void filter_constructor(FILTER * iir, filter_type_enum ftype, Eq_type_enum eq, uint8_t band_num, float freq, float gain, float bw)
{
    
    float    *coef;
    float   fs;     /* Sampling frequency*/
    
    fs = 96000; 
    
    unsigned nInd;
    double   a0, a1, a2, b0, b1, b2;
    float sn1[3], sd1[3];
    float sn2[3] = {1,0,0};
    float sd2[3] = {0,0,1};
    
    if(eq != Pmetric){
        switch (ftype) {
           
           case L_R: {
                
                if(band_num == 2){
                            
                          printf("nth == 2L_R\n");
                            sd1[0] = 0;
                            sd1[1] = 1;
                            sd1[2] = 1;
                            sd2[0] = 0;
                            sd2[1] = 1;
                            sd2[2] = 1;
                          if(eq == LPF){
                                sn1[0] = 1;
                                sn1[1] = 0;
                                sn1[2] = 0;
                                sn2[0] = 1;
                                sn2[1] = 0;
                                sn2[2] = 0;
                               }
                          if(eq == HPF){
                                sn1[0] = 0;
                                sn1[1] = 1;
                                sn1[2] = 0;
                                sn2[0] = 0;
                                sn2[1] = 1;
                                sn2[2] = 0;
                               }
                            }
                if(band_num == 4){
                            
                          printf("nth == 4L_R\n");
                            sd1[0] = 1;
                            sd1[1] = sqrt(2);
                            sd1[2] = 1;
                            sd2[0] = 1;
                            sd2[1] = sqrt(2);
                            sd2[2] = 1;
                            
                            if(eq == LPF){
                                sn1[0] = 1;
                                sn1[1] = 0;
                                sn1[2] = 0;
                                sn2[0] = 1;
                                sn2[1] = 0;
                                sn2[2] = 0;
                               }
                          if(eq == HPF){
                                sn1[0] = 0;
                                sn1[1] = 1;
                                sn1[2] = 0;
                                sn2[0] = 0;
                                sn2[1] = 1;
                                sn2[2] = 0;
                               }
                            }
                            break;
                }
                
           case BWORTH: {
                     
                     if(band_num%2 == 1){
                          sd1[0] = 0;
                          sd1[1] = 1;
                          sd1[2] = 1;
                          
                          if(eq == LPF){
                                sn1[0] = 1;
                                sn1[1] = 0;
                                sn1[2] = 0;
                               }
                          if(eq == HPF){
                                sn1[0] = 0;
                                sn1[1] = 1;
                                sn1[2] = 0;
                                
                               }
                          
                          if(band_num == 3){
                                      
                          printf("nth == 3\n");
                                      sd2[0] = 1;
                                      sd2[1] = 1;
                                      sd2[2] = 1;
                          
                                      if(eq == LPF){
                                            sn1[0] = 1;
                                            sn1[1] = 0;
                                            sn1[2] = 0;
                                            sn2[0] = 1;
                                            sn2[1] = 0;
                                            sn2[2] = 0;
                                            }
                                      if(eq == HPF){
                                            sn1[0] = 0;
                                            sn1[1] = 0;
                                            sn1[2] = 1;
                                            sn2[0] = 0;
                                            sn2[1] = 0;
                                            sn2[2] = 1;
                                           }
                                      }         
                          }
                          
                     if(band_num == 2){
                                 
                          printf("Butterworth nth == 2\n");
                          sd1[0] = 1;
                          sd1[1] = sqrt(2);
                          sd1[2] = 1;
                          
                          if(eq == LPF){
                                sn1[0] = 1;
                                sn1[1] = 0;
                                sn1[2] = 0;
                               }
                          if(eq == HPF){
                               sn1[0] = 0;
                               sn1[1] = 0;
                               sn1[2] = 1;
                               }
                          }
                     if(band_num == 4) {
                          printf("band_num == 4Bworth\n");
                          sd1[0] = 1;
                          sd1[1] = 1.8477590650226;
                          sd1[2] = 1;
                          sd2[0] = 1;
                          sd2[1] = .76536686473018;
                          sd2[2] = 1;
                          if(eq == LPF) {
                               sn1[0] = 1;
                               sn1[1] = 0;
                               sn1[2] = 0;
                               sn2[0] = 1;
                               sn2[1] = 0;
                               sn2[2] = 0;
                               }
                          if (eq == HPF){
                               sn1[0] = 0;
                               sn1[1] = 0;
                               sn1[2] = 1;
                               sn2[0] = 0;
                               sn2[1] = 0;
                               sn2[2] = 1;
                               }
                          
                          }
                          break;
                          }
           case BESSEL:{
                        if(band_num%2 == 1){
                                      
                          printf("1st order Bessel\n");
                          sd1[0] = 0;
                          sd1[1] = 1;
                          sd1[2] = 1;
                          
                          if(eq == LPF){
                                sn1[0] = 1;
                                sn1[1] = 0;
                                sn1[2] = 0;
                               }
                          if(eq == HPF){
                                sn1[0] = 0;
                                sn1[1] = 1;
                                sn1[2] = 0;
                               }
                          
                          if(band_num == 3){
                          printf("3rd order Bessel-not complete\n");
                                     
                                     //Place numerator here
                          
                                      if(eq == LPF){
                                          
                                            }
                                      if(eq == HPF){
                                           
                                           }
                                      }         
                          }
                if(band_num == 2){
                                 
                          printf("2nd order Bessel\n");
                          sd1[0] = 1;
                          sd1[1] = 3;
                          sd1[2] = 3;
                          
                          if(eq == LPF){
                                sn1[0] = 3;
                                sn1[1] = 0;
                                sn1[2] = 0;
                               }
                          if(eq == HPF){
                               sn1[0] = 0;
                               sn1[1] = 0;
                               sn1[2] = 3;
                               }
                          }
                     if(band_num == 4) {
                          printf("4th order Bessel-not complete\n");
                          
                          if(eq == LPF) {
                                
                               }
                          if (eq == HPF){
                                 
                               }
                          
                          }
                break;
                }
             }
           }
    if(eq == Pmetric){
                float pi = 3.141592653897;
                
               float temp = 20*log10(gain); 
               float k = pow(10,abs(temp)/20.);
                
                sn1[2] = 1;
                sd1[0] = 1;
                sn1[0] = 1;
                sd1[2] = sn1[0];
                
                //If gain is positive then boost
                if(gain >= 1){
                      printf("parametric Boost\nGain(in db):%15.10f\nk:%15.10f\n",gain,k);
                      
                
                      
                      sn1[1] = 2*pi*k/bw;
                      sd1[1] = 2*pi/(bw);
                      
                      //Correct fc, not correct gain
                      //Its works fine with k = 0, bw = 2
                      //sn1[1] = gain/(bw);
                      //sd1[1] = 1/(bw);
                      
                      }
                      
                //If gain is less than 1 dB then cut
                if(gain < 1){
                     
                      printf("parametric Cut\nGain(in db):%15.10f\nk:%15.10f\n",gain,k);
                    
                      //sd1[1] = gain*2*pi*freq*freq/bw;
                      //sn1[1] = 2*pi*freq*freq/bw;
                      sn1[1] = 2*pi/(bw);
                      sd1[1] = k*2*pi/(bw);
                     
                      }


               // printf("\n%15.10f\n%15.10f\n",sn1[1],sd1[1]);
               // system("PAUSE");
                       
//                      
//                  float bc = -cos(2*pi*freq/fs);
//                  float ac = (1-tan(pi*bw/fs))/(1+tan(pi*bw/fs));    
//                      
//                     float beta0 = (1+ac+gain-gain*ac)*.5;
//                     float beta1 = bc+bc*ac;
//                     float beta2 = (1+ac-gain+gain*ac)*.5;
//                     float alpha1 = beta1;
//                     float alpha2=ac;
        
      //printf("Parametric Z domain coeff.");
      //printf("\n%15.10f\n%15.10f\n%15.10f\n%15.10f\n%15.10f\n", beta0,beta1,beta2,alpha1,alpha2);
            
                // system("PAUSE");
          }
    
  //printf("\n%15.10f\n%15.10f\n%15.10f\n%15.10f\n%15.10f\n%15.10f\n", sn1[0],sn1[1],sn1[2],sd1[0],sd1[1],sd1[2]);
    
    //printf("\n%15.10f\n%15.10f\n",sd1[1],sd2[1]);
    /*
     * Setup filter s-domain coefficients
     */
    /* Section 1 */
    ProtoCoef[0].a0 = sn1[0]; 
    ProtoCoef[0].a1 = sn1[1];
    ProtoCoef[0].a2 = sn1[2];
    ProtoCoef[0].b0 = sd1[2];
    ProtoCoef[0].b1 = sd1[1];
    ProtoCoef[0].b2 = sd1[0];
      
    /* Section 2 */
    ProtoCoef[1].a0 = sn2[0];
    ProtoCoef[1].a1 = sn2[1];
    ProtoCoef[1].a2 = sn2[2];
    ProtoCoef[1].b0 = sd2[2];
    ProtoCoef[1].b1 = sd2[1];
    ProtoCoef[1].b2 = sd2[0];
    
    iir->length = FILTER_SECTIONS;         /* Number of filter sections */
    
    /*
     * Allocate array of z-domain coefficients for each filter section
     * plus filter gain variable
     */
    iir->coef = (float *) calloc(4 * iir->length + 1, sizeof(float));
    if (!iir->coef)
    {
        printf("Unable to allocate coef array, exiting\n");
        exit(1);
    }
    
    coef = iir->coef + 1;     /* Skip k, or gain */
                         /* Sampling frequency (Hz) */
    
    /*
     * Compute z-domain coefficients for each biquad section
     * for new Cutoff Frequency and Resonance
     */
    for (nInd = 0; nInd < iir->length; nInd++)
    {
        a0 = ProtoCoef[nInd].a0;
        a1 = ProtoCoef[nInd].a1;
        a2 = ProtoCoef[nInd].a2;
        
        b0 = ProtoCoef[nInd].b0;
        b1 = ProtoCoef[nInd].b1 / bw;      /* Divide by resonance or Q
                                           */
        b2 = ProtoCoef[nInd].b2;
        szxform(&a0, &a1, &a2, &b0, &b1, &b2, freq, fs, &gain, coef);
        coef += 4;                       /* Point to next filter
                                          section */
    }
    
    /* Update overall filter gain in coef array */
    iir->coef[0] = gain;
    
    /* Display filter coefficients */
    for (nInd = 0; nInd < (iir->length * 4 + 1); nInd++)
        printf("C[%d] = %15.10f\n", nInd, iir->coef[nInd]);
    
    printf("\n[%15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f]\n\n", iir->coef[0], iir->coef[1], 
           iir->coef[2], iir->coef[3], iir->coef[4], iir->coef[5], iir->coef[6], 
           iir->coef[7], iir->coef[8]);
    /*
     * To process audio samples, call function iir_filter()
     * for each audio sample
     */
/*    
    iir_filter(0.5000, &iir);
    iir_filter(0.4989, &iir);
    iir_filter(0.4957, &iir);
    iir_filter(0.4904, &iir);
    iir_filter(0.4830, &iir);
    
    printf("output: %lf\n", iir_filter(0.4735, &iir));
*/
 //return 0;
 system("PAUSE");
}
