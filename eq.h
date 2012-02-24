/*
 * eq.h
 *
 *  Created on: Dec 20, 2011
 *      Author: ParallelsWin7
 */

#ifndef EQ_H_
#define EQ_H_
#include <stdint.h>
#include <float.h>
#include "common.h"


typedef enum Eq_type_enum {
	LPF = 0,
	HPF = 1,
	Pmetric = 2
} Eq_type_enum;

typedef enum filter_type_enum{
	BESSEL = 0,
	BWORTH =1,
	L_R = 2
}filter_type_enum;

typedef struct EqBand {
	uint8_t band_num;
	Enable_enum enable;
	Eq_type_enum type;
	float freq;
	float bw;
	float gain;
} EqBand;

//Constructors
void eqband_ctor(EqBand *eq, uint8_t band_num, Enable_enum enable, Eq_type_enum type,
                 float freq, float bw, float gain);


//Getters
uint8_t eqband_get_band_num(EqBand *eq);
Enable_enum eqband_is_enabled(EqBand *eq);
Eq_type_enum eqband_get_type(EqBand *eq);
float eqband_get_freq(EqBand *eq);
float eqband_get_bw(EqBand *eq);
float eqband_get_gain(EqBand *eq);

//Setters
void eqband_set_band_num(EqBand *eq, uint8_t num);
void eqband_enable(EqBand *eq);
void eqband_disable(EqBand *eq);
void eqband_set_type(EqBand *eq, Eq_type_enum type);
void eqband_set_freq(EqBand *eq, float value);
void eqband_set_bw(EqBand *eq, float value);
void eqband_set_gain(EqBand *eq, float value);

//Inspect
void eqband_inspect(EqBand *eq);
#endif /* EQ_H_ */
