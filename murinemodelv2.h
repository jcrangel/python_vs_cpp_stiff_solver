/**********************************************************************************************/ /**
 * Copyright 2018 Julio C. Rangel
 * @file	murinemodelv2.h
 *
 * @brief	This is the murine model for the system on the Dendritic Immunotherapy Improvement for an Optimal
Control Murine Model. 2017 Hindawi
 **************************************************************************************************/

#ifndef MURINEMODEL_H
#define MURINEMODEL_H

#include "util.h"
#include <cmath>
class MurineModelv2
{
private:
	std::vector<Doub> parameters;
	const int systemSize = 8; // number of states
	int controlIndex = 3;	 // The index of the state where the control happens
	const int numParameters = 28;
	Doub aH, MuH, rH, KH, aC, MuC, rC, KC, KT, MuD, rI, MuIC, MuI;
	Doub rT, aT, eT, hT, aTBeta, eTBeta, rTBeta, MuBeta, aGammaC;
	Doub MuGamma, gMl, aMlGamma, eMlGamma, MuMl;

	void loadParameters()
	{
		aH = parameters[0];
		MuH = parameters[1];
		rH = parameters[2];
		KH = parameters[3];
		aC = parameters[4];
		MuC = parameters[5];
		rC = parameters[6];
		KC = parameters[7];
		KT = parameters[8];
		MuD = parameters[9];
		rI = parameters[10];
		MuIC = parameters[11];
		MuI = parameters[12];
		rT = parameters[13];
		aT = parameters[14];
		eT = parameters[15];
		hT = parameters[16];
		aTBeta = parameters[17];
		eTBeta = parameters[18];
		rTBeta = parameters[19];
		MuBeta = parameters[20];
		aGammaC = parameters[21];
		MuGamma = parameters[22];
		gMl = parameters[23];
		aMlGamma = parameters[24];
		eMlGamma = parameters[25];
		MuMl = parameters[26];
	}

	void changeParameters(std::vector<Doub> newParameters)
	{
		parameters = newParameters;
		loadParameters();
	}
public:
	MurineModelv2(std::vector<Doub> parameters) : parameters(parameters)
	{
		loadParameters();
	}

	void operator()(const StateType &x, StateType &dxdt, Doub /*t*/)
	{
		Doub T = x[0];
		Doub H = x[1];
		Doub CTL = x[2];
		Doub Den = x[3];
		Doub IL2 = x[4];
		Doub FBeta = x[5];
		Doub FGamma = x[6];
		Doub Ml = x[7];

		dxdt[0] = rT * T * std::log(KT / T) - (aT * T * CTL * Ml / (eT + Ml)) * ((aTBeta * FBeta + eTBeta) / (FBeta + eTBeta));
		dxdt[1] = aH - MuH * H + rH * Den*(H*(1 - H / KH));
		dxdt[2] = aC - MuC * CTL + rC * IL2*(CTL*(1 - CTL / KC));
		dxdt[3] = -MuD * Den * CTL;
		dxdt[4] = rI * H * Den - MuIC * CTL * IL2 - MuI * IL2;
		dxdt[5] = rTBeta * T - MuBeta * FBeta;
		dxdt[6] = aGammaC * CTL - MuGamma * FGamma;
		dxdt[7] = gMl + (aMlGamma * FGamma) / (FGamma + eMlGamma) - MuMl * Ml;
	}
	/**
 * This is needed for the Stiff integrator
 *
*/
	void operator()(const VectorBoost &x, VectorBoost &dxdt, Doub /*t*/)
	{
		Doub T = x[0];
		Doub H = x[1];
		Doub CTL = x[2];
		Doub Den = x[3];
		Doub IL2 = x[4];
		Doub FBeta = x[5];
		Doub FGamma = x[6];
		Doub Ml = x[7];

		dxdt[0] = rT * T * std::log(KT / T) - (aT * T * CTL * Ml / (eT + Ml)) * ((aTBeta * FBeta + eTBeta) / (FBeta + eTBeta));
		dxdt[1] = aH - MuH * H + rH * Den * (H * (1 - H / KH));
		dxdt[2] = aC - MuC * CTL + rC * IL2 * (CTL * (1 - CTL / KC));
		dxdt[3] = -MuD * Den * CTL;
		dxdt[4] = rI * H * Den - MuIC * CTL * IL2 - MuI * IL2;
		dxdt[5] = rTBeta * T - MuBeta * FBeta;
		dxdt[6] = aGammaC * CTL - MuGamma * FGamma;
		dxdt[7] = gMl + (aMlGamma * FGamma) / (FGamma + eMlGamma) - MuMl * Ml;
	}

	/**
	 * The jacobian of the system
	*/
	void operator()(const VectorBoost &x, MatrixBoost &M, Doub /*t*/, VectorBoost &dfdt) {
		Doub T = x[0];
		Doub H = x[1];
		Doub CTL = x[2];
		Doub Den = x[3];
		Doub IL2 = x[4];
		Doub FBeta = x[5];
		Doub FGamma = x[6];
		Doub Ml = x[7];
		//row 1
		M(0, 0) = -rT + rT * std::log(KT / T) - (aT*CTL*(eTBeta + aTBeta * FBeta) * Ml) / ((eTBeta + FBeta)* (eT + Ml));
		M(0, 1) = 0;
		M(0, 2) = -(aT * (eTBeta + aTBeta * FBeta)* Ml * T) / ((eTBeta + FBeta) * (eT + Ml));
		M(0, 3) = 0;
		M(0, 4) = 0;
		M(0, 5) = -((aT*aTBeta*CTL*Ml*T) / ((eTBeta + FBeta) *(eT + Ml))) + (aT*CTL *(eTBeta + aTBeta * FBeta)* Ml*T) / (std::pow(eTBeta + FBeta, 2) * (eT + Ml));
		M(0, 6) = 0;
		M(0, 7) = (aT * CTL*(eTBeta + aTBeta * FBeta)* Ml * T) / ((eTBeta + FBeta)* std::pow(eT + Ml, 2)) - (aT * CTL*(eTBeta + aTBeta * FBeta)* T) / ((eTBeta + FBeta)*(eT + Ml));
		//row 2
		M(1, 0) = 0;
		M(1, 1) = -MuH - (rH * Den * H) / KH + rH * Den* (1 - H / KH);
		M(1, 2) = 0;
		M(1, 3) = rH * H*(1 - H / KH);
		M(1, 4) = 0;
		M(1, 5) = 0;
		M(1, 6) = 0;
		M(1, 7) = 0;
		//row 3
		M(2, 0) = 0;
		M(2, 1) = 0;
		M(2, 2) = -MuC - (rC * CTL * IL2) / KC + rC * (1 - CTL / KC)*IL2;
		M(2, 3) = 0;
		M(2, 4) = rC * CTL*(1 - CTL / KC);
		M(2, 5) = 0;
		M(2, 6) = 0;
		M(2, 7) = 0;
		//row 4
		M(3, 0) = 0;
		M(3, 1) = 0;
		M(3, 2) = -MuD * Den;
		M(3, 3) = -MuD * CTL;
		M(3, 4) = 0;
		M(3, 5) = 0;
		M(3, 6) = 0;
		M(3, 7) = 0;
		//row 5
		M(4, 0) = 0;
		M(4, 1) = rI * Den;
		M(4, 2) = -MuIC * IL2;
		M(4, 3) = rI * H;
		M(4, 4) = -MuI - MuIC * CTL;
		M(4, 5) = 0;
		M(4, 6) = 0;
		M(4, 7) = 0;
		//row 6
		M(5, 0) = rTBeta;
		M(5, 1) = 0;
		M(5, 2) = 0;
		M(5, 3) = 0;
		M(5, 4) = 0;
		M(5, 5) = -MuBeta;
		M(5, 6) = 0;
		M(5, 7) = 0;

		//row 7
		M(6, 0) = 0;
		M(6, 1) = 0;
		M(6, 2) = aGammaC;
		M(6, 3) = 0;
		M(6, 4) = 0;
		M(6, 5) = 0;
		M(6, 6) = -MuGamma;
		M(6, 7) = 0;

		//row 8
		M(7, 0) = 0;
		M(7, 1) = 0;
		M(7, 2) = 0;
		M(7, 3) = 0;
		M(7, 4) = 0;
		M(7, 5) = 0;
		M(7, 6) = -((aMlGamma * FGamma) / std::pow(eMlGamma + FGamma, 2)) + aMlGamma / (eMlGamma + FGamma);
		M(7, 7) = -MuMl;

		dfdt[0] = 0;
		dfdt[1] = 0;
		dfdt[2] = 0;
		dfdt[3] = 0;
		dfdt[4] = 0;
		dfdt[5] = 0;
		dfdt[6] = 0;
		dfdt[7] = 0;
	}

	int getSystemSize() {
		return systemSize;
	}
	int getControlIndex()
	{
		return controlIndex;
	}
};

#endif
