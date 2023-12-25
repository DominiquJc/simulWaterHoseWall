/*
    simulWaterHoseWall performs a thermal simulation of a thick wall, 
    several time crossed by hoses where water circulates in a closed loop.
    Copyright (C) 2023  Dominique PAUL

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
- Purpose of the program
It allows the study of the thermal behavior of this device, for various periods, and various wall and hose configurations.



- 2D temperature map
When there are hoses, the map is 2D, with the depth axis and the lateral axis. The hoses are considered a repetitive pattern.
For the convenience of expression, I called "line" the cells of the same depth.

The map elements are square cells, whose geographic boundaries are integers, and whose centers of gravity are shifted by 0.5.
The temperature is that of the center of gravity. Each cell has 4 neighbors.

The depth line "-1" is the air of the living room. Its temperature is imposed, it is a sinusoidal function.
The other cells have their temperatures calculated.

The last line is facing an insulator. It has no deeper neighbor. If a neighbor is required for programming purposes, it is itself.
A cell of the first column has a cell of the same temperature for symmetry. The first column is in the same case that last line.
Same for the last column.

The hoses are circular. Since there are only square cells, the shape of the circle will be stepped.
A cell belongs to a circle if its center of gravity is at a distance from the center of the circle less than the radius.


- 1D temperature map
When there is no hose, the map is 1D, with only the depth axis. Each cell has 2 neighbors.



- Principles for calculating temperatures
-- equation: dθ/dt = -diffusivity * (d²θ/ddeep² + d²θ/dlateral²)
this translates into: (tempPrev[i,j] - temp[i,j])/diffusivity + ((temp[i+1,j] - temp[i,j]) + (temp[i-1,j] - temp[i,j]) + (temp[i,j-1] - temp[i,j]) + (temp[i,j+1] - temp[i,j])) = 0
tempPrev being the temperature at previous time, and diffusivity = lambda / volumetric heat capacity.

-- Convergence process
To obtain values close to equilibrium, the imbalance is calculated, then the temperature is changed in proportion to the imbalance; this results in other, smaller imbalances.
Then we start again, until the corrections become insignificant.

-- Coefficient of correction
In this example:
tempDelta = ((temp[i,j-1] - temp[i,j]) + (temp[i,j+1] - temp[i,j])) + ((temp[i-1,j] - temp[i,j]) + (temp[i+1,j] - temp[i,j])) + ((tempPrev[i,j] - temp[i,j])/diffusivity)
temp[i,j] += tempDelta / (4 + 1/diffusivity)
"1 / (4 + 1/diffusivity)" is the most rapid correction coefficient.
If it is smaller, it should cause slower convergence.
If it is larger, it should cause oscillations during convergence iterations, and overall slower convergence.
If it is twice too large, there is divergence instead of convergence.

Warning: the correction coefficients must be consistent with each other.
For example, if a cell has only 3 neighbors (it is next to the insulation), it would be tempting to replace 4 by 3 in the denominator;
convergence would be faster for this cell, but would distort equilibrium temperatures.
And in the case of water, which has a different heat capacity?
The rule would be: the proportion of the flux imbalance devoted to heat inertia must be the same throughout the map.
Or again: the correction rate, corrected for the inverse of the heat capacity, must be the same throughout the map.

-- Synchronicity of input data
A correction calculation pass should use the values from the previous iteration.
However, a deviation to this rule may be made: use the just calculated temperature of the shallower neighbor instead of that of the previous iteration.
This accelerates the convergence process; this save 1/4 of the exec time.
This cause a small error on equilibrium values; but by adding some additional iterations to the 2 algorithms, this difference decreases significantly.
Here, I choose the fast algorithm.


- Approximations
-- air-earth interface
In reality, between the first line of the wall and that of the air of the room, there is a thermal resistance of contact "surface - air" (Rsi, about 0.13 m²K/W),
in addition to that of the earth cell (1/(2 lambda)) between its center of gravity and one side.
I did not retain this additional resistance to stay close to the mathematical model of the delay line.

-- earth-water interface
In reality, the border has a circular shape, which I replaced with a stair shape.
In reality, the border includes a hose of unknown thickness and unknown lambda, which I did not take into account.
In reality, the thermal resistance between a water cell and earth cell should be the sum of the water and earth cell resistance.
Instead, I took twice the thermal resistance of an earth cell. This being so, the lambda of the water is close to that of the earth.
In reality, the water cell temperatures are not identical, but I acted as if they had been instantly homogenized.

*/

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <stdbool.h>
#include <signal.h>
#include <errno.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <float.h>

/* constants */

// moderation of correction
/* if the equilibrium value of a temperature is 1 + epsilon/2,
we can avoid, during the convergence iterations, the sequence: 1, 1 + epsilon, 1, 1 + epsilon... */
#define CORR_MODERATION					(127.0/128)

// cycle number before steady state
#define CYCLENB_BEFORE_1D_STEADY			19
#ifndef NDEBUG
#define CYCLENB_BEFORE_1D_STEADY			2
#endif
#define CYCLENB_BEFORE_2D_STEADY			2
// cycle number for amplitude increasing
#define CYCLENB_AMPLITUDE_INCR				1

/* macros */

// conversion
#define FLOAT_TO_INT(x) 				((x)>=0?(int)((x)+0.5):(int)((x)-0.5))
// square
#define SQUARE(val)					(val * val)
// acceptable limit of quadratic errors sum during convergence
#define CONV_ERR_LIM(val, cellNb, denum)		(SQUARE(val * DBL_EPSILON) * cellNb / denum)
// limit number of iterations for convergence
#define CONV_ITER_NUMBER(cellNb)			(3 * cellNb)

/* debug/verbose/supervision */
#define DEBUG_SUPERV_HOSE1DEEP				hose1Deep


/* structures */

typedef struct {
	unsigned	period;
	double		omega;
	double		angle;
	double		sinus;
	double		cosinus;
} Trig;

typedef struct {
	const Trig	*trig;
	const char	*name;
	double		sinus;
	double		cosinus;
} TrigSum;

typedef struct {
	int		begIx;		// first index of the line
	int		limIx;		// end index of the line + 1
} EarthLineDesc;

typedef struct {
	int		deepIx;		// deep index
	bool		right;		// left/right
} HosePos;

typedef struct {
	int		deepIx;		// deep index
	int		lateralIx;	// left/right index (= position in the line)
} Pos;

/* globals */

static	const char *	_progname;
static	bool		_verbose = false;
static	bool		_stop = false;
static	FILE		*_resf = NULL;


/* functions */

static	const char *	_basename(const char *path) {
	const char	*str = strrchr(path, '/' );       // Find end of path.
	return (!str) ? path : ++str;
}

static	void	_lineSet(double *line, const unsigned size, const double val) {
	double		*p = line;
	double		*pLim = line + size;

	while (p < pLim)
		*(p++) = val;
}

// average on 1 line
static	double	_lineAverage(const double *line, const unsigned size) {
	const double	*p;
	const double	*pLim;
	double		sum = 0;

	sum = 0;
	for (p = line, pLim = line + size; p < pLim; p++)
		sum += *p;
	return sum / size;
}

// average of difference between 2 lines
static	double	_lineDiffAverage(const double *line, const int offset, const unsigned size) {
	const double	*p;
	const double	*pLim;
	double		sum = 0;

	sum = 0;
	for (p = line, pLim = line + size; p < pLim; p++)
		sum += (*p - p[offset]);
	return sum / size;
}

// average on 2 lines
static	double	_lineMidAverage(const double *line, const int offset, const unsigned size) {
	const double	*p;
	const double	*pLim;
	double		sum = 0;

	sum = 0;
	for (p = line, pLim = line + size; p < pLim; p++)
		sum += (*p + p[offset]);
	return sum / (2 * size);
}

// set tempCT with water temperature
static	void	_tempWaterSet(const EarthLineDesc *earthLineDescT, double *tempCT, double *lineCPLim, const unsigned lineSz, const double tempCWater) {
	const EarthLineDesc	*descP;
	double			*lineCP;		// current time line temperature pointer
	int			limMax = lineSz -1;	// = width

	for (descP = earthLineDescT, lineCP = tempCT; lineCP < lineCPLim; descP++, lineCP += lineSz) {
		if (descP->begIx > 0) {
			double	*p = lineCP;
			double	*pLim = lineCP + descP->begIx;
			do {
				*(p++) = tempCWater;
			} while (p < pLim);
		}
		if (descP->limIx < limMax) {
			double	*p = lineCP + descP->limIx;
			double	*pLim = lineCP + limMax;
			do {
				*(p++) = tempCWater;
			} while (p < pLim);
		}
	}
}

static	void	_trigSet(Trig *trig, const unsigned timeCpt) {
	double		angle;

	if (timeCpt < trig->period/2) {
		angle = trig->omega * timeCpt;
		trig->angle = angle;
		trig->sinus = sin(angle);
		trig->cosinus = cos(angle);
	} else {
		angle = trig->omega * (timeCpt - trig->period/2);
		trig->angle = M_PI + angle;
		trig->sinus = -sin(angle);
		trig->cosinus = -cos(angle);
	}
}

static	void	_trigSumReset(TrigSum *trigSum) {
	trigSum->sinus = 0.0;
	trigSum->cosinus = 0.0;
}

static	void	_trigSumInit(TrigSum *trigSum, const Trig *trig, const char *name) {
	trigSum->trig = trig;
	trigSum->name = name;
	_trigSumReset(trigSum);
}

static	void	_trigSumAdd(TrigSum *trigSum, const double val) {
	trigSum->sinus += trigSum->trig->sinus * val;
	trigSum->cosinus += trigSum->trig->cosinus * val;
}

static	void	_trigSumResult(const TrigSum *trigSum, double *module, double *phase) {
	if (!!module)
		*module = sqrt(trigSum->sinus * trigSum->sinus + trigSum->cosinus * trigSum->cosinus) * 2 / trigSum->trig->period;
	if (!!phase)
		*phase = atan2(trigSum->cosinus, trigSum->sinus);
}

static	void	_trigPrint(FILE *resf, const TrigSum *trigSum) {
	double		module;
	double		phase;

	_trigSumResult(trigSum, &module, &phase);
	fprintf(_resf, "%s: %.16e, %f°\n", trigSum->name, module, phase * 180 / M_PI);
}

static	void	_onStop(int sig_number) {
	_stop = true;
	if (_verbose) {
		printf("_onStop: signal %d received!\n", sig_number);
	}
}

static	void	_printTime(FILE *resf) {
	time_t timestamp = time( NULL );
	struct tm * pTime = localtime(& timestamp);
	char buffer[60];

	strftime(buffer, sizeof(buffer), "%F %T", pTime );
	fprintf(resf,"time=%s\n", buffer);
}

static	void	_dumpTemp(FILE *resf, const double *tempT, const unsigned sz) {
	for(unsigned ix = 0; ix < sz; ix++) {
		fprintf(_resf, "%u;%.16e\n", ix, tempT[ix]);
	}
}

static	void	_dumpTempAverage(FILE *resf, const double *tempT, double *lineCPLim, const unsigned lineSz) {
	const double		*lineP;		// time line temperature pointer
	unsigned		ix;

	for (lineP = tempT, ix = 0; lineP < lineCPLim; lineP += lineSz, ix++) {
		fprintf(_resf, "%u;%.16e\n", ix, _lineAverage(lineP, lineSz -1));
	}
}

static	int	_simul(const unsigned width, const unsigned deep, const unsigned period, const double flowDensityU, const double earthLambda, const double earthVHeatCap, const double waterVHeatCap, const unsigned hoseNb, const unsigned hoseRadius, const unsigned hose1Deep, const unsigned hoseInterval) {
	const double	earthDiffusivity = earthLambda / earthVHeatCap;	// earth diffusivity, in (cell size)^2/time unit
	const double	waterDiffusivity = earthLambda / waterVHeatCap;	// water (in earth border) diffusivity, in (cell size)^2/time unit
	const double	invEarthDiffusivity = 1 / earthDiffusivity;	// inverse of earth diffusivity
	const double	omega =  M_PI * 2 / period;			// the pulsation
	double		*temp1DT;					// cell temperature 1D (widthless) table
	Trig		timeTrig;					// trig data in the cycle related to timeCpt
	TrigSum		flowDensityTrigSum;				// for flow density analysis
	TrigSum		waterTempTrigSum;				// for water temperature analysis
	TrigSum		midDeepTempTrigSum;				// for middle deep temperature analysis
	TrigSum		eighthDeepTempTrigSum;				// for eighth deep temperature analysis
#ifdef DEBUG_SUPERV_HOSE1DEEP
	TrigSum		pos_hose1deep_av_TempTrigSum;			// for position hose1deep,average temperature analysis
	TrigSum		pos_hose1deep_width_m1_TempTrigSum;		// for position hose1deep,width-1 temperature analysis
	TrigSum		pos_hose1deep_mid_TempTrigSum;			// for position hose1deep,middle temperature analysis
	TrigSum		pos_hose1deep_nh_TempTrigSum;			// for position hose1deep,near hose temperature analysis
	TrigSum		pos_hose1deep_d2_av_TempTrigSum;		// for position hose1deep/2,average temperature analysis
#endif
	{
		__sighandler_t  oldhandler1, oldhandler2;

		oldhandler1 = signal(SIGTERM, _onStop);
		oldhandler2 = signal(SIGINT, _onStop);
		if (_verbose) {
			printf("signal oldhandler1=%p, oldhandler2=%p\n", oldhandler1, oldhandler2);
		}
	}
	timeTrig.period = period;
	timeTrig.omega = omega;
	_trigSumInit(&flowDensityTrigSum, &timeTrig, "flow density");
	_trigSumInit(&waterTempTrigSum, &timeTrig, "water temperature");
	_trigSumInit(&midDeepTempTrigSum, &timeTrig, "middle deep temperature");
	_trigSumInit(&eighthDeepTempTrigSum, &timeTrig, "eighth deep temperature");
#ifdef DEBUG_SUPERV_HOSE1DEEP
	_trigSumInit(&pos_hose1deep_av_TempTrigSum, &timeTrig, "pos hose1deep,av temperature");
	_trigSumInit(&pos_hose1deep_width_m1_TempTrigSum, &timeTrig, "pos hose1deep,width-1 temperature");
	_trigSumInit(&pos_hose1deep_mid_TempTrigSum, &timeTrig, "pos hose1deep,middle temperature");
	_trigSumInit(&pos_hose1deep_nh_TempTrigSum, &timeTrig, "pos hose1deep,near hose temperature");
	_trigSumInit(&pos_hose1deep_d2_av_TempTrigSum, &timeTrig, "pos hose1deep/2,av temperature");
#endif
	_printTime(_resf);
	if (_verbose)
		fprintf(_resf, "deep=%u, period=%u, flowDensityU=%.16e, earthLambda=%.16e, earthVHeatCap=%.16e, earthDiffusivity=%.16e\n", 
				deep, period, flowDensityU, earthLambda, earthVHeatCap, earthDiffusivity);
	{
		/* 1D simulation		*/
		const double	tempEvolCoef2 = 1.0 / (2.0 + invEarthDiffusivity) * CORR_MODERATION;	// coefficient of delta temperature to compute next temperature in earth cell
		unsigned	temp1DTSz;		// cell temperature table size
		double		*tempCT;		// current time cell temperature table
		double		*tempPT;		// previous time cell temperature table
		double		*tempCPMax;		// point to the last item of current time cell temperature
		unsigned	cycleIterLim = CYCLENB_BEFORE_1D_STEADY + (!hoseNb);

		/* init */
		temp1DTSz = deep + 1;			// +1: room cell before wall
		temp1DT = (double *) calloc(temp1DTSz, sizeof(double));
		assert(temp1DT != NULL);
		tempCT = temp1DT + 1;
		tempPT = (double *) calloc(deep, sizeof(double));
		assert(tempPT != NULL);
		tempCPMax = tempCT + deep -1;

		/* simulation loops */
		for (unsigned cycleIter = 0; cycleIter < cycleIterLim; cycleIter++) {
			unsigned	timeCpt;		// time counter in the cycle
			unsigned	convCpt;		// counter during convergence
			unsigned	convCptLim = CONV_ITER_NUMBER(deep);	// counter limit
			double		convErr;		// quadratic error during convergence
			double		convErrLim = CONV_ERR_LIM(1.0, deep, 512);	// errors limit
			double		signalAmplitude;	// before steady, it is gradually increasing

			_trigSumReset(&flowDensityTrigSum);
			_trigSumReset(&midDeepTempTrigSum);
			_trigSumReset(&eighthDeepTempTrigSum);
#ifdef DEBUG_SUPERV_HOSE1DEEP
			_trigSumReset(&pos_hose1deep_av_TempTrigSum);
			_trigSumReset(&pos_hose1deep_d2_av_TempTrigSum);
#endif

			for (timeCpt = 0; timeCpt < period; timeCpt++) {
				memcpy(tempPT, tempCT, sizeof(double) * deep);
				/* new angle */
				_trigSet(&timeTrig, timeCpt);

				/*
				There is a side effect of the starting of the sinusoid, resulting in a temperature offset in the depth of the wall.
				The ways to reduce this side effect:
				- increase the number of cycles before steady
				- attack the sinusoid with a ramp.
				*/
				signalAmplitude = 1.0;
				if (cycleIter < CYCLENB_AMPLITUDE_INCR) {
					double frac = (1.0 * (cycleIter * period + timeCpt)) / (CYCLENB_AMPLITUDE_INCR * period);
					/* If we found a better algorithm... write this here
					 * Already tried: 0.5 * (1 - cos(frac * M_PI)), sin(frac * M_PI / 2)
					*/
					signalAmplitude = frac;
				}
				tempCT[-1] = signalAmplitude * timeTrig.sinus;

				/* convergence */
				convCpt = 0;
				do {
					double		*tempCP;	// current time cell temperature pointer
					double		*tempPP;	// previous time cell temperature pointer
					double		tempDelta;	// sum of diff with neighbors
					double		tempCurr;	// current time cell temperature
					double		corr;		// noted temperature correction

					convErr = 0.0;
					int	cDeep = 1;		// deeper cell offset
					tempCP = tempCT;
					tempPP = tempPT;

					// first cell has distance half with its shallower neighbor
					tempCurr = *tempCP;
					tempDelta = 	((tempCP[-1] - tempCurr) *2 + (tempCP[cDeep] - tempCurr)) 
							+ ((*tempPP - tempCurr) * invEarthDiffusivity);
					*tempCP = (tempDelta * tempEvolCoef2) + tempCurr;
					corr = *tempCP - tempCurr;
					convErr += (corr * corr);

					// other cells
					for (tempCP++, tempPP++; tempCP <= tempCPMax; tempCP++, tempPP++) {
						if (tempCP >= tempCPMax)
							// last cell, with only one neighbor
							// the second neighbors will be itself => no delta temperature
							cDeep = 0;
						tempCurr = *tempCP;
						tempDelta = 	((tempCP[-1] - tempCurr) + (tempCP[cDeep] - tempCurr)) 
								+  ((*tempPP - tempCurr) * invEarthDiffusivity);
						*tempCP = (tempDelta * tempEvolCoef2) + tempCurr;
						corr = *tempCP - tempCurr;
						convErr += (corr * corr);
					}

					if (_stop)
						return 1;
				} while (++convCpt < convCptLim && convErr >= convErrLim);
				if (_verbose) {
					if (convCpt >= convCptLim) {
						printf("cycleIter=%u, timeCpt=%u/%u, convCpt=%u\n", cycleIter, timeCpt, period, convCpt);
						fprintf(_resf, "convCpt>=lim, cycleIter=%u, timeCpt=%u/%u\n", cycleIter, timeCpt, period);
						_dumpTemp(_resf, tempCT, deep);
					}
					if (!(timeCpt % 100))
						printf("cycleIter=%u, timeCpt=%u/%u, convCpt=%u\n", cycleIter, timeCpt, period, convCpt);
					if (timeCpt == period/2-1) {
						fprintf(_resf, "\n");
						_dumpTemp(_resf, tempCT, deep);
					}
				}
				// results for this angle
				_trigSumAdd(&flowDensityTrigSum, earthLambda/flowDensityU * (tempCT[-1] - *tempCT) *2);
				_trigSumAdd(&midDeepTempTrigSum, (tempCT[deep/2] + tempCT[deep/2 -1]) /2);
				_trigSumAdd(&eighthDeepTempTrigSum, (tempCT[deep/8] + tempCT[deep/8 -1]) /2);
#ifdef DEBUG_SUPERV_HOSE1DEEP
				if (!!DEBUG_SUPERV_HOSE1DEEP) {
					_trigSumAdd(&pos_hose1deep_av_TempTrigSum, (tempCT[DEBUG_SUPERV_HOSE1DEEP] + tempCT[DEBUG_SUPERV_HOSE1DEEP -1]) /2);
					_trigSumAdd(&pos_hose1deep_d2_av_TempTrigSum, (tempCT[DEBUG_SUPERV_HOSE1DEEP /2] + tempCT[DEBUG_SUPERV_HOSE1DEEP /2 -1]) /2);
				}
#endif
			}

			/* exploitation of the results for the whole cycle (one period) */
			if (_verbose || cycleIter >= CYCLENB_BEFORE_1D_STEADY) {
				/* write in the result file */
				fprintf(_resf, "\n");
				_printTime(_resf);
				fprintf(_resf, "cycleIter: %u\n", cycleIter);
				// module and phase
				_trigPrint(_resf, &flowDensityTrigSum);
				_trigPrint(_resf, &midDeepTempTrigSum);
				_trigPrint(_resf, &eighthDeepTempTrigSum);
#ifdef DEBUG_SUPERV_HOSE1DEEP
				if (!!DEBUG_SUPERV_HOSE1DEEP) {
					_trigPrint(_resf, &pos_hose1deep_av_TempTrigSum);
					_trigPrint(_resf, &pos_hose1deep_d2_av_TempTrigSum);
				}
#endif
				// dump tempCT
				_dumpTemp(_resf, tempCT, deep);
			}
		}

		/* free function data */
		free(tempPT);
	}
	if (!hoseNb) {
		/* free function data */
		free(temp1DT);
		return 0;
	}
	if (_verbose)
		fprintf(_resf, "\n\nwidth=%u, waterVHeatCap=%.16e, waterDiffusivity=%.16e, hoseNb=%u, hoseRadius=%u, hose1Deep=%u, hoseInterval=%u\n", 
				width, waterVHeatCap, waterDiffusivity, hoseNb, hoseRadius, hose1Deep, hoseInterval);
	{
		/* 2D simulation		*/
		const double	tempEvolCoef4 = 1.0 / (4.0 + invEarthDiffusivity) * CORR_MODERATION;	// coefficient of delta temperature to compute next temperature in earth cell
		double		tempCWater;		// current time water temperature
		double		tempPWater;		// previous time water temperature
		double		*temp2DT;		// cell temperature 2D table
		unsigned	temp2DTSz;		// temp2DT size
		double		*tempCT;		// current time cell temperature table
		double		*tempPT;		// previous time cell temperature table
		unsigned	tempTSz;		// tempPT size
		double		*lineCPLim;		// point after the last line of current time cell temperature
		double		*lineCPMax;		// point to the last line of current time cell temperature
		// one line for one deep
		unsigned	lineSz;			// line (in cell temperature table) size
		EarthLineDesc	*earthLineDescT;	// EarthLineDesc table
		EarthLineDesc	*earthLineDescTLim;	// point after EarthLineDesc table
		HosePos		*hosePosT;		// HosePos table
		Pos		*hoseNeighborT;		// neighbor position of hoses table
		Pos		*hoseNeighborPLim;	// point after hoseNeighborT
		unsigned	hoseNeighborTSz;	// hoseNeighborT size
		unsigned	waterCellNb;		// number of cells for water, used for calculate convErr

		unsigned	cycleIterLim = CYCLENB_BEFORE_2D_STEADY +1;

		/* init */
		lineSz = width + 1;			// +1: for temporary variable (the neighbor of an extreme)
		tempTSz = deep * lineSz;
		temp2DTSz = (deep + 1) * lineSz;	// +1: one line for room cells, before the wall
		temp2DT = (double *) calloc(temp2DTSz, sizeof(double));
		assert(temp2DT != NULL);
		tempCT = temp2DT + lineSz;		// this table starts after room cells
		tempPT = (double *) calloc(tempTSz, sizeof(double));
		assert(tempPT != NULL);
		lineCPLim = tempCT + tempTSz;
		lineCPMax = lineCPLim - lineSz;
		// init 2D tempCT values with 1D tempCT
		for (unsigned i = 0; i < deep; i++)
			_lineSet(tempCT + i * lineSz, width, temp1DT[i+1]);
		free(temp1DT);				// no need of temp1DT

		earthLineDescT = (EarthLineDesc *) calloc(deep, sizeof(EarthLineDesc));
		assert(earthLineDescT != NULL);
		earthLineDescTLim = earthLineDescT + deep;
		for (EarthLineDesc *descP = earthLineDescT; descP < earthLineDescTLim; descP++) {
			descP->begIx = 0;
			descP->limIx = width;
		}

		hosePosT = (HosePos *) calloc(hoseNb, sizeof(HosePos));
		assert(hosePosT != NULL);
		hosePosT[0].deepIx = hose1Deep;
		hosePosT[0].right = false;
		for (unsigned i = 1; i < hoseNb; i++) {
			hosePosT[i].deepIx = hosePosT[i-1].deepIx + hoseInterval;
			hosePosT[i].right = !hosePosT[i-1].right;
		}

		hoseNeighborTSz = hoseNb * (4 * hoseRadius);
		hoseNeighborT = (Pos *) calloc(hoseNeighborTSz, sizeof(Pos));
		assert(hoseNeighborT != NULL);
		hoseNeighborPLim = hoseNeighborT + hoseNeighborTSz;
		{
			unsigned	*roundLimT;	// round limit, line by line, table

			roundLimT = (unsigned *) calloc(hoseRadius, sizeof(unsigned));
			assert(roundLimT != NULL);
			{
				unsigned	sqRadius = hoseRadius * hoseRadius *4;
				unsigned	pos;

				for (unsigned i = 0; i < hoseRadius; i++) {
					pos = i + i +1;
					roundLimT[i] = ((unsigned)sqrt(sqRadius - pos * pos) +1) /2;
				}
			}

			waterCellNb = 0;
			// number of cells for one hose within the "width" window
			for (unsigned i = 0; i < hoseRadius; i++)
				waterCellNb += roundLimT[i] *2;
			// total number
			waterCellNb *= hoseNb;

			int		deepIx;
			int		deepRel;		// relative pos
			int		lateralAbs, deepAbs;	// absolute relative pos
			Pos		*hoseNeighborP = hoseNeighborT;
			for (unsigned i = 0; i < hoseNb; i++) {
				// less deep neighbors
				for (lateralAbs = 0; lateralAbs < hoseRadius; lateralAbs++) {
					deepAbs = roundLimT[lateralAbs] +1;
					hoseNeighborP->deepIx = hosePosT[i].deepIx - deepAbs;
					hoseNeighborP->lateralIx = (hosePosT[i].right) ? width - lateralAbs -1 : lateralAbs;
					hoseNeighborP++;
				}
				// lateral neighbors
				for (deepRel = -hoseRadius; deepRel < (int)hoseRadius; deepRel++) {
					deepAbs = (deepRel < 0) ? -deepRel -1 : deepRel;
					lateralAbs = roundLimT[deepAbs];
					hoseNeighborP->deepIx = deepIx = hosePosT[i].deepIx + deepRel;
					hoseNeighborP->lateralIx = (hosePosT[i].right) ? width - lateralAbs -1 : lateralAbs;
					hoseNeighborP++;
					// correcting earthLineDescT
					if (hosePosT[i].right)
						earthLineDescT[deepIx].limIx -= lateralAbs;
					else
						earthLineDescT[deepIx].begIx += lateralAbs;
				}
				// deeper neighbors
				for (lateralAbs = hoseRadius -1; lateralAbs >= 0; lateralAbs--) {
					deepAbs = roundLimT[lateralAbs];
					hoseNeighborP->deepIx = hosePosT[i].deepIx + deepAbs;
					hoseNeighborP->lateralIx = (hosePosT[i].right) ? width - lateralAbs -1 : lateralAbs;
					hoseNeighborP++;
				}
			}

			free(roundLimT);
		}
		const double	tempEvolCoefWater = tempEvolCoef4 * earthVHeatCap / (waterVHeatCap * waterCellNb);	// coefficient of delta temperature to compute next temperature in the hoses
		const double	invWaterDiffusivity = waterCellNb / waterDiffusivity;	// inverse of water diffusivity by water cell number

		// water temperature init
		tempCWater = 0.0;
		/* warning: simplified computing of water temperature
		we content ourselves with the temperature of the center of the hose
		*/
		for (unsigned i = 0; i < hoseNb; i++) {
			tempCWater += tempCT[hosePosT[i].deepIx * lineSz];
		}
		tempCWater /= hoseNb;
		_tempWaterSet(earthLineDescT, tempCT, lineCPLim, lineSz, tempCWater);

		/* simulation loops */
		for (unsigned cycleIter = 0; cycleIter < cycleIterLim; cycleIter++) {
			unsigned	timeCpt;		// time counter in the cycle
			unsigned	convCpt;		// counter during convergence
			int		convExtra;		// additional convergence iterations
			unsigned	convCptLim = CONV_ITER_NUMBER(tempTSz);	// counter limit
			double		convErr;		// quadratic error during convergence
			double		convErrLim = CONV_ERR_LIM(1.0, tempTSz, 2048);	// errors limit

			_trigSumReset(&flowDensityTrigSum);
			_trigSumReset(&waterTempTrigSum);
			_trigSumReset(&midDeepTempTrigSum);
			_trigSumReset(&eighthDeepTempTrigSum);
#ifdef DEBUG_SUPERV_HOSE1DEEP
			_trigSumReset(&pos_hose1deep_av_TempTrigSum);
			_trigSumReset(&pos_hose1deep_width_m1_TempTrigSum);
			_trigSumReset(&pos_hose1deep_mid_TempTrigSum);
			_trigSumReset(&pos_hose1deep_nh_TempTrigSum);
			_trigSumReset(&pos_hose1deep_d2_av_TempTrigSum);
#endif

			for (timeCpt = 0; timeCpt < period; timeCpt++) {
				tempPWater = tempCWater;
				memcpy(tempPT, tempCT, sizeof(double) * tempTSz);
				/* new angle */
				_trigSet(&timeTrig, timeCpt);
				_lineSet(tempCT -lineSz, width, timeTrig.sinus);

				/* convergence */
				convCpt = 0;
				convExtra = 1+5;
				do {
					Pos		*hoseNeighborP;
					EarthLineDesc	*descP;
					double		*lineCP;	// current time line temperature pointer
					double		*linePP;	// previous time line temperature pointer
					double		*tempCP;	// current time cell temperature pointer
					double		*tempPP;	// previous time cell temperature pointer
					double		*tempCPLim;
					double		tempDelta;	// sum of diff between neighbors
					double		tempCurr;	// current time cell temperature
					double		tempSv;		// current time left neighbor temperature, before update
					double		corr;		// noted temperature correction
					int		begIx;

					convErr = 0.0;

					// temperature for water
					tempDelta = 0.0;
					for (hoseNeighborP = hoseNeighborT; hoseNeighborP < hoseNeighborPLim; hoseNeighborP++) {
						tempDelta += tempCT[hoseNeighborP->deepIx * lineSz + hoseNeighborP->lateralIx] - tempCWater;
					}
					tempDelta += (tempPWater - tempCWater) * invWaterDiffusivity;
					corr = -tempCWater;
					tempCWater = (tempDelta * tempEvolCoefWater) + tempCWater;
					corr += tempCWater;
					convErr += (corr * corr) * waterCellNb;

					// temperature for earth
					int	lShal = -(int)lineSz;	// shallower line offset
					int	lDeep = lineSz;		// deeper line offset
					descP = earthLineDescT;
					lineCP = tempCT;
					linePP = tempPT;

					// first line has distance half with its shallower neighbor
					begIx = descP->begIx;
					lineCP[width] = lineCP[width -1];	// extremum neighbor
					lineCP[-1] = *lineCP;			// extremum neighbor
					tempSv = lineCP[begIx -1];
					for (tempCP = lineCP + begIx, tempCPLim = lineCP + descP->limIx, tempPP = linePP + begIx;
							tempCP < tempCPLim;
							tempCP++, tempPP++) {
						tempCurr = *tempCP;
						tempDelta = 	((tempSv - tempCurr) + (tempCP[1] - tempCurr)) 
								+ ((tempCP[lShal] - tempCurr) *2 + (tempCP[lDeep] - tempCurr)) 
								+ ((*tempPP - tempCurr) * invEarthDiffusivity);
						*tempCP = (tempDelta * tempEvolCoef4) + tempCurr;
						tempSv = tempCurr;
						corr = *tempCP - tempCurr;
						convErr += (corr * corr);
					}

					// other lines
					for (descP++, lineCP += lineSz, linePP += lineSz; lineCP < lineCPLim; descP++, lineCP += lineSz, linePP += lineSz) {
						if (lineCP >= lineCPMax)
							// last line, with only 3 neighbors
							// the fourth neighbors will be itself => no delta temperature
							lDeep = 0;
						begIx = descP->begIx;
						lineCP[width] = lineCP[width -1];	// extremum neighbor
						lineCP[-1] = *lineCP;			// extremum neighbor
						tempSv = lineCP[begIx -1];
						for (tempCP = lineCP + begIx, tempCPLim = lineCP + descP->limIx, tempPP = linePP + begIx;
								tempCP < tempCPLim;
								tempCP++, tempPP++) {
							tempCurr = *tempCP;
							tempDelta = 	((tempSv - tempCurr) + (tempCP[1] - tempCurr))
									+ ((tempCP[lShal] - tempCurr) + (tempCP[lDeep] - tempCurr)) 
									+ ((*tempPP - tempCurr) * invEarthDiffusivity);
							*tempCP = (tempDelta * tempEvolCoef4) + tempCurr;
							tempSv = tempCurr;
							corr = *tempCP - tempCurr;
							convErr += (corr * corr);
#ifndef NDEBUG
							if (convCpt > 50 && !!corr)
								printf("convCpt=%u, line=%ld, lateral=%ld, corr=%e\n", convCpt, (lineCP-tempCT)/lineSz, tempCP-lineCP, corr);
#endif
						}
					}

					// set tempCT with water temperature
					_tempWaterSet(earthLineDescT, tempCT, lineCPLim, lineSz, tempCWater);

					if (_stop)
						return 1;
					if (convErr < convErrLim) {
						if (!convErr)
							convExtra = 0;
						else
							convExtra--;
					}
				} while (++convCpt < convCptLim && convExtra > 0);
				if (_verbose) {
					if (convCpt >= convCptLim) {
						printf("cycleIter=%u, timeCpt=%u/%u, convCpt=%u\n", cycleIter, timeCpt, period, convCpt);
						fprintf(_resf, "convCpt>=lim, cycleIter=%u, timeCpt=%u/%u\n", cycleIter, timeCpt, period);
						_dumpTempAverage(_resf, tempCT, lineCPLim, lineSz);
					}
					if (!(timeCpt % 100))
						printf("cycleIter=%u, timeCpt=%u/%u, convCpt=%u\n", cycleIter, timeCpt, period, convCpt);
					if (timeCpt == period/2-1) {
						fprintf(_resf, "\n");
						_dumpTempAverage(_resf, tempCT, lineCPLim, lineSz);
#if 0
						//fprintf(_resf, "\n");
						//_dumpTemp(_resf, tempCT + (hosePosT[0].deepIx - hoseRadius -1) * lineSz, width);
						fprintf(_resf, "\n");
						_dumpTemp(_resf, tempCT + hosePosT[0].deepIx * lineSz, width);
#endif
					}
				}
				// results for this angle
				_trigSumAdd(&flowDensityTrigSum, earthLambda/flowDensityU * _lineDiffAverage(tempCT -lineSz, lineSz, width) *2);
				_trigSumAdd(&waterTempTrigSum, tempCWater);
				_trigSumAdd(&midDeepTempTrigSum, _lineMidAverage(tempCT + (deep/2 -1) * lineSz, lineSz, width));
				_trigSumAdd(&eighthDeepTempTrigSum, _lineMidAverage(tempCT + (deep/8 -1) * lineSz, lineSz, width));
#ifdef DEBUG_SUPERV_HOSE1DEEP
				_trigSumAdd(&pos_hose1deep_av_TempTrigSum, _lineMidAverage(tempCT + (DEBUG_SUPERV_HOSE1DEEP -1) * lineSz, lineSz, width));
				_trigSumAdd(&pos_hose1deep_width_m1_TempTrigSum, (tempCT[(DEBUG_SUPERV_HOSE1DEEP -1) * lineSz + width -1] + tempCT[(DEBUG_SUPERV_HOSE1DEEP) * lineSz + width -1])/2);
				_trigSumAdd(&pos_hose1deep_mid_TempTrigSum, (tempCT[(DEBUG_SUPERV_HOSE1DEEP -1) * lineSz + width/2] + tempCT[(DEBUG_SUPERV_HOSE1DEEP) * lineSz + width/2])/2);
				_trigSumAdd(&pos_hose1deep_nh_TempTrigSum, (tempCT[(DEBUG_SUPERV_HOSE1DEEP -1) * lineSz + hoseRadius +2] + tempCT[(DEBUG_SUPERV_HOSE1DEEP) * lineSz + hoseRadius +2])/2);
				_trigSumAdd(&pos_hose1deep_d2_av_TempTrigSum, _lineMidAverage(tempCT + (DEBUG_SUPERV_HOSE1DEEP/2 -1) * lineSz, lineSz, width));
#endif
			}

			/* exploitation of the results for the whole cycle (one period) */
			if (_verbose || cycleIter >= CYCLENB_BEFORE_2D_STEADY) {
				/* write in the result file */
				fprintf(_resf, "\n");
				_printTime(_resf);
				fprintf(_resf, "cycleIter: %u\n", cycleIter);
				// module and phase
				_trigPrint(_resf, &flowDensityTrigSum);
				_trigPrint(_resf, &waterTempTrigSum);
				_trigPrint(_resf, &midDeepTempTrigSum);
				_trigPrint(_resf, &eighthDeepTempTrigSum);
#ifdef DEBUG_SUPERV_HOSE1DEEP
				_trigPrint(_resf, &pos_hose1deep_av_TempTrigSum);
				_trigPrint(_resf, &pos_hose1deep_width_m1_TempTrigSum);
				_trigPrint(_resf, &pos_hose1deep_mid_TempTrigSum);
				_trigPrint(_resf, &pos_hose1deep_nh_TempTrigSum);
				_trigPrint(_resf, &pos_hose1deep_d2_av_TempTrigSum);
#endif
				// dump tempCT
				_dumpTempAverage(_resf, tempCT, lineCPLim, lineSz);
#if 0
				//fprintf(_resf, "\n");
				//_dumpTemp(_resf, tempCT + (hosePosT[0].deepIx - hoseRadius -1) * lineSz, width);
				fprintf(_resf, "\n");
				_dumpTemp(_resf, tempCT + hosePosT[0].deepIx * lineSz, width);
#endif
			}


		}

		/* free function data */
		free(hoseNeighborT);
		free(hosePosT);
		free(earthLineDescT);
		free(tempPT);
		free(temp2DT);
	}

	return 0;
}

static	void	_print_option_err(const char *msg) {
	fprintf(stderr, "option: incorrect value for %s.\n", msg);
}

static	void	_print_param_err(const char *msg) {
	fprintf(stderr, "param.: incorrect value for %s.\n", msg);
}

static	void	_print_usage() {
	printf("Usage:  %s -h\n", _progname);
	printf("      | %s OPTIONS... --hoseNb=0 result_file\n", _progname);
	printf("      | %s OPTIONS... HOSE_OPTIONS... result_file\n", _progname);
	printf("OPTIONS:\n");
	printf(" -h                       help\n");
	printf(" -v                       verbose\n");
	printf(" --lengthU=real_num       length internal unit = cell size\n");
	printf(" --timeU=real_num         time internal unit  = simulation time step,\n");
	printf("                          one second, or a multiple\n");
	printf(" --width=real_num         width of the area, multiple of lengthU, > 0, \n");
	printf("                          significant if hoseNb >0\n");
	printf(" --deep=real_num          deep of the area, multiple of lengthU, > 0\n");
	printf(" --period=real_num        period of heat wave, multiple of timeU, > 0\n");
	printf(" --earthLambda=real_num   earth lambda, > 0.0\n");
	printf(" --earthVHeatCap=real_num earth volumetric heat capacity, > 0\n");
	printf(" --waterVHeatCap=real_num water volumetric heat capacity, > 0, \n");
	printf("                          significant if hoseNb >0\n");
	printf(" --hoseNb=integer         number of hoses\n");
	printf("HOSE_OPTIONS, significant if hoseNb >0:\n");
	printf(" --hoseRadius=real_num    radius of hose, multiple of lengthU, > 0\n");
	printf(" --hose1Deep=real_num     deep of the first hose, multiple of lengthU,\n");
	printf("                          > hoseRadius\n");
	printf(" --hoseInterval=real_num  deep interval between hoses, multiple of lengthU, \n");
	printf("                       	  > 2*hoseRadius, significant if hoseNb >1\n");
	printf("\n");
	printf("                          hose1Deep + hoseInterval*(hoseNb-1) + hoseRadius\n");
	printf("                          < deep\n");
	exit(1);
}

int	main(int argc, char *argv[]) {
	double		lengthU = 1.0;		// length internal unit = cell size (less than a meter)
	double		timeU = 1.0;		// time internal unit (one second, or a multiple)
	double		width = 0.0;		// width of the map, multiple of lengthU
	double		deep = 0.0;		// deep of the map, multiple of lengthU
	double		period = 0.0;		// period of heat wave, multiple of timeU
	double		earthLambda = 0.0;	// earth lambda
	double		earthVHeatCap = 0.0;	// earth volumetric heat capacity
	double		waterVHeatCap = 0.0;	// water volumetric heat capacity
	uint32_t	hoseNb = 0;		// number of hoses
	double		hoseRadius = 0.0;	// radius of hose, multiple of lengthU
	double		hose1Deep = 0.0;	// deep of the first hose, multiple of lengthU
	double		hoseInterval = 0.0;	// deep interval between hoses, multiple of lengthU

	const char * 	resfile = NULL;		// result file name
	int		ret;			// return of main function

	_progname = _basename(argv[0]);
	/* parameters getting */
	while (1) {
		int		option;
		int		option_index = 0;
		static struct option long_options[] = {
			{"lengthU",		required_argument,	0, 'L' },
			{"timeU",		required_argument,	0, 'T' },
			{"width",		required_argument,	0, 'w' },
			{"deep",		required_argument,	0, 'd' },
			{"period",		required_argument,	0, 'p' },
			{"earthLambda",		required_argument,	0, 'l' },
			{"earthVHeatCap",	required_argument,	0, 'c' },
			{"waterVHeatCap",	required_argument,	0, 'C' },
			{"hoseNb",		required_argument,	0, 'N' },
			{"hoseRadius",		required_argument,	0, 'R' },
			{"hose1Deep",		required_argument,	0, 'D' },
			{"hoseInterval",	required_argument,	0, 'I' },
			{0,		0,			0, 0 }
		};
		option = getopt_long(argc, argv, "hv", long_options, &option_index);
		if (option == -1)
			break;
		switch (option) {
		case 'h' :
			_print_usage();
			break;
		case 'v' :
			_verbose = true;
			break;
		case 'L' :
			if (sscanf(optarg, "%lf", &lengthU) != 1)
				_print_option_err("lengthU");
			break;
		case 'T' :
			if (sscanf(optarg, "%lf", &timeU) != 1)
				_print_option_err("timeU");
			break;
		case 'w' :
			if (sscanf(optarg, "%lf", &width) != 1)
				_print_option_err("width");
			break;
		case 'd' :
			if (sscanf(optarg, "%lf", &deep) != 1)
				_print_option_err("deep");
			break;
		case 'p' :
			if (sscanf(optarg, "%lf", &period) != 1)
				_print_option_err("period");
			break;
		case 'l' :
			if (sscanf(optarg, "%lf", &earthLambda) != 1)
				_print_option_err("earthLambda");
			break;
		case 'c' :
			if (sscanf(optarg, "%lf", &earthVHeatCap) != 1)
				_print_option_err("earthVHeatCap");
			break;
		case 'C' :
			if (sscanf(optarg, "%lf", &waterVHeatCap) != 1)
				_print_option_err("waterVHeatCap");
			break;
		case 'N' :
			if (sscanf(optarg, "%u", &hoseNb) != 1)
				_print_option_err("hoseNb");
			break;
		case 'R' :
			if (sscanf(optarg, "%lf", &hoseRadius) != 1)
				_print_option_err("hoseRadius");
			break;
		case 'D' :
			if (sscanf(optarg, "%lf", &hose1Deep) != 1)
				_print_option_err("hose1Deep");
			break;
		case 'I' :
			if (sscanf(optarg, "%lf", &hoseInterval) != 1)
				_print_option_err("hoseInterval");
			break;
		default:
			fprintf(stderr, "%s: bad option '%c'\n", _progname, option);
		}
	}
	if (optind < argc)
		resfile = argv[optind];
	if (!resfile) {
		fprintf(stderr, "%s: missing result file.\n", _progname);
		_print_usage();
	}
	if (_verbose) {
		printf("%s lengthU=%.16e, timeU=%.16e, width=%.16e, deep=%.16e, period=%.16e, earthLambda=%.16e, earthVHeatCap=%.16e, waterVHeatCap=%.16e, hoseNb=%u, hoseRadius=%.16e, hose1Deep=%.16e, hoseInterval=%.16e, result_file=%s\n",
		       _progname, lengthU, timeU, width, deep, period, earthLambda, earthVHeatCap, waterVHeatCap, hoseNb, hoseRadius, hose1Deep, hoseInterval, resfile);
	}
	if (deep <= lengthU) {
		_print_param_err("deep");
		_print_usage();
	}
	if (period <= timeU) {
		_print_param_err("period");
		_print_usage();
	}
	if (earthLambda <= 0.0) {
		_print_param_err("earthLambda");
		_print_usage();
	}
	if (earthVHeatCap <= 0.0) {
		_print_param_err("earthVHeatCap");
		_print_usage();
	}
	if (!!hoseNb) {
		if (width <= lengthU) {
			_print_param_err("width");
			_print_usage();
		}
		if (waterVHeatCap <= 0.0) {
			_print_param_err("waterVHeatCap");
			_print_usage();
		}
		if (hoseRadius <= lengthU) {
			_print_param_err("hoseRadius");
			_print_usage();
		}
		if (hose1Deep <= hoseRadius) {
			_print_param_err("hose1Deep");
			_print_usage();
		}
		if (hoseNb > 1 && hoseInterval <= 2*hoseRadius) {
			_print_param_err("hoseInterval");
			_print_usage();
		}
		if (hose1Deep + hoseInterval*(hoseNb-1) + hoseRadius >= deep) {
			_print_param_err("hose1Deep + hoseInterval*(hoseNb-1) + hoseRadius");
			_print_usage();
		}
	}
	/* result file */
	if ((_resf = fopen(resfile,"w")) == NULL) {
		fprintf(stderr,"Cannot open result file\n");
		return -1;
	}
	/* simulation */
	ret = _simul(
		      FLOAT_TO_INT(width/lengthU),
		      FLOAT_TO_INT(deep/lengthU),
		      FLOAT_TO_INT(period/timeU),
		      lengthU * lengthU * timeU,
		      earthLambda * lengthU * timeU,
		      earthVHeatCap * lengthU * lengthU * lengthU,
		      waterVHeatCap * lengthU * lengthU * lengthU,
		      hoseNb,
		      FLOAT_TO_INT(hoseRadius/lengthU),
		      FLOAT_TO_INT(hose1Deep/lengthU),
		      FLOAT_TO_INT(hoseInterval/lengthU));
	/* end */
	fclose(_resf);
	return ret;
}
