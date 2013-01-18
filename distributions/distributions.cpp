//************************************************************************************************************************************************
//
//    Copyright (c) 2011 David Briant - see https://github.com/DangerMouseB
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Lesser General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//************************************************************************************************************************************************


#include <ymath.h>			// _Nan, _Inf
#include <float.h>				// _isnan, _isfinite
#include <math.h>			// fabs, exp, sqrt, log
//#include "stdlib.h"				// __min
#include <algorithm>			// std::min

#define NOMINMAX			// see http://www.suodenjoki.dk/us/archive/2010/min-max.htm
using namespace std;

#define M_SQRT2PI 2.50662827463100050242		// see http://www.jquantlib.org/sites/jquantlib/apidocs/0.1.2/src-html/org/jquantlib/math/Constants.html

extern "C" double __stdcall CN_Hart(double x) {
    // from BETTER APPROXIMATIONS TO CUMULATIVE NORMAL FUNCTIONS By GRAEME WEST
    double xabs, result;
	if (_isnan(x)) return _Nan._Double;
	if (!_finite(x)) return (x < 0 ? 0.0 : 1.0);
    xabs = fabs(x);
    if (xabs > 37.0) return 0.0;
	if (xabs < 7.07106781186547) {
		result = exp(-xabs * xabs * 0.5)
			*  (((((((3.52624965998911E-02 * xabs + 0.700383064443688)) * xabs + 6.37396220353165) * xabs + 33.912866078383) * xabs + 112.079291497871) * xabs + 221.213596169931) * xabs + 220.206867912376)
			/ (((((((8.83883476483184E-02 * xabs + 1.75566716318264) * xabs + 16.064177579207) * xabs + 86.7807322029461) * xabs + 296.564248779674) * xabs + 637.333633378831) * xabs + 793.826512519948) * xabs + 440.413735824752);
	} else {
		result = exp(-xabs * xabs * 0.5) / (xabs + 1.0 / (xabs + 2.0 / (xabs + 3.0 / (xabs + 4.0 / (xabs + 0.65))))) / 2.506628274631;
	}
    if (x < 0.0) {
		return result;
	} else {
		return 1.0 - result;
	}
}

extern "C" int __stdcall CN_Harts(double const *ix, double *op, int size) {
	for (int i = 0; i < size; i++) {
		*op++ = CN_Hart(*ix++);
	}
	return 0;
}

extern "C" double __stdcall InvCN_Acklam(double p) {
	// http://home.online.no/~pjacklam/notes/invnorm/index.html
	// http://home.online.no/~pjacklam/notes/invnorm/impl/herrero/inversecdf.txt

	const double a[6] = {
		-3.969683028665376e+01,  2.209460984245205e+02,
		-2.759285104469687e+02,  1.383577518672690e+02,
		-3.066479806614716e+01,  2.506628277459239e+00
	};
	const double b[5] = {
		-5.447609879822406e+01,  1.615858368580409e+02,
		-1.556989798598866e+02,  6.680131188771972e+01,
		-1.328068155288572e+01
	};
	const double c[6] = {
		-7.784894002430293e-03, -3.223964580411365e-01,
		-2.400758277161838e+00, -2.549732539343734e+00,
		4.374664141464968e+00,  2.938163982698783e+00
	};
	const double d[4] = {
		7.784695709041462e-03,  3.224671290700398e-01,
		2.445134137142996e+00,  3.754408661907416e+00
	};

	double q, t, u;

	if (_isnan(p) || p > 1.0 || p < 0.0) return _Nan._Double;
	if (p == 0.0) return -_Inf._Double;
	if (p == 1.0) return  _Inf._Double;
	q = min(p, 1-p);
	if (q > 0.02425) {
		// Rational approximation for central region.
		u = q - 0.5;
		t = u*u;
		u = u*(((((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4])*t+a[5])
		/(((((b[0]*t+b[1])*t+b[2])*t+b[3])*t+b[4])*t+1);
	} else {
		// Rational approximation for tail region.
		t = sqrt(-2*log(q));
		u = (((((c[0]*t+c[1])*t+c[2])*t+c[3])*t+c[4])*t+c[5])
		/((((d[0]*t+d[1])*t+d[2])*t+d[3])*t+1);
	}
	// The relative error of the approximation has absolute value less
	// than 1.15e-9.  One iteration of Halley's rational method (third
	// order) gives full machine precision...
	t = CN_Hart(u) - q;
	t = t * M_SQRT2PI * exp(u * u * 0.5);			// f(u)/df(u)
	u = u - t / (1 + u * t * 0.5);						// Halley's method

	return (p > 0.5 ? -u : u);
}

extern "C" int __stdcall InvCN_Acklams(double const *ip, double *ox, int size) {
	for (int i = 0; i < size; i++) {
		*ox++ = InvCN_Acklam(*ip++);
	}
	return 0;
}

