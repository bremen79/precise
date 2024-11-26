#include <cmath>
#include <limits>
#include <cstdint>
#include <iostream>
#include <algorithm>

using namespace std;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

float objective(const float *g,
		const uint32_t length,
		const float bet) {
    float total = 0;

    for (int i=0; i < length; ++i) {
	total += log(1.0 + g[i]*bet);
    }
    
    return total;
}

float d_objective(const float *g,
		  const uint32_t length,
		  const float bet) {
    float total = 0;
    for (int i=0; i < length; ++i) {
	total += (g[i] / (1.0 + g[i]*bet));
    }
    return total;
}

float d2_objective(const float *g,
		   const uint32_t length,
		   const float bet) {
    float total = 0;
    for (int i=0; i < length; ++i) {
	total += - pow(g[i] / (1.0 + g[i]*bet), 2);
    }
    return total;
}

 
float newton_1d_bnd(const float *g,
		    const uint32_t length,
		    const float ax,
		    const float bx,
		    const float x_init,
		    float &x_out) {
    float deriv = numeric_limits<float>::max();
    float x_old = numeric_limits<float>::max();
    float x = x_init;

    while ( (abs(deriv) > 0.001) && (abs(x-x_old) > 0.01) ) {
	x_old = x;
	deriv = d_objective(g, length, x);
	if ((x == ax) && (deriv < 0)) {
	    break;
	}

	if ((x == bx) & (deriv > 0)) {
	    break;
	}

	float update = 0;
	if (abs(deriv) > 1e3) {
	    update = 0.01 * sgn(deriv);       	
	}
	else {
	    const float deriv2 = d2_objective(g, length, x);
	    update = -deriv/deriv2;
	}
	
	x += update;
	x = max(min(x, bx), ax);
    }

    const float fval = objective(g, length, x);
    x_out = x;

    return fval;
}

extern "C" float find_mean_max_log_wealth_constrained(const float *X,
						      const uint32_t n,
						      const float bmin,
						      const float bmax,
						      float &x_out) {

    float x_init = 0;
    float sum_log_W_star = 0;

    float update = 0;

    sum_log_W_star = newton_1d_bnd(X, n, bmin, bmax, x_init, x_out);
    return sum_log_W_star;
}
