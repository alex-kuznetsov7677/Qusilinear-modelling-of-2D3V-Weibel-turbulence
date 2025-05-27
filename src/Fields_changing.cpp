#include "my_variables.h"
#include "Fields_changing.h"
#include <vector>
#include <cmath>
#include <complex>
using std::complex;
using std::vector;

complex <double> SpeedIntegral(std::vector<complex<double> >& flka, double p11) {


	complex<double> k(0, 0);
	double BBnomer = 0;
	double BBreal = 0;
	double BBYrl = 0;
	for (int z = 0; z < setkaBB; z++) {
		double BBZnomerk = setkaBBkvadr * z;
		for (double i = 0; i < setkaBB; i++) {
			BBnomer = p11 * setkaBBkub + BBZnomerk + i * setkaBB;
			BBreal = -BBXmax + i * dBBX;
			for (double j = 0; j < setkaBB; j++) {
				BBYrl = -BBYmax + j * dBBY;
				k = k + C4 * flka[BBnomer + j] * BBYrl * dBBY * dBBX * dBBZ;
			}
		}
	}
	return k;
}

void WEIBEL_FIELDS(int& oblast, int& size, int& Ngarmonik,std::vector<double>& Kxvector,std::vector<double>& Kzvector, std::vector<complex<double> >& IEylast, std::vector<complex<double> >& bzlast,std::vector<complex<double> >& bxlast,
	std::vector<complex<double> >& IEy, std::vector<complex<double> >& bz,std::vector<complex<double> >& bx, std::vector<complex<double> >& IEy1, std::vector<complex<double> >& bz1,std::vector<complex<double> >& bx1,
std::vector<complex<double> >& IEy1last, std::vector<complex<double> >& bz1last, std::vector<complex<double> >& bx1last,std::vector<complex<double> >& fk, std::vector<complex<double> >& fkcel) {

	for (int p = 0; p < Ngarmonik/ size; p++) {
		int preal = p + oblast * Ngarmonik / size;
		bz1[p] = bz1last[p] -Kxvector[preal]*IEylast[p] * dt;
		bx1[p] = bx1last[p] +Kzvector[preal]*IEylast[p] * dt;
		complex <double> Integ = SpeedIntegral(fk, p);
		IEy[p] = IEylast[p] + dt * (Kxvector[preal] * bz1[p] -Kzvector[preal] * bx1[p]- I * Integ);
		IEy1[p] = (IEy[p] + IEylast[p]) / 2.;
		bz[p] = bzlast[p] - Kxvector[preal] * IEy1[p] * dt;
		bx[p] = bxlast[p] + Kzvector[preal] * IEy1[p] * dt;
	}
}
