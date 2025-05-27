#include "my_variables.h"
#include "Integrals.h"

#include <vector>
#include <cmath>

using std::complex;
using std::vector;

complex <double> ENERGY_Y(const std::vector<complex<double> >& fka) {

	complex<double> k(0, 0);
	double BBnomer = 0;
	double BBreal = 0;
	double BBYrl = 0;
	for (int z = 0; z < setkaBB; z++) {
		double BBZnomerk = setkaBBkvadr * z;
		for (double i = 0; i < setkaBB; i++) {
			BBnomer = BBZnomerk + i * setkaBB;
			double BBXrl = -BBXmax + i * dBBX;
			for (double j = 0; j < setkaBB; j++) {
				BBYrl = -BBYmax + j * dBBY;
				k = k + fka[BBnomer + j] * BBYrl * BBYrl * dBBY * dBBX * dBBZ;
			}
		}
	}
	return k;
}
complex <double> ENERGY_X(const std::vector<complex<double> >& fka) {
	complex<double> k(0, 0);
	double BBnomer = 0;
	double BBreal = 0;
	double BBYrl = 0;
	for (int z = 0; z < setkaBB; z++) {
		double BBZnomerk = setkaBBkvadr * z;
		for (double i = 0; i < setkaBB; i++) {
			BBnomer = BBZnomerk + i * setkaBB;
			double BBXrl = -BBXmax + i * dBBX;
			for (double j = 0; j < setkaBB; j++) {
				BBYrl = -BBYmax + j * dBBY;
				k = k + fka[BBnomer + j] * BBXrl * BBXrl * dBBY * dBBX * dBBZ;
			}
		}
	}
	return k;
}
complex <double> ENERGY_Z(const std::vector<complex<double> >& fka) {
	complex<double> k(0, 0);
	double BBnomer = 0;
	double BBreal = 0;
	double BBYrl = 0;
	for (int z = 0; z < setkaBB; z++) {
		double BBZnomerk = setkaBBkvadr * z;
		double BBZrl = -BBZmax + z * dBBZ;


		for (double i = 0; i < setkaBB; i++) {
			BBnomer = BBZnomerk + i * setkaBB;
			double BBXrl = -BBXmax + i * dBBX;
			for (double j = 0; j < setkaBB; j++) {
				BBYrl = -BBYmax + j * dBBY;
				k = k + fka[BBnomer + j] * BBZrl * BBZrl * dBBY * dBBX * dBBZ;
			}
		}
	}
	return k;

}
complex <double> ENERGY_PERP(const std::vector<complex<double> >& fka) {
	complex<double> k(0, 0);
	double BBnomer = 0;
	double BBreal = 0;
	double BBYrl = 0;
	for (int z = 0; z < setkaBB; z++) {
		double BBZnomerk = setkaBBkvadr * z;
		double BBZrl = -BBZmax + z * dBBZ;


		for (double i = 0; i < setkaBB; i++) {
			BBnomer = BBZnomerk + i * setkaBB;
			double BBXrl = -BBXmax + i * dBBX;
			for (double j = 0; j < setkaBB; j++) {

				k = k + fka[BBnomer + j] * (BBZrl * BBZrl + BBXrl * BBXrl) * dBBY * dBBX * dBBZ;
			}
		}
	}
	return k;
}
complex <double> DENSITY(const std::vector<complex<double> >& fka) {
	complex<double> k(0, 0);
	double BBnomer = 0;
	double BBreal = 0;
	double BBYrl = 0;
	for (int z = 0; z < setkaBB; z++) {
		double BBZnomerk = setkaBBkvadr * z;

		for (double i = 0; i < setkaBB; i++) {
			BBnomer = BBZnomerk + i * setkaBB;
			double BBXrl = -BBXmax + i * dBBX;
			for (double j = 0; j < setkaBB; j++) {
				BBYrl = -BBYmax + j * dBBY;
				k = k + fka[BBnomer + j] * dBBY * dBBX * dBBZ;
			}
		}
	}
	return k;
}

