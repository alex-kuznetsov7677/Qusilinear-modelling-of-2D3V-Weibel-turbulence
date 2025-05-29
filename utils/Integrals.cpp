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
			for (double j = 0; j < setkaBB; j++) {
				BBYrl =Vminy + j * Vstepy;
				k = k + fka[BBnomer + j] * BBYrl * BBYrl * Vstepy * Vstepx * Vstepz;
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
			double BBXrl = Vminx + i * Vstepx;
			for (double j = 0; j < setkaBB; j++) {
				k = k + fka[BBnomer + j] * BBXrl * BBXrl * Vstepy * Vstepx * Vstepz;
			}
		}
	}
	return k;
}
complex <double> ENERGY_Z(const std::vector<complex<double> >& fka) {
	complex<double> k(0, 0);
	double BBnomer = 0;
	for (int z = 0; z < setkaBB; z++) {
		double BBZnomerk = setkaBBkvadr * z;
		double BBZrl =  Vminz + z * Vstepz;
		for (double i = 0; i < setkaBB; i++) {
			BBnomer = BBZnomerk + i * setkaBB;
			for (double j = 0; j < setkaBB; j++) {
				k = k + fka[BBnomer + j] * BBZrl * BBZrl * Vstepy * Vstepx * Vstepz;
			}
		}
	}
	return k;

}
complex <double> ENERGY_PERP(const std::vector<complex<double> >& fka) {
	complex<double> k(0, 0);
	double BBnomer = 0;
	for (int z = 0; z < setkaBB; z++) {
		double BBZnomerk = setkaBBkvadr * z;
		double BBZrl =  Vminz + z * Vstepz;


		for (double i = 0; i < setkaBB; i++) {
			BBnomer = BBZnomerk + i * setkaBB;
			double BBXrl = Vminx + i * Vstepx;
			for (double j = 0; j < setkaBB; j++) {
				k = k + fka[BBnomer + j] * (BBZrl * BBZrl + BBXrl * BBXrl) *Vstepy * Vstepx * Vstepz;
			}
		}
	}
	return k;
}
complex <double> DENSITY(const std::vector<complex<double> >& fka) {
	complex<double> k(0, 0);
	double BBnomer = 0;
	for (int z = 0; z < setkaBB; z++) {
		double BBZnomerk = setkaBBkvadr * z;
		for (double i = 0; i < setkaBB; i++) {
			BBnomer = BBZnomerk + i * setkaBB;
			for (double j = 0; j < setkaBB; j++) {
				k = k + fka[BBnomer + j] * Vstepy * Vstepx * Vstepz;
			}
		}
	}
	return k;
}

