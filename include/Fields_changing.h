#pragma once 
#include <vector>
#include <cmath>
#include <complex>

std::complex <double> SpeedIntegral(std::vector<std::complex<double> >& flka, double p11); 
void WEIBEL_FIELDS(int& oblast,
	int& size,
	int& Ngarmonik,
	std::vector<double>& Kxvector,std::vector<double>& Kzvector,
	std::vector<std::complex<double> >& IEylast, std::vector<std::complex<double> >& bzlast,std::vector<std::complex<double> >& bxlast,
	std::vector<std::complex<double> >& IEy, std::vector<std::complex<double> >& bz,std::vector<std::complex<double> >& bx,
	std::vector<std::complex<double> >& IEy1, std::vector<std::complex<double> >& bz1,std::vector<std::complex<double> >& bx1,
	std::vector<std::complex<double> >& IEy1last, std::vector<std::complex<double> >& bz1last, std::vector<std::complex<double> >& bx1last,
	std::vector<std::complex<double> >& fk, std::vector<std::complex<double> >& fkcel);