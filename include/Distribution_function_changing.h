#pragma once
#include <complex>
#include <vector>

std::complex<double>DifY(int& fk, std::vector<std::complex<double> >& fka);
std::complex<double>DifZ(int& fk, std::vector<std::complex<double> >& fka, int& p);
std::complex<double>DifX(int& N, std::vector<std::complex<double> >& fka);
void PERTURBATION_OF_MODES_DISTRIBUTION(int oblast, int size, int Ngarmonik,  std::vector<std::complex<double> >& BEXTvectorx,std::vector<std::complex<double> >& BEXTvectorz, std::vector<double>& Kxvector,std::vector<double>& Kzvector, std::vector<std::complex<double> >& IEy, std::vector<std::complex<double> >& bx, std::vector<std::complex<double> >& bz, std::vector<std::complex<double> >& fk,
	std::vector<std::complex<double> >& fk2cel, std::vector<std::complex<double> >& f0k2cel, double NUU, std::vector<std::complex<double> >& parallel_mas1);
void PERTURBATION_OF_UNIFORM_DISTRIBUTION(int oblast, int size, int Ngarmonik, std::vector<double>& Kxvector,std::vector<double>& Kzvector,  std::vector<std::complex<double> >& IEy1, std::vector<std::complex<double> >& bx1,
	 std::vector<std::complex<double> >& bz1, std::vector<std::complex<double> >& fk2,  std::vector<std::complex<double> >& parallel_mas2) ;
void  PERTURBATION_OF_UNIFORM_DISTRIBUTION_MP(int oblast, int size, std::vector<std::complex<double> >& f0k2, std::vector<std::complex<double> >& f0k2cel, std::vector<std::complex<double> >& f0k_MAXW, double NUU);