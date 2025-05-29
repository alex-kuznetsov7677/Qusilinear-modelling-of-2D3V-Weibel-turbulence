#include "my_variables.h"
#include "Distribution_function_changing.h"

#include <vector>
#include <cmath>
#include <complex> 

using std::complex;
using std::vector;

complex<double>DifY(int& fk, std::vector<complex<double> >& fka) {
	complex<double> k;
	if (fk - setkaBB != 0 && (fk - setkaBBminus) % setkaBB != 0) {
		return k = (fka[fk + 1] - fka[fk - 1]) / Vstepy2;
	}
	else {
		return k = 0;
	}
}
complex<double>DifZ(int& fk, std::vector<complex<double> >& fka, int& p) {
	complex<double> k;

	int N = fk - fk % setkaBBkvadr;
	if (N % setkaBBkub != 0 && (N - setkaBBminus * setkaBBkvadr) % setkaBBkub != 0) {
		return k = (fka[fk + setkaBBkvadr] - fka[fk - setkaBBkvadr]) / Vstepz2;
	}
	else {
		return k = 0;

	}

}

complex<double>DifX(int& N, std::vector<complex<double> >& fka) {
	complex<double> k;
	int fk = (N - N % setkaBB);
	int Nz = (fk - fk % setkaBBkvadr);

	if (fk - Nz != 0 && (fk - setkaBBminus * setkaBB) != Nz) {
		return k = (fka[N + setkaBB] - fka[N - setkaBB]) / Vstepx2;
	}
	else {
		return k = 0;
	}
}


void PERTURBATION_OF_MODES_DISTRIBUTION(int oblast, int size, int Ngarmonik,  std::vector<complex<double> >& BEXTvectorx,std::vector<complex<double> >& BEXTvectorz, std::vector<double>& Kxvector,std::vector<double>& Kzvector, std::vector<complex<double> >& IEy, std::vector<complex<double> >& bx, std::vector<complex<double> >& bz, std::vector<complex<double> >& fk,
	std::vector<complex<double> >& fk2cel, std::vector<complex<double> >& f0k2cel, double NUU, std::vector<complex<double> >& parallel_mas1) {

	for (int p = 0; p < Ngarmonik / size; p++)
	{
		int zero = 0;
		int preal = p + oblast * Ngarmonik / size;

		complex<double> Bx = bx[p];
		complex<double> Bz = bz[p];
		complex<double> IIEy = I * IEy[p];
		for (int z = 0; z < setkaBB; z++) {
			int BBZnomerk = setkaBBkvadr * z;
			double BBZrl = Vminz + z * Vstepz;
			for (int j = 0; j < setkaBB; j++) {
				int BBXnomer = j * setkaBB;
				double BBXrl = Vminx + j * Vstepx;
				for (int k = 0; k < setkaBB; k++) {
					int BBXnomerk = p * setkaBBkub + BBZnomerk + BBXnomer + k;
					int BBXnomerkF0 = BBXnomer + BBZnomerk + k;

					double BBYrl =Vminy + k * Vstepy; 
					

					complex<double> df0x = DifX(BBXnomerkF0, f0k2cel);
					complex<double> df0y = DifY(BBXnomerkF0, f0k2cel);
					complex<double> df0z = DifZ(BBXnomerkF0, f0k2cel, zero);
					complex<double> df0g = (BBZrl * df0y - BBYrl * df0z);
					complex<double> df0a = (BBYrl * df0x - BBXrl * df0y);

					complex<double> dfa = (BBXrl * DifZ(BBXnomerk, fk2cel, p) - BBZrl * DifX(BBXnomerk, fk2cel));


					parallel_mas1[BBXnomerk] = fk[BBXnomerk] -dt*NUU*fk2cel[BBXnomerk]- dt * (I * BBXrl * Kxvector[preal] * fk2cel[BBXnomerk]+I * BBZrl * Kzvector[preal] * fk2cel[BBXnomerk]  + C3 * (2. * MAGNIT * dfa - 2. * IIEy * df0y + 2. * (BEXTvectorz[preal]+Bz) * df0a+ 2. * (BEXTvectorx[preal]+Bx) * df0g));	

				}

			}
		}

	}
}

void PERTURBATION_OF_UNIFORM_DISTRIBUTION(int oblast, int size, int Ngarmonik, std::vector<double>& Kxvector,std::vector<double>& Kzvector,  std::vector<complex<double> >& IEy1, std::vector<complex<double> >& bx1,
	 std::vector<complex<double> >& bz1, std::vector<complex<double> >& fk2,  std::vector<complex<double> >& parallel_mas2) {

	for (int p = 0; p < Ngarmonik / size; p++) {
		int zero = 0;

		complex<double> Bx = bx1[p];
		complex<double> Bz = bz1[p];
		complex<double> IIEy = I * IEy1[p];

		for (int z = 0; z < setkaBB; z++) {

			int BBZnomerk = setkaBBkvadr * z;
			double BBZrl = Vminz + z * Vstepz;

			for (int j = 0; j < setkaBB; j++) {
				int BBXnomer = j * setkaBB;
				double BBXrl = Vminx + j * Vstepx;
				for (int k = 0; k < setkaBB; k++) {
					int BBXnomerkF0 = BBZnomerk + BBXnomer + k;
					int BBXnomerk = p * setkaBBkub + BBXnomerkF0;

					double BBYrl =Vminy + k * Vstepy; 	
					complex<double> dfx = DifX(BBXnomerk, fk2);
					complex<double> dfy = DifY(BBXnomerk, fk2);
					complex<double> dfz = DifZ(BBXnomerk, fk2, p);
					complex<double> dfg = (BBZrl * dfy - BBYrl * dfz);
					complex<double> dfa = (BBYrl * dfx - BBXrl * dfy);

					complex<double> yh1 =  Bz * conj(dfa)-IIEy * conj(dfy) + Bx * conj(dfg);
					parallel_mas2[BBXnomerkF0] = parallel_mas2[BBXnomerkF0] - dt * 2 * C3 * real(yh1);


				}

			}
		}
	}
}
void  PERTURBATION_OF_UNIFORM_DISTRIBUTION_MP(int oblast, int size, std::vector<complex<double> >& f0k2, std::vector<complex<double> >& f0k2cel, std::vector<complex<double> >& f0k_MAXW, double NUU) {
	int zero = 0;
	for (int z = 0; z < setkaBB; z++) {
		int BBZnomerk = setkaBBkvadr * z;
		double BBZrl = Vminz + z * Vstepz;
		for (int j = 0; j < setkaBB; j++) {
			int BBXnomer = j * setkaBB;
			double BBXrl = Vminx + j * Vstepx;
			for (int k = 0; k < setkaBB; k++) {
				int BBXnomerkF0 = BBZnomerk + BBXnomer + k;
				complex<double> df0a = (BBXrl * DifZ(BBXnomerkF0, f0k2, zero) - BBZrl * DifX(BBXnomerkF0, f0k2));	
				f0k2cel[BBXnomerkF0] = f0k2cel[BBXnomerkF0] - dt * 2 * C3 * real(MAGNIT * df0a)-dt *NUU*(f0k2[BBXnomerkF0]-f0k_MAXW[BBXnomerkF0]);
			}
		}
	}

}