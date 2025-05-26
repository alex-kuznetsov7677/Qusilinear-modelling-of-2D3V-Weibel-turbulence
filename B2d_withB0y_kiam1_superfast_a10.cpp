#include <cmath>
#include <complex>
#include <vector>
#include <cstdlib>
#include <random>
#include <fstream>
#include <algorithm> 

#include <mpi.h>


double dk = 0.2;
double MAGNIT =2.;//1./1.269;//2./sqrt(6.);//4./sqrt(6.)/1.39;///6./sqrt(6.)/1.5;//5./sqrt(6.)/1.25;//2./sqrt(6.)/1.29;//1./sqrt(6.)*sqrt(10./(7.));// 1./sqrt(6.)*sqrt(10./(1.5*7.));
//double Beta_prod1 = 0.1;
double A = -10;
double Beta_perp = 0.1;

double Pi = 3.14159265;
double Tmax = 650.;				// время работы программы
double dt = 0.04;
int setkaBB = 80;
int setkaBBkvadr = setkaBB * setkaBB;
int setkaBBminus = setkaBB - 1;
int setkaBBkub = setkaBB * setkaBB * setkaBB;

double dBBX = Beta_perp*7./setkaBB;
double dBBY = dBBX * sqrt(1. - A);
double dBBZ = dBBX;
double dBBY2 = dBBY * 2;
double dBBYkvadr = dBBY * dBBY;
double dBBXkvadr = dBBX * dBBX;
double dBBXdBBY = 4 * dBBY * dBBX;
double BBXmax = Beta_perp * 3.5 - (dBBX / 2.);
double BBYmax = BBXmax * sqrt(1. - A);
double BBZmax = BBXmax;
double dBBX2 = dBBX * 2;
double dBBZ2 = dBBZ * 2;

using namespace std;
double c = 299792458;
double Ngarmonikconst = 160;
complex<double> nul(0., 0.);
complex<double> I(0., 1.);
complex<double> b0(0., 0.);
complex<double> Vb0(0., 0.);
complex<double> fk0(0., 0.);
double BetaperpsqrtA = Beta_perp * sqrt(1 - A);
double C1 = 2. * 2. / ((1. - A) * pow(Beta_perp, 3.) * pow(2. * Pi, 3. / 2.));
double C2 = Beta_perp * (1. - A) * Beta_perp / (8. * Pi); // отношение квадратов циклотронной и плазменной частоты, деленное на 4
double C3 = Beta_perp * sqrt(1. - A) / (2 * pow(2., 3. / 2.) * sqrt(Pi));
// лучевая сетка. Сохраняет аксиальную симметрию, но нарушает спектральную плотность гармоник. 
/*void  Postanovka1(std::vector<double>& Kxvector,std::vector<double>& Kzvector){
	vector <double> Kmod;
	vector <double> Kphi;
	double Ngarmonik_mod = 14;
	double Ngarmonik_phi = 6;
	for (double i = 0; i <Ngarmonik_mod ; i++) {
		Kmod.push_back(0.3+i*0.15);
	}
	for (double i = 0; i < Ngarmonik_phi; i++) {
		Kphi.push_back(i*Pi/6.);
	}
	for (double i = 0; i <Ngarmonik_mod ; i++) {
		for (double j = 0; j < Ngarmonik_phi; j++) {
			Kxvector[j*Ngarmonik_mod+i]=Kmod[i]*cos(Kphi[j]);
			Kzvector[j*Ngarmonik_mod+i]=Kmod[i]*sin(Kphi[j]);
		}
	}



}*/

void  Postanovka1(std::vector<double>& Kxvector, std::vector<double>& Kzvector) {
	vector <double> Kmod;
	vector <double> Kphi;
	double Ngarmonik_mod = 32;
	double Ngarmonik_phi = 5;
	//double Ngarmonik_phi = 10;
	for (double i = 0; i < Ngarmonik_mod; i++) {
		Kmod.push_back(0.2 + i * 0.045);
	}
	for (double i = 0; i < Ngarmonik_phi; i++) {
		Kphi.push_back(i*Pi/5);
	}
	for (double i = 0; i < Ngarmonik_mod; i++) {
		for (double j = 0; j < Ngarmonik_phi; j++) {
			Kxvector[j * Ngarmonik_mod + i] = Kmod[i] * cos(Kphi[j]);
			Kzvector[j * Ngarmonik_mod + i] = Kmod[i] * sin(Kphi[j]);
		}
	}



}






complex<double>DifY(long int& fk, std::vector<complex<double> >& fka) {
	complex<double> k;
	//long int fk = N;
	if (fk - setkaBB != 0 && (fk - setkaBBminus) % setkaBB != 0) {
		return k = (fka[fk + 1] - fka[fk - 1]) / dBBY2;


	}
	else {
		return k = 0;

	}

}
complex<double>DifZ(long int& fk, std::vector<complex<double> >& fka, int& p) {
	complex<double> k;

	long int N = fk - fk % setkaBBkvadr;
	if (N % setkaBBkub != 0 && (N - setkaBBminus * setkaBBkvadr) % setkaBBkub != 0) {
		return k = (fka[fk + setkaBBkvadr] - fka[fk - setkaBBkvadr]) / dBBZ2;


	}
	else {
		return k = 0;

	}

}

complex<double>DifX(long int& N, std::vector<complex<double> >& fka) {
	complex<double> k;
	long int fk = (N - N % setkaBB);
	long int Nz = (fk - fk % setkaBBkvadr);

	if (fk - Nz != 0 && (fk - setkaBBminus * setkaBB) != Nz) {
		return k = (fka[N + setkaBB] - fka[N - setkaBB]) / dBBX2;


	}
	else {
		return k = 0;
	}

}


std::vector<complex<double> >& Oblast2(int oblast, int size, int Ngarmonik, double dt1,  std::vector<complex<double> >& BEXTvectorx,std::vector<complex<double> >& BEXTvectorz, std::vector<double>& Kxvector,std::vector<double>& Kzvector, std::vector<complex<double> >& IEy, std::vector<complex<double> >& bx, std::vector<complex<double> >& bz, std::vector<complex<double> >& fk,
	std::vector<complex<double> >& fk2cel, std::vector<complex<double> >& f0k2cel, double NUU, std::vector<complex<double> >& parallel_mas1) {
	//#pragma omp parallel for
	for (int p = 0; p < Ngarmonik / size; p++)
	{
		int zero = 0;
		int preal = p + oblast * Ngarmonik / size;
		complex<double> Bx = bx[preal];
		complex<double> Bz = bz[preal];
		complex<double> IIEy = I * IEy[preal];
		for (int z = 0; z < setkaBB; z++) {
			int BBZnomerk = setkaBBkvadr * z;
			double BBZrl = -BBZmax + z * dBBZ;
			for (int j = 0; j < setkaBB; j++) {
				int BBXnomer = j * setkaBB;
				double BBXrl = -BBXmax + j * dBBX;
				for (int k = 0; k < setkaBB; k++) {
					int BBXnomerk = p * setkaBBkub + BBZnomerk + BBXnomer + k;
					long int BBXnomerkF0 = BBXnomer + BBZnomerk + k;
					long int BBXnomerkreal = preal * setkaBBkub + BBZnomerk + BBXnomer + k;
					double BBYrl = -BBYmax + k * dBBY; //настоящий полярный угол
					

					complex<double> df0x = DifX(BBXnomerkF0, f0k2cel);
					complex<double> df0y = DifY(BBXnomerkF0, f0k2cel);
					complex<double> df0z = DifZ(BBXnomerkF0, f0k2cel, zero);
					complex<double> df0g = (BBZrl * df0y - BBYrl * df0z);
					complex<double> df0a = (BBYrl * df0x - BBXrl * df0y);

					complex<double> dfa = (BBXrl * DifZ(BBXnomerkreal, fk2cel, preal) - BBZrl * DifX(BBXnomerkreal, fk2cel));


					parallel_mas1[BBXnomerk] = fk[BBXnomerkreal] -dt1*NUU*fk2cel[BBXnomerkreal]- dt1 * (I * BBXrl * Kxvector[preal] * fk2cel[BBXnomerkreal]+I * BBZrl * Kzvector[preal] * fk2cel[BBXnomerkreal]  + C3 * (2. * MAGNIT * dfa - 2. * IIEy * df0y + 2. * (BEXTvectorz[preal]+Bz) * df0a+ 2. * (BEXTvectorx[preal]+Bx) * df0g));		// основной цикл

				}

			}
		}

	}
}

void  Oblast3(int oblast, int size, int Ngarmonik, double dt1,  std::vector<complex<double> >& BEXTvectorx,std::vector<complex<double> >& BEXTvectorz, std::vector<double>& Kxvector,std::vector<double>& Kzvector,  std::vector<complex<double> >& IEy1, std::vector<complex<double> >& bx1, std::vector<complex<double> >& bz1, std::vector<complex<double> >& fk2,
	std::vector<complex<double> >& f0k2, std::vector<complex<double> >& f0k2cel, std::vector<complex<double> >& f0kcel, std::vector<complex<double> >& f0k_MAXW, double NUU, std::vector<complex<double> >& parallel_mas2) {
	complex<double> dfx;
	complex<double> dfy;
	complex<double> dfz;
	complex<double> dfa;
	
	complex<double> dfg;
	int BBZnomerk;
	long int BBXnomerk;
	int BBZnomerkMAS;
	for (int p = 0; p < Ngarmonik; p++)
	{
		int zero = 0;

		complex<double> Bx = bx1[p];
		complex<double> Bz = bz1[p];
		complex<double> IIEy = I * IEy1[p];

		//#pragma omp parallel for
		for (int z = oblast * setkaBB / size; z < (oblast + 1) * setkaBB / size; z++) {
			BBZnomerkMAS = setkaBBkvadr * (z - oblast * setkaBB / size);
			BBZnomerk = setkaBBkvadr * z;
			double BBZrl = -BBZmax + z * dBBZ;


			for (int j = 0; j < setkaBB; j++) {
				long int BBXnomer = j * setkaBB;
				double BBXrl = -BBXmax + j * dBBX;
				for (int k = 0; k < setkaBB; k++) {
					BBXnomerk = p * setkaBBkub + BBZnomerk + BBXnomer + k;
					long int BBXnomerkF0 = BBZnomerk + BBXnomer + k;
					long int BBXnomerkF0MAS = BBZnomerkMAS + BBXnomer + k;

					double BBYrl = -BBYmax + k * dBBY;
	
					dfx = DifX(BBXnomerk, fk2);
					dfy = DifY(BBXnomerk, fk2);
					dfz = DifZ(BBXnomerk, fk2, p);

					dfg = (BBZrl * dfy - BBYrl * dfz);
					dfa = (BBYrl * dfx - BBXrl * dfy);

					complex<double> df0a = (BBXrl * DifZ(BBXnomerkF0, f0k2, zero) - BBZrl * DifX(BBXnomerkF0, f0k2));


					if (p > 0.5) {
						complex<double> yh1 =  Bz * conj(dfa)-IIEy * conj(dfy) + Bx * conj(dfg);
						parallel_mas2[BBXnomerkF0MAS] = parallel_mas2[BBXnomerkF0MAS] - dt1 * 2 * C3 * real(yh1);


					}
					else {

						complex<double> yh1 = MAGNIT * df0a - IIEy  * conj(dfy) + Bz * conj(dfa)+ Bx * conj(dfg);
						parallel_mas2[BBXnomerkF0MAS] = f0kcel[BBXnomerkF0] - dt1 * 2 * C3 * real(yh1)-dt1 *NUU*(f0k2[BBXnomerkF0]-f0k_MAXW[BBXnomerkF0]);
					}
				}

			}
		}
	}


}

complex <double> Vysrednekvadr(std::vector<complex<double> >& fka) {

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
complex <double> Vxsrednekvadr(std::vector<complex<double> >& fka) {
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
complex <double> Vzsrednekvadr(std::vector<complex<double> >& fka) {
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
complex <double> Vperpsrednekvadr(std::vector<complex<double> >& fka) {
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

double C4 = 2. * sqrt(2.) * Pi / (sqrt(Pi) * Beta_perp * sqrt(1. - A));
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


complex <double> Konc(std::vector<complex<double> >& fka) {
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


int main(int argc, char** argv)
{

	double endtime;
	//double itime = omp_get_wtime( );
	double starttime;

	//omp_set_num_threads(10);


	double Bext = 0.;
	int df = 6230023;
	if ((df - df % 100) % 10000 == 0) {

	}
	int Ngarmonik = Ngarmonikconst;

	vector <double> Kxvector;
	vector <double> Kzvector;

	vector <double> Tvector;
	vector <complex<double> > BEXTvectorz;
	vector <complex<double> > BEXTvectorx;

	for (double i = 0; i < Ngarmonik; i++) {
		Kxvector.push_back(0.);
		Kzvector.push_back(0.);

	}
	Postanovka1(Kxvector,Kzvector);


	vector <double> rand_nu;
	for (double i = 0; i < Ngarmonik; i++) {
		double t = i;
		double phi(sqrt(2) * i / 5.);
		complex< double > z(exp(0));
		BEXTvectorz.push_back(z * 0.001*Kxvector[i]/sqrt(Kxvector[i]*Kxvector[i]+Kzvector[i]*Kzvector[i]));
		BEXTvectorx.push_back(-z * 0.001*Kzvector[i]/sqrt(Kxvector[i]*Kxvector[i]+Kzvector[i]*Kzvector[i]));

	}

	for (double i = 0; i < Ngarmonik; i++) {
		Tvector.push_back(0);

	}
	double T0 = Tvector[0];
	for (double i = 0; i < Ngarmonik; i++) {
		Tvector[i] = Tvector[i] - T0;

	}

	
	vector <complex<double> > parallel_mas1;
	vector <complex<double> > parallel_mas2;
	vector <complex<double> > bz;
	vector <complex<double> > bz1; // половинное мп	
	vector <complex<double> > bzlast;
	vector <complex<double> > bz1last; // половинное мп	
	vector <complex<double> > bx;
	vector <complex<double> > bx1; // половинное мп	
	vector <complex<double> > bxlast;
	vector <complex<double> > bx1last; // половинное мп
	

	vector <complex<double> > IEy;
	vector <complex<double> > IEy1; //половинная скорость мп
	vector <complex<double> > IEylast;
	vector <complex<double> > IEy1last; //половинная скорость мп

	vector <complex<double> > fk;
	
	vector <complex<double> > f0_MAXW;
	vector <complex<double> > f0k;
	vector <complex<double> > f0k2;
	vector <complex<double> > f0kcel;
	vector <complex<double> > f0k2cel;
	vector <complex<double> > fkcel;
	vector <double> f0_proect;
	for (int i = 0; i < setkaBBkvadr; i++) {
		f0_proect.push_back((0., 0.));
	}


	vector <double> Areal;
	ofstream CMFt100000("DESCRIPTIONRe[B]__168kiam.txt");

	CMFt100000 <<"2DWEIBEL__withB0y;A=10; dt=0.04(0.04x5),Beta_perp=0.1;160garmoniki((31x5)POSTANOVKA_1_Kmod(0.2:0.045:0.2+31.*0.045)_Kphi(i*pi/5))" << '\n';
	CMFt100000 << "Bext=0;"<< '\n'<< "1 / (Beta_perp * pow(Pi, 1.5) * Beta_perp * BetaperpsqrtA) * exp(-pow(BBXrl / Beta_perp, 2.) - pow(BBYrl / BetaperpsqrtA, 2.) - pow(BBZrl / Beta_perp, 2.));";
	CMFt100000 << "NU=double NU_effective = 3 / 64. * BetaperpsqrtA * BetaperpsqrtA * B_MOD * B_MOD / KB_MOD/sqrt_VV0;;";
	


	ofstream CMFt2000("abs[Bz]__168kiamMATLAB.txt");
	ofstream CMFt2001("abs[Bx]__168kiamMATLAB.txt");
	ofstream CMFt("ImFUNC000__168kiam.txt");
	ofstream CMFt0("ReFUNC000_168kiam.txt");
	ofstream CMFt1301("Re,Wy0]__168kiam.txt");
	ofstream CMFt1303("Re,Wx0]__168kiam.txt");
	ofstream CMFt1302("Im,Areal]__168kiam.txt");
	ofstream CMFt30("Re[BkSREDN]__168kiam.txt");
	ofstream CMFt31("Re[AREALWyWx]_168kiam.txt");
	ofstream CMFt32("Re[Time]__168kiam.txt");


	//необходимые константы



	//первую итерацию удобно вынести за цикл в силу ее упрощенности в сравнении со всем последующим говном


	for (int i = 0; i < Ngarmonik + 10; i++) {
		bx.push_back((0., 0.));
		bx1.push_back((0., 0.));
		bxlast.push_back((0., 0.));
		bx1last.push_back((0., 0.));

		bz.push_back((0., 0.));
		bz1.push_back((0., 0.));
		bzlast.push_back((0., 0.));
		bz1last.push_back((0., 0.));

		IEy.push_back((0., 0.));
		IEy1.push_back((0., 0.));
		IEylast.push_back((0., 0.));
		IEy1last.push_back((0., 0.));
	}

	for (int i = 0; i < setkaBBkub * Ngarmonik; i++) {

		fk.push_back((0., 0.));
		fkcel.push_back((0., 0.));
		

	
	}
	for (int i = 0; i < setkaBBkub; i++) {
		f0k.push_back((0., 0.));
		f0k2.push_back((0., 0.));
		f0kcel.push_back((0., 0.));
		f0k2cel.push_back((0., 0.));
		f0_MAXW.push_back((0., 0.));

	}

	double Beta_potok=0.6;
	double Beta_2=0.06;
	double normirovka_potokov=1;
	double Beta_2potok=0.06*normirovka_potokov;

	for (int z = 0; z < setkaBB; z++) {
		double BBZnomer = z * setkaBBkvadr;
		double BBZrl = -BBZmax + z * dBBZ;
		for (double j = 0; j < setkaBB; j++) {
			double BBXnomer = j * setkaBB;
			double BBXrl = -BBXmax + j * dBBX;
			for (double k = 0; k < setkaBB; k++) {
				double BBYrl = -BBYmax + k * dBBY;
				//сосисочный бимаксвелл
				f0k2[BBZnomer + BBXnomer + k] = 1 / (Beta_perp * pow(Pi, 1.5) * Beta_perp * BetaperpsqrtA) * exp(-pow(BBXrl / Beta_perp, 2.) - pow(BBYrl / BetaperpsqrtA, 2.) - pow(BBZrl / Beta_perp, 2.));
				f0k2cel[BBZnomer + BBXnomer + k] = 1 / (Beta_perp * pow(Pi, 1.5) * Beta_perp * BetaperpsqrtA) * exp(-pow(BBXrl / Beta_perp, 2.) - pow(BBYrl / BetaperpsqrtA, 2.) - pow(BBZrl / Beta_perp, 2.));

				f0_MAXW[BBZnomer + BBXnomer + k] = 1 / (Beta_perp * pow(Pi, 1.5) * Beta_perp * Beta_perp) * exp(-pow(BBXrl / Beta_perp, 2.) - pow(BBYrl / Beta_perp, 2.) - pow(BBZrl / Beta_perp, 2.));

				//ненормированное тримаксвелловское распределение
				//f0k2[BBZnomer + BBXnomer + k] = 2.* exp(-pow(BBXrl / Beta_2, 2.) - pow(BBYrl / Beta_2, 2.) - pow(BBZrl / Beta_2, 2.))+exp(-pow(BBXrl / Beta_2potok, 2.) - pow((BBYrl-Beta_potok) / Beta_2potok, 2.) - pow(BBZrl / Beta_2potok, 2.))/normirovka_potokov+exp(-pow(BBXrl / Beta_2potok, 2.) - pow((BBYrl+Beta_potok) /Beta_2potok, 2.) - pow(BBZrl / Beta_2potok, 2.))/normirovka_potokov;
				//f0k2cel[BBZnomer + BBXnomer + k] =2.* exp(-pow(BBXrl / Beta_2, 2.) - pow(BBYrl / Beta_2, 2.) - pow(BBZrl / Beta_2, 2.))+exp(-pow(BBXrl / Beta_2potok, 2.) - pow((BBYrl-Beta_potok) / Beta_2potok, 2.) - pow(BBZrl / Beta_2potok, 2.))/normirovka_potokov+exp(-pow(BBXrl / Beta_2potok, 2.) - pow((BBYrl+Beta_potok) /Beta_2potok, 2.) - pow(BBZrl / Beta_2potok, 2.))/normirovka_potokov;
			}
		}
	}

//нужно для трехмаксвелловского распределения (или для любого другого с мутной нормировкой)
	double KONNNCC=real(Konc(f0k2));

	for (int z = 0; z < setkaBB; z++) {
		double BBZnomer = z * setkaBBkvadr;
		for (double j = 0; j < setkaBB; j++) {
			double BBXnomer = j * setkaBB;
			for (double k = 0; k < setkaBB; k++) {
				f0k2[BBZnomer + BBXnomer + k] =f0k2[BBZnomer + BBXnomer + k] /KONNNCC;
				f0k2cel[BBZnomer + BBXnomer + k] =f0k2cel[BBZnomer + BBXnomer + k] /KONNNCC;

			}
		}
	}

	int rank, size;
	MPI_Init(&argc, &argv);
	starttime = MPI_Wtime();
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	for (int i = 0; i < setkaBBkub * Ngarmonik / size; i++) {
		parallel_mas1.push_back((0., 0.));
	}
	for (int i = 0; i < setkaBBkub / size; i++) {
		parallel_mas2.push_back((0., 0.));
	}

	
	for (int i = 2; i < Tmax / dt; i++) {
		if (i == 50) {
			for (double ks = 0; ks < Ngarmonik; ks++) {
				BEXTvectorz[ks] = 0.;
				BEXTvectorx[ks] = 0.;
			}
		}


		
		f0kcel = f0k2cel;
	

		double sqrt_VVy0 = sqrt(real(Vysrednekvadr(f0k2cel))+real(Vxsrednekvadr(f0k2cel)));

		double B_MOD = 0.; double KB_MOD = 0.; double K2B_MOD = 0.;

		for (int p = 0; p < Ngarmonik; p++) {
				B_MOD = B_MOD + (pow(abs(bxlast[p] + 0.0000000001), 2)+pow(abs(bzlast[p] + 0.0000000001), 2));
				KB_MOD = KB_MOD + sqrt(Kxvector[p] *Kxvector[p]+Kzvector[p] *Kzvector[p])*(pow(abs(bxlast[p] + 0.0000000001), 2)+pow(abs(bzlast[p] + 0.0000000001), 2));
				K2B_MOD = K2B_MOD + (Kxvector[p] *Kxvector[p]+Kzvector[p] *Kzvector[p])*(pow(abs(bxlast[p] + 0.0000000001), 2)+pow(abs(bzlast[p] + 0.0000000001), 2));
		}
		


		double NU_effective = 3 / 64. * BetaperpsqrtA * BetaperpsqrtA * B_MOD* B_MOD/ KB_MOD/sqrt_VVy0;




		Oblast3(rank, size, Ngarmonik, dt, BEXTvectorx,BEXTvectorz,  Kxvector,Kzvector, IEy1last, bx1last, bz1last,fk, f0k2, f0k2cel, f0kcel,f0_MAXW,NU_effective, parallel_mas2);

		MPI_Allgather(parallel_mas2.data(), setkaBBkub / size, MPI_DOUBLE_COMPLEX, f0k2cel.data(), setkaBBkub / size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);

		Oblast2(rank, size, Ngarmonik, dt, BEXTvectorx,BEXTvectorz,  Kxvector,Kzvector, IEylast, bxlast, bzlast,fk, fkcel, f0k2cel,NU_effective, parallel_mas1);

		MPI_Allgather(parallel_mas1.data(), Ngarmonik * setkaBBkub / size, MPI_DOUBLE_COMPLEX, fk.data(), Ngarmonik * setkaBBkub / size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);


		if (rank == 0) {
			for (int p = 0; p < Ngarmonik; p++) {
				bz1[p] = bz1last[p] -Kxvector[p]*IEylast[p] * dt;
				bx1[p] = bx1last[p] +Kzvector[p]*IEylast[p] * dt;
				complex <double> Integ = SpeedIntegral(fk, p);
				IEy[p] = IEylast[p] + dt * (Kxvector[p] * bz1[p] -Kzvector[p] * bx1[p]- I * Integ);
				IEy1[p] = (IEy[p] + IEylast[p]) / 2.;
				bz[p] = bzlast[p] - Kxvector[p] * IEy1[p] * dt;
				bx[p] = bxlast[p] + Kzvector[p] * IEy1[p] * dt;
			}
			int pk = i;
			if (pk % 5 == 0) {
				double VVx0 = real(Vxsrednekvadr(f0k2cel));
				double VVy0 = real(Vysrednekvadr(f0k2cel));
				double VVperp0 = real(Vperpsrednekvadr(f0k2cel));
				double ArealWydelWx = (VVy0) / (VVx0)-1;
				double Areal = 2. * (VVy0) / (VVperp0)-1;
				double Bmodkvadr = 0.;
				for (double it = 0; it < Ngarmonikconst; it++) {
					Bmodkvadr = Bmodkvadr + pow(abs(bx[it]), 2)+ pow(abs(bz[it]), 2);

				}

				double Bmod = sqrt(Bmodkvadr);


				for (double p = 0; p < Ngarmonikconst; p++) {
					double ip = p * Tmax + i;
					CMFt2000 << real(bz[p]) << " ";
					CMFt2001 << real(bx[p]) << " ";

				}

				CMFt2000 << ";";
				CMFt2001 << ";";
				CMFt30 << Bmod << ",";
	
				CMFt1301 << VVy0 << ",";
				CMFt1301 << real(Vzsrednekvadr(f0k2cel)) << ",";
				CMFt1303 << VVx0 << ",";
				CMFt1302 << Areal << ",";
				CMFt31 << ArealWydelWx << ",";
				//CMFt32 << ArealSrednKVADR << ",";
			}
			if (pk - 10 == 0) {


				
				CMFt0 << '\n' << pk << '\n';
				for (double z = 0; z < setkaBB; z++) {

					for (double j = 0; j < setkaBB; j++) {
						double BBXnomer = j * setkaBB;
						for (double k = 0; k < setkaBB; k++) {
							double BBXnomerk = z * setkaBBkvadr + BBXnomer + k;
							double BBXnomerk2 = z * setkaBB + k;
 							f0_proect[BBXnomerk2]=f0_proect[BBXnomerk2]+real(f0k2cel[BBXnomerk])*dBBX;
						}


					}
				}

				for (double z = 0; z < setkaBB; z++) {
					for (double j = 0; j < setkaBB; j++) {											
						double BBXnomerk2 = z * setkaBB + j;							
 						CMFt0 << f0_proect[BBXnomerk2] << " ";
						f0_proect[BBXnomerk2]=0;
					}
					CMFt0 << ";";
				}
				CMFt0 << '\n' << '\n' << '\n';


			}

			if (pk - 700 == 0) {


				CMFt0 << '\n' << pk << '\n';

				for (double z = 0; z < setkaBB; z++) {

					for (double j = 0; j < setkaBB; j++) {
						double BBXnomer = j * setkaBB;
						for (double k = 0; k < setkaBB; k++) {
							double BBXnomerk = z * setkaBBkvadr + BBXnomer + k;
							double BBXnomerk2 = z * setkaBB + k;
 							f0_proect[BBXnomerk2]=f0_proect[BBXnomerk2]+real(f0k2cel[BBXnomerk])*dBBX;
						}


					}
				}

				for (double z = 0; z < setkaBB; z++) {
					for (double j = 0; j < setkaBB; j++) {											
						double BBXnomerk2 = z * setkaBB + j;							
 						CMFt0 << f0_proect[BBXnomerk2] << " ";
						f0_proect[BBXnomerk2]=0;
					}
					CMFt0 << ";";
				}
				CMFt0 << '\n' << '\n' << '\n';


			}
			if (pk - 850 == 0) {



				CMFt0 << '\n' << pk << '\n';
				for (double z = 0; z < setkaBB; z++) {

					for (double j = 0; j < setkaBB; j++) {
						double BBXnomer = j * setkaBB;
						for (double k = 0; k < setkaBB; k++) {
							double BBXnomerk = z * setkaBBkvadr + BBXnomer + k;
							double BBXnomerk2 = z * setkaBB + k;
 							f0_proect[BBXnomerk2]=f0_proect[BBXnomerk2]+real(f0k2cel[BBXnomerk])*dBBX;
						}


					}
				}

				for (double z = 0; z < setkaBB; z++) {
					for (double j = 0; j < setkaBB; j++) {											
						double BBXnomerk2 = z * setkaBB + j;							
 						CMFt0 << f0_proect[BBXnomerk2] << " ";
						f0_proect[BBXnomerk2]=0;
					}
					CMFt0 << ";";
				}
				CMFt0 << '\n' << '\n' << '\n';


			}


				if (pk - 1000 == 0) {



				CMFt0 << '\n' << pk << '\n';
				for (double z = 0; z < setkaBB; z++) {

					for (double j = 0; j < setkaBB; j++) {
						double BBXnomer = j * setkaBB;
						for (double k = 0; k < setkaBB; k++) {
							double BBXnomerk = z * setkaBBkvadr + BBXnomer + k;
							double BBXnomerk2 = z * setkaBB + k;
 							f0_proect[BBXnomerk2]=f0_proect[BBXnomerk2]+real(f0k2cel[BBXnomerk])*dBBX;
						}


					}
				}

				for (double z = 0; z < setkaBB; z++) {
					for (double j = 0; j < setkaBB; j++) {											
						double BBXnomerk2 = z * setkaBB + j;							
 						CMFt0 << f0_proect[BBXnomerk2] << " ";
						f0_proect[BBXnomerk2]=0;
					}
					CMFt0 << ";";
				}
				CMFt0 << '\n' << '\n' << '\n';


			}
	if (pk - 1150 == 0) {



				CMFt0 << '\n' << pk << '\n';
				for (double z = 0; z < setkaBB; z++) {

					for (double j = 0; j < setkaBB; j++) {
						double BBXnomer = j * setkaBB;
						for (double k = 0; k < setkaBB; k++) {
							double BBXnomerk = z * setkaBBkvadr + BBXnomer + k;
							double BBXnomerk2 = z * setkaBB + k;
 							f0_proect[BBXnomerk2]=f0_proect[BBXnomerk2]+real(f0k2cel[BBXnomerk])*dBBX;
						}


					}
				}

				for (double z = 0; z < setkaBB; z++) {
					for (double j = 0; j < setkaBB; j++) {											
						double BBXnomerk2 = z * setkaBB + j;							
 						CMFt0 << f0_proect[BBXnomerk2] << " ";
						f0_proect[BBXnomerk2]=0;
					}
					CMFt0 << ";";
				}
				CMFt0 << '\n' << '\n' << '\n';


			}
	if (pk - 1300 == 0) {



				CMFt0 << '\n' << pk << '\n';
				for (double z = 0; z < setkaBB; z++) {

					for (double j = 0; j < setkaBB; j++) {
						double BBXnomer = j * setkaBB;
						for (double k = 0; k < setkaBB; k++) {
							double BBXnomerk = z * setkaBBkvadr + BBXnomer + k;
							double BBXnomerk2 = z * setkaBB + k;
 							f0_proect[BBXnomerk2]=f0_proect[BBXnomerk2]+real(f0k2cel[BBXnomerk])*dBBX;
						}


					}
				}

				for (double z = 0; z < setkaBB; z++) {
					for (double j = 0; j < setkaBB; j++) {											
						double BBXnomerk2 = z * setkaBB + j;							
 						CMFt0 << f0_proect[BBXnomerk2] << " ";
						f0_proect[BBXnomerk2]=0;
					}
					CMFt0 << ";";
				}
				CMFt0 << '\n' << '\n' << '\n';


			}

		
		}

		
		f0k = f0k2;

		


		
		MPI_Bcast(bx1.data(), Ngarmonik, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		MPI_Bcast(IEy1.data(), Ngarmonik, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		MPI_Bcast(bx.data(), Ngarmonik, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		MPI_Bcast(IEy.data(), Ngarmonik, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		MPI_Bcast(bz.data(), Ngarmonik, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		MPI_Bcast(bz1.data(), Ngarmonik, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

		B_MOD = 0.; KB_MOD = 0.; K2B_MOD = 0.;

		for (int p = 0; p < Ngarmonik; p++) {
				
				B_MOD = B_MOD + (pow(abs(bx1[p] + 0.0000000001), 2)+pow(abs(bz1[p] + 0.0000000001), 2));
				KB_MOD = KB_MOD + sqrt(Kxvector[p] *Kxvector[p]+Kzvector[p] *Kzvector[p])*(pow(abs(bx1[p] + 0.0000000001), 2)+pow(abs(bz1[p] + 0.0000000001), 2));
				K2B_MOD = K2B_MOD + (Kxvector[p] *Kxvector[p]+Kzvector[p] *Kzvector[p])*(pow(abs(bx1[p] + 0.0000000001), 2)+pow(abs(bz1[p] + 0.0000000001), 2));
		}
		sqrt_VVy0 = sqrt(real(Vysrednekvadr(f0k2))+real(Vxsrednekvadr(f0k2)));
		NU_effective =3 / 64. * BetaperpsqrtA * BetaperpsqrtA * B_MOD* B_MOD/ KB_MOD/sqrt_VVy0;

	
		Oblast3(rank, size, Ngarmonik, dt, BEXTvectorx,BEXTvectorz, Kxvector,Kzvector, IEylast, bxlast, bzlast, fkcel, f0k2cel, f0k2, f0k, f0_MAXW, NU_effective, parallel_mas2);

		MPI_Allgather(parallel_mas2.data(), setkaBBkub / size, MPI_DOUBLE_COMPLEX, f0k2.data(), setkaBBkub / size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);

		Oblast2(rank, size, Ngarmonik, dt, BEXTvectorx,BEXTvectorz,Kxvector,Kzvector, IEy1, bx1, bz1, fkcel, fk, f0k2,NU_effective, parallel_mas1);
		MPI_Allgather(parallel_mas1.data(), Ngarmonik * setkaBBkub / size, MPI_DOUBLE_COMPLEX, fkcel.data(), Ngarmonik * setkaBBkub / size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);




		bx1last = bx1;
		IEy1last = IEy1;
		bxlast = bx;
		IEylast = IEy;
		bz1last = bz1;
		bzlast = bz;	


		//MPI_Bcast(b1last.data(), Ngarmonik, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		//MPI_Bcast(Vb1last.data(), Ngarmonik, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		//MPI_Bcast(blast.data(), Ngarmonik, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		//MPI_Bcast(Vblast.data(), Ngarmonik, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

	}

	endtime = MPI_Wtime();
	CMFt32 << endtime - starttime;

	MPI_Finalize();

	CMFt100000.close();

	CMFt2000.close();


	CMFt2001.close();
	CMFt30.close();
	CMFt.close();
	CMFt0.close();
	CMFt31.close();



	CMFt1301.close();

	CMFt1303.close();
	CMFt1302.close();






	//ftime = omp_get_wtime( );
		//exec_time = ftime - itime;



	CMFt32.close();
}
