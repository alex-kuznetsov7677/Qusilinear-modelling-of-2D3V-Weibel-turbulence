#include <cmath>
#include <complex>
#include <vector>
#include <cstdlib>
#include <random>
#include <fstream>
#include <algorithm> 
#include <iostream>
#include <mpi.h>



#include "my_variables.h" 
#include "Preparation.h" 
#include "Fields_changing.h" 
#include "Distribution_function_changing.h" 
#include "Integrals.h" 
#include "ConfigReader.h"
#include "OutputWriter.h"
using namespace std;

auto config = ReadConfig("config.txt");


int Ngarmonik = 160;

double dk = 0.2;
double MAGNIT =2.;
double A = -10;
double Beta_perp = 0.1;
double BetaperpsqrtA = Beta_perp * sqrt(1 - A);
double Pi = 3.14159265;

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

double Tmax = 10.;				// время работы программы
double dt = 0.04;



complex<double> I(0., 1.);

double C1 = 2. * 2. / ((1. - A) * pow(Beta_perp, 3.) * pow(2. * Pi, 3. / 2.));
double C2 = Beta_perp * (1. - A) * Beta_perp / (8. * Pi); // отношение квадратов циклотронной и плазменной частоты, деленное на 4
double C3 = Beta_perp * sqrt(1. - A) / (2 * pow(2., 3. / 2.) * sqrt(Pi));
double C4 = 2. * sqrt(2.) * Pi / (sqrt(Pi) * Beta_perp * sqrt(1. - A));



int main(int argc, char** argv)
{

	double endtime;
	double starttime;


	double Bext = 0.;
	int df = 6230023;
	if ((df - df % 100) % 10000 == 0) {

	}

	vector <double> Kxvector;
	vector <double> Kzvector;

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


	vector <complex<double> > Bx;
	vector <complex<double> > Bz;
		
	vector <complex<double> > parallel_mas1;
	vector <complex<double> > parallel_mas2;
	vector <complex<double> > parallel_mas3;
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
	cout<< Beta_potok;
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
	double KONNNCC=real(DENSITY(f0k2));

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
	for (int i = 0; i < Ngarmonik; i++) {
		Bx.push_back((0., 0.));
		Bz.push_back((0., 0.));
	}
	int rank, size;
	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	OutputFiles out;
	if (rank==0){
		out.log_number_MPI_processes(size);
		InitializeOutputFiles(out,config);
	}
	for (int i = 0; i < setkaBBkub * Ngarmonik / size; i++) {
		parallel_mas1.push_back((0., 0.));
	}
	for (int i = 0; i < setkaBBkub; i++) {
		parallel_mas2.push_back((0., 0.));
		parallel_mas3.push_back((0., 0.));
	}
	if (rank == 0) {
		starttime = MPI_Wtime();
	}

	for (int i = 2; i < Tmax / dt; i++) {
		if (i == 50) {
			for (double ks = 0; ks < Ngarmonik; ks++) {
				BEXTvectorz[ks] = 0.;
				BEXTvectorx[ks] = 0.;
			}
		}


		
		f0kcel = f0k2cel;
	
		double NU_effective =0;
		
		fill(parallel_mas2.begin(), parallel_mas2.end(), complex<double>(0., 0.));
		PERTURBATION_OF_UNIFORM_DISTRIBUTION(
			rank, size, 
			Ngarmonik, Kxvector,Kzvector, IEy1last, bx1last, bz1last,fk, parallel_mas2
		);
		MPI_Barrier(MPI_COMM_WORLD); 
		MPI_Reduce(parallel_mas2.data(), parallel_mas3.data(), setkaBBkub, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			transform(f0k2cel.begin(), f0k2cel.end(), parallel_mas3.begin(), f0k2cel.begin(), std::plus<complex<double> >());
		}

		PERTURBATION_OF_UNIFORM_DISTRIBUTION_MP(
			rank, size, 
			f0k2, f0k2cel,f0_MAXW,NU_effective
		);
		MPI_Barrier(MPI_COMM_WORLD); 
		MPI_Bcast(f0k2cel.data(), setkaBBkub, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);	
		PERTURBATION_OF_MODES_DISTRIBUTION(rank, size,Ngarmonik, BEXTvectorx,BEXTvectorz,  Kxvector,Kzvector, IEylast, bxlast, bzlast,fk, fkcel, f0k2cel,NU_effective, parallel_mas1);
		fk = parallel_mas1;	

		WEIBEL_FIELDS(rank, size, Ngarmonik,Kxvector,Kzvector, IEylast,bzlast, bxlast,IEy,bz, bx,IEy1,bz1, bx1,IEy1last,bz1last, bx1last,fk,fkcel);
		MPI_Gather(bx.data(), Ngarmonik / size, MPI_DOUBLE_COMPLEX, Bx.data(), Ngarmonik / size, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		MPI_Gather(bz.data(), Ngarmonik / size, MPI_DOUBLE_COMPLEX, Bz.data(), Ngarmonik / size, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

		int pk = i;
		if (rank == 0 && pk % 5 == 0) {
			double time = pk*dt;
            out.log_fields(time, Bx, Bz,Ngarmonik, f0k2cel, NU_effective);
		}
		if (rank == 0 && pk - 10 == 0) {
            out.log_distribution(f0k2cel,f0_proect);
		}

		f0k = f0k2;
		NU_effective=0;

		fill(parallel_mas2.begin(), parallel_mas2.end(), complex<double>(0., 0.));
		PERTURBATION_OF_UNIFORM_DISTRIBUTION(
			rank, size, 
			Ngarmonik, Kxvector,Kzvector, IEylast, bxlast, bzlast,fkcel,parallel_mas2
		);
		MPI_Barrier(MPI_COMM_WORLD); 
		MPI_Reduce(parallel_mas2.data(), parallel_mas3.data(), setkaBBkub, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			transform(f0k2.begin(), f0k2.end(), parallel_mas3.begin(), f0k2.begin(), std::plus<complex<double> >());
		}
		PERTURBATION_OF_UNIFORM_DISTRIBUTION_MP(
			rank, size, 
			f0k2cel, f0k2,f0_MAXW,NU_effective
		);
		MPI_Barrier(MPI_COMM_WORLD);  
		MPI_Bcast(f0k2.data(), setkaBBkub, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);	
		PERTURBATION_OF_MODES_DISTRIBUTION(rank, size,Ngarmonik, BEXTvectorx,BEXTvectorz,  Kxvector,Kzvector, IEy1, bx1, bz1,fkcel, fk, f0k2,NU_effective, parallel_mas1);
		fkcel = parallel_mas1;	

		bx1last = bx1;
		IEy1last = IEy1;
		bxlast = bx;
		IEylast = IEy;
		bz1last = bz1;
		bzlast = bz;	
	}

	if (rank == 0) {
		endtime = MPI_Wtime();
		double time_ex=endtime - starttime;
		out.log_execution_time(time_ex);
	}
	void close_all();

	MPI_Finalize();


}
