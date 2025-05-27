#include "Preparation.h"

#include <vector>
#include <cmath>

#include "my_variables.h"
using std::vector;
// лучева€ сетка. —охран€ет аксиальную симметрию, но нарушает спектральную плотность гармоник. 
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
		Kphi.push_back(i*3.1415926535/5);
	}
	for (double i = 0; i < Ngarmonik_mod; i++) {
		for (double j = 0; j < Ngarmonik_phi; j++) {
			Kxvector[j * Ngarmonik_mod + i] = Kmod[i] * cos(Kphi[j]);
			Kzvector[j * Ngarmonik_mod + i] = Kmod[i] * sin(Kphi[j]);
		}
	}
}