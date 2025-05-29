#include "Preparation.h"

#include <vector>
#include <cmath>

#include "my_variables.h"
#include "ConfigReader.h"
using std::vector;

void  Postanovka1(std::vector<double>& Kxvector, std::vector<double>& Kzvector) {
	auto config = ReadConfig("config.txt");
	vector <double> Kmod;
	vector <double> Kphi;
	int Ngarmonik_r = static_cast<int>(config["Ngarmonik_r"]);
	int Ngarmonik_phi =static_cast<int>(config["Ngarmonik_phi"]);
	//double Ngarmonik_phi = 10;
	for (double i = 0; i < Ngarmonik_r; i++) {
		Kmod.push_back(config["kmin_r"] + i * config["kstep_r"]);
	}
	for (double i = 0; i < Ngarmonik_phi; i++) {
		Kphi.push_back(i*3.1415926535/Ngarmonik_phi);
	}
	for (double i = 0; i < Ngarmonik_r; i++) {
		for (double j = 0; j < Ngarmonik_phi; j++) {
			Kxvector[j * Ngarmonik_r + i] = Kmod[i] * cos(Kphi[j]);
			Kzvector[j * Ngarmonik_r + i] = Kmod[i] * sin(Kphi[j]);
		}
	}
}