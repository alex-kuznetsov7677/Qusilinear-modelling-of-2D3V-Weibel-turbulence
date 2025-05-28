// 1. Главный заголовок (должен быть первым)
#include "OutputWriter.h"

// 2. Стандартные библиотеки (только то, что реально используется)
#include <cmath>       // для std::sqrt, std::norm, std::abs
#include <complex>     // для std::complex, std::real

// 3. Локальные заголовки
#include "Integrals.h"  // для ENERGY_X, ENERGY_Y, DENSITY
#include "my_variables.h" 

void InitializeOutputFiles(OutputFiles &files, const std::map<std::string, double>& config) {
    files.CMFt100000.open("output/DESCRIPTIONRe[B]__170kiam.txt");
    
    files.CMFt100000 << "=== CONFIGURATION PARAMETERS ===\n";
  	for (const auto& pair : config) {
        files.CMFt100000 << pair.first << " = " << pair.second << "\n";
    }
	
	files.CMFt2000.open("output/abs[Bz]__170kiamMATLAB.txt");
	files.CMFt2001.open("output/abs[Bx]__170kiamMATLAB.txt");
	files.CMFt.open("output/output__170kiam.txt");
	files.CMFt0.open("output/ReFUNC000_170kiam.txt");
	files.CMFt1301.open("output/Re,Wy0]__170kiam.txt");
	files.CMFt1303.open("output/Re,Wx0]__170kiam.txt");
	files.CMFt1302.open("output/Im,Areal]__170kiam.txt");
	files.CMFt30.open("output/Re[BkSREDN]__170kiam.txt");
	files.CMFt31.open("output/Re[AREALWyWx]_170kiam.txt");
	files.CMFt32.open("output/Re[Time]__170kiam.txt");
}
void OutputFiles::log_execution_time(double time){
	CMFt << "Execution time:" << time << '\n';
}
void OutputFiles::log_number_MPI_processes(int number){
	CMFt << "Number of MPIprocesses:" << number << '\n';
}
void OutputFiles::log_fields(
    double time,
    const std::vector<std::complex<double>>& Bx,
    const std::vector<std::complex<double>>& Bz,
    int Ngarmonik,
    const std::vector<std::complex<double>>& f0k2cel,
    double NU_effective
) {
    CMFt << time << '\t';
    double Bmodkvadr = 0.;
    for (int it = 0; it < Ngarmonik; ++it) {
        Bmodkvadr = Bmodkvadr + pow(abs(Bx[it]), 2)+ pow(abs(Bz[it]), 2);
    }
    double Bmod = std::sqrt(Bmodkvadr);
    double VVx0 = real(ENERGY_X(f0k2cel));
    double VVy0 = real(ENERGY_Y(f0k2cel));
    double VVperp0 = real(ENERGY_PERP(f0k2cel));

	double ArealWydelWx = (VVy0) / (VVx0)-1;
	double Areal = 2. * (VVy0) / (VVperp0)-1;
    CMFt << real(DENSITY(f0k2cel)) << '\t' << Areal << '\t' << Bmod << '\t' << NU_effective << '\n';
	for (double p = 0; p < Ngarmonik; p++) {
		CMFt2000 << real(Bz[p]) << " ";
		CMFt2001 << real(Bx[p]) << " ";

	}
	CMFt2000 << ";";
	CMFt2001 << ";";
    CMFt30 << Bmod << ",";
	CMFt31 << ArealWydelWx << ",";
    CMFt1302 << Areal << ",";
    CMFt1301 << real(ENERGY_Z(f0k2cel)) << ",";
    CMFt1301 << VVy0 << ",";
    CMFt1303 << VVx0 << ",";
}

void OutputFiles::log_distribution(
    const std::vector<std::complex<double>>& f0k2cel,std::vector<double>& f0_proect
) {
	CMFt0 << '\n' << time << '\n';
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
void OutputFiles::close_all() {
    CMFt100000.close();
    CMFt2000.close();
    CMFt2001.close();
    CMFt1301.close();
    CMFt1303.close();
    CMFt1302.close();
    CMFt30.close();
    CMFt31.close();
    CMFt32.close();
    CMFt.close();
    CMFt0.close();
}
