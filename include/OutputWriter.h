#pragma once

#include <fstream>   
#include <map>       
#include <string>  
#include <vector>    
#include <complex>   
struct OutputFiles {
    std::ofstream CMFt100000;
    std::ofstream CMFt2000,CMFt2001;
    std::ofstream CMFt1301, CMFt1303, CMFt1302;
    std::ofstream CMFt30, CMFt31;
    std::ofstream CMFt32, CMFt, CMFt0;

    void log_number_MPI_processes(int number);
    void log_execution_time(double time);
    void close_all();
    void log_fields(
    		double time,
    		const std::vector<std::complex<double>>& Bx,
    		const std::vector<std::complex<double>>& Bz,
    		int Ngarmonik,
    		const std::vector<std::complex<double>>& f0k2cel,
    		double NU_effective);

    void log_distribution(const std::vector<std::complex<double>>& f0k2cel,std::vector<double>& f0_proect);
};


void InitializeOutputFiles(OutputFiles &files, const std::map<std::string, double>& config);