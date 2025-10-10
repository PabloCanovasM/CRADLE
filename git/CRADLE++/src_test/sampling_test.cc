#include "CRADLE/Polarisation.hh"

#include <sstream>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <complex>
#include <iomanip>

namespace CRADLE{
  namespace test{

    void RunDefaultVNRSampling(double a, double A, double B, double D, double E, int nDecays, std::string fileName){
      /*Function that performs Von Neumann Rejection Sampling with the given Decay
	parameters (a, A, B, D, E) to generate nDecays and outputs the kinematical
        variables cos(theta_e), cos(theta_nu), phi in the file fileName*/
      
      
      //compute the maximum of Angular Correlation Factor = F_max once
      std::ofstream fileStream;
      double F_max = 1;
      if (((a == 0) && (A == 0)) && ((B == 0) && (D == 0))){
	std::cout << "Trivial Angular Distribution" << std::endl;
      } else {
	std::cout << "Computing Maximum" << std::endl;
	F_max = polarisation::MaximumAngCorrFactor(a, A, B, D, E);
      }
	
      std::random_device rd;
      std::mt19937 generator(rd()); 
      std::uniform_real_distribution<double> distribution_cos(-1.0, 1.0); 
      std::uniform_real_distribution<double> distribution_phi(0,2*polarisation::PI);
      std::uniform_real_distribution<double> distribution_F(0.0, F_max);

      fileStream.open(fileName);
      int logRate = nDecays/20;
      int count = 0;
      while (count < nDecays){
	double angles[3] = {distribution_cos(generator), distribution_cos(generator), distribution_phi(generator)};
	double F = distribution_F(generator);
	double F_point = polarisation::CalculateAngularCorrelationFactor(a, A, B, D, E, angles[0], angles[1], angles[2]);
	if (F < F_point){
	  count++;
	  if (count%logRate == 0){
	    std::cout << "Decays so far: " << count << std::endl;
	  }
	  fileStream << std::fixed << std::setprecision(6) << angles[0] << '\t' << angles[1] << '\t' << angles[2] << '\n';
	}
      }
      fileStream.flush();
      fileStream.close();
    }
  }
}
