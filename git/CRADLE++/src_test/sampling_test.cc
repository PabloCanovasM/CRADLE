#include "CRADLE/Polarisation.hh"

#include <sstream>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <complex>
#include <iomanip>
#include <boost/numeric/ublas/vector.hpp>

namespace CRADLE{
  namespace test{
    inline ublas::vector<double> RunDefaultVNRSampling(double a, double A, double B, double D, double E, double F_max){
      /*Function that performs Von Neumann Rejection Sampling with the given Decay
	parameters (a, A, B, D, E) to generate one decay and outputs the kinematical
        variables cos(theta_e), cos(theta_nu), phi in a list*/

      std::random_device rd;
      std::mt19937 generator(rd()); 
      std::uniform_real_distribution<double> distribution_cos(-1.0, 1.0); 
      std::uniform_real_distribution<double> distribution_phi(0,2*polarisation::PI);
      std::uniform_real_distribution<double> distribution_F(0.0, F_max);

      ublas::vector<double> decayVars(3);
      
      double F = F_max;
      double F_point = 0;
      while (F > F_point){
	double angles[3] = {distribution_cos(generator), distribution_cos(generator), distribution_phi(generator)};
	double F = distribution_F(generator);
	double F_point = polarisation::CalculateAngularCorrelationFactor(a, A, B, D, E, angles[0], angles[1], angles[2]);
	if (F < F_point){
	  for (ublas::vector<double>::iterator i = decayVars.begin(); i != decayVars.end(); i++){
	    int pos = i-decayVars.begin();
	    *i = angles[pos];
	  }
	  return decayVars;
	}
      } 
    }

    void RunDefaultVNRSamplingTest(double a, double A, double B, double D, double E, int nDecays, std::string fileName){
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
      fileStream.open(fileName);
      int logRate = nDecays/20;
      for(int count = 0; count < nDecays; count++){
	ublas::vector<double> angles = RunDefaultVNRSampling(a, A, B, D, E, F_max);
	if (count%logRate == 0){
	    std::cout << "Decays so far: " << count << std::endl;
	  }
	fileStream << std::fixed << std::setprecision(6) << angles[0] << '\t' << angles[1] << '\t' << angles[2] << '\n';
	}
      fileStream.flush();
      fileStream.close();
    }

    void RunDoubleVarBVNRSamplingTest(){
      double varList[3];
      std::string varNames[3] = {"a","A","D"};
      std::stringstream fileNameSS;
      for (int i = 0; i < 3; i++){
	for (int j = 0; j < 3; j++)
	  varList[j] = j == i ? 1 : 0;
	std::cout << "Non Zero " << varNames[i] << " and B" << std::endl;  
	std::cout << "Positive " << varNames[i] << ", positive B" << std::endl;
	std::cout << "Low Energy" << std::endl;
	fileNameSS << "sample_pos" << varNames[i] << "_posB_lowE.txt";
	RunDefaultVNRSamplingTest(varList[0], varList[1], 1, varList[2], 520, 300000, fileNameSS.str());
	fileNameSS.str("");
	std::cout << "Medium Energy" << std::endl;
	fileNameSS << "sample_pos" << varNames[i] << "_posB_medE.txt";
	RunDefaultVNRSamplingTest(varList[0], varList[1], 1, varList[2], 800, 300000, fileNameSS.str()); 
	fileNameSS.str("");  
	std::cout << "High Energy" << std::endl;
	fileNameSS << "pos" << varNames[i] << "_posB_hiE.txt"; 
	RunDefaultVNRSamplingTest(varList[0], varList[1], 1, varList[2], 5000, 300000, fileNameSS.str());
	fileNameSS.str("");
	std::cout << "High " << varNames[i] << std::endl;
	fileNameSS << "sample_pos" << varNames[i] << "_posB_hi" << varNames[i]  << ".txt";
	varList[i] = 2;
	RunDefaultVNRSamplingTest(varList[0], varList[1], 1, varList[2], 5000, 300000, fileNameSS.str());
	fileNameSS.str("");
	varList[i] = 5;
	std::cout << "Very High " << varNames[i] << std::endl;
	fileNameSS << "sample_pos" << varNames[i] << "_posB_vhi" << varNames[i]  << ".txt";
	RunDefaultVNRSamplingTest(varList[0], varList[1], 1, varList[2], 5000, 300000, fileNameSS.str());
	fileNameSS.str("");     
	varList[i] = 1; //reset for the negative part
	std::cout << "Positive " << varNames[i] << ", negative B" << std::endl;
	std::cout << "Low Energy" << std::endl; 
	fileNameSS << "sample_pos" << varNames[i] << "_negB_lowE.txt";
	RunDefaultVNRSamplingTest(varList[0], varList[1], -1, varList[2], 5000, 300000, fileNameSS.str());
	fileNameSS.str("");
	std::cout << "High Energy" << std::endl;
	fileNameSS << "sample_pos" << varNames[i] << "_negB_hiE.txt"; 
	RunDefaultVNRSamplingTest(varList[0], varList[1], -1, varList[2], 5000, 300000, fileNameSS.str());
	fileNameSS.str("");     
      }
    }

    void RunDoubleVarAVNRSamplingTest(){
      double varList[2];
      std::string varNames[2] = {"a","D"};
      std::stringstream fileNameSS;
      for (int i = 0; i < 2; i++){
	for (int j = 0; j < 2; j++)
	  varList[j] = j == i ? 1 : 0;
	double E = 100000; //very high, to make beta almost 1
	int N = 300000;
	std::cout << "Non Zero " << varNames[i] << " and A" << std::endl;
	double A_vals[5] = {0.1,0.5,1,2,5};
	std::string A_names[5] = {"xs","s","eq","l","xl"};
	for (int j = 0; j < 5; j++){
	  double A = A_vals[j];
	  std::cout << "A = +" << A << std::endl;  
	  fileNameSS << "pos" << varNames[i] << "_" << A_names[j]  << "posA.txt";
	  RunDefaultVNRSamplingTest(varList[0], A, 0, varList[1], E, N, fileNameSS.str());
	  fileNameSS.str("");
	  std::cout << "A = -" << A << std::endl;  
	  fileNameSS << "pos" << varNames[i] << "_" << A_names[j]  << "negA.txt";
	  RunDefaultVNRSamplingTest(varList[0], -A, 0, varList[1], E, N, fileNameSS.str());
	  fileNameSS.str("");
	}
      }
    }

    void RunDoubleVarDVNRSamplingTest(){
      std::stringstream fileNameSS;
      std::cout << "Non Zero a and D" << std::endl;
      double D_vals[5] = {0.25,0.5,1,2,4};
      std::string D_names[5] = {"xs","s","eq","l","xl"};
      double E = 100000;
      int N = 300000;
      for (int j = 0; j < 5; j++){
	double D = D_vals[j];
	std::cout << "D = +" << D << std::endl;  
	fileNameSS << "posa_" << D_names[j]  << "posD.txt";
	RunDefaultVNRSamplingTest(1, 0, 0, D, E, N, fileNameSS.str());
	fileNameSS.str("");
	std::cout << "D = -" << D << std::endl;  
	fileNameSS << "posa_" << D_names[j]  << "negD.txt";
	RunDefaultVNRSamplingTest(1, 0, 0, -D, E, N, fileNameSS.str());
	fileNameSS.str("");
      }
    }
  }
}
