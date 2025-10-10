#include "CRADLE/Utilities.hh"
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
    std::complex<double>* cConst = new std::complex<double>[8]; //{cs,csp,ct,ctp,cv,cvp,ca,cap}
    
    //example of a betaDecay parameters
    int Z = 14; //close to 13.7 = 0.1*alpha^-1
    int betaType = 1;
    double Q = 2000;  
   
    
    std::string LambdaTest(double j_max){
      std::stringstream file_content;
      for (double j_i = 0.; j_i <= j_max; j_i+=0.5){
	for (double j_f = 0.; j_f <= j_max; j_f+=0.5){
	  double lambda = polarisation::lambdaJiJf_factor(j_i,j_f);
	  file_content << std::setprecision(4) << lambda << "\t";
	}
	std::cout << j_i << std::endl;
	file_content << '\n';
      }
      return file_content.str();
    }

    std::string DoubleVarTest(int c_const1, int c_const2, bool complexValued, double mf, double mgt, double j_in, double j_f) {
      std::stringstream file_content;
      for (int i = 0; i < 8; i++){
	if (complexValued){
	  //using 
	  (cConst + i)->real(i == c_const1 ? 0.8 : (i == c_const2 ? 0.6 : 0));
	  (cConst + i)->imag(i == c_const1 ? 0.6 : (i == c_const2 ? -0.8 : 0));
	} else {
	  (cConst + i)->imag(0);
	  (cConst + i)->real(i == c_const1 ? 1 : (i == c_const2 ? 1 : 0));
	}
      }
      if ((mf == 0) || (mgt == 0)){
	double xi = utilities::CalculateXiBetaDecay(cConst[0], cConst[1], cConst[2], cConst[3], cConst[4], cConst[5], cConst[6], cConst[7], mf, mgt);
	std::cout << "Xi : " << xi << std::endl;
      }
      double* max_A = polarisation::MaximumA(cConst[0], cConst[1], cConst[2], cConst[3], cConst[4], cConst[5], cConst[6], cConst[7], mf, mgt, j_in, j_f, betaType, Z, Q);
      std::cout << "\tMaximum of A: " << *(max_A+1) << " at " << *(max_A) << " kev " << std::endl;
      

      double* max_B = polarisation::MaximumB(cConst[0], cConst[1], cConst[2], cConst[3], cConst[4], cConst[5], cConst[6], cConst[7], mf, mgt, j_in, j_f, betaType, Z, Q);
      std::cout << "\tMaximum of B: " << *(max_B+1) << " at " << *(max_B) << " kev " << std::endl;
      

      double* max_D = polarisation::MaximumD(cConst[0], cConst[1], cConst[2], cConst[3], cConst[4], cConst[5], cConst[6], cConst[7], mf, mgt, j_in, j_f, betaType, Z, Q);
      std::cout << "\tMaximum of D: " << *(max_D+1) << " at " << *(max_D) << " kev " << std::endl;
      
      bool changeA = false;
      bool changeB = false;
      bool changeD = false;
      for (int kin_en = 0; kin_en < Q; kin_en++){
	double en = polarisation::EMASSC2 + kin_en + 0.1; //last is to avoid 1/0 errors
	file_content << std::setprecision(5) << en << '\t';
	double A = polarisation::CalculateBetaAssymetry(cConst[0], cConst[1], cConst[2], cConst[3], cConst[4], cConst[5], cConst[6], cConst[7], mf, mgt, j_in, j_f, betaType, Z, en);
	if ((*(max_A+1) != std::abs(A)) && (not changeA)) changeA = true;
	file_content << std::setprecision(5) << A << '\t';
	double B = polarisation::CalculateNeutrinoAssymetry(cConst[0], cConst[1], cConst[2], cConst[3], cConst[4], cConst[5], cConst[6], cConst[7], mf, mgt, j_in, j_f, betaType, Z, en);
	if ((*(max_B+1) != std::abs(B)) && (not changeB)) changeB = true;
	file_content << std::setprecision(5) << B << '\t';
	double D = polarisation::CalculateDTripleCorrelation(cConst[0], cConst[1], cConst[2], cConst[3], cConst[4], cConst[5], cConst[6], cConst[7], mf, mgt, j_in, j_f, betaType, Z, en);
	if ((*(max_D+1) != std::abs(D)) && (not changeD)) changeD = true;
	file_content << std::setprecision(5) << D << '\n';
			}
      std::cout << "\tConstant A: " << (changeA ? "No" : "Yes") << std::endl;
      std::cout << "\tConstant B: " << (changeB ? "No" : "Yes") << std::endl;
      std::cout << "\tConstant D: " << (changeD ? "No" : "Yes") << std::endl;
      delete max_A;
      delete max_B;      
      delete max_D;
      return file_content.str();
    }

    void RunLambdaTest(){
      std::ofstream fileStream;
      fileStream.open("lambda_output.txt");
      fileStream << CRADLE::test::LambdaTest(10);
      fileStream.flush();
      fileStream.close();
    }

    std::string cteNames[8] = {"cs","csp","ct","ctp","cv","cvp","ca","cap"};
    
    void RunCConstGamovTellerTest(){
      std::ofstream fileStream;
      std::cout << "\t--------------------" << std::endl;
      std::cout << "\t  Real Coefficients " << std::endl;
      std::cout << "\t--------------------" << std::endl;
      for (int i = 0; i < 8; i++){
	for (int j = i+1; j < 8; j++){
	  std::stringstream fileNameSS;
	  std::cout << "Non-zero variables: " << cteNames[i] << ", " << cteNames[j] << std::endl;  
	  fileNameSS << cteNames[i] << cteNames[j] << "_real_gt_12.txt";
	  fileStream.open(fileNameSS.str());
	  fileStream << DoubleVarTest(i, j, false, 0., 1., 1., 2.);
	  fileStream.flush();
	  fileStream.close();
	  fileNameSS.str("");
	}
      }
      std::cout << "\t--------------------" << std::endl;
      std::cout << "\tComplex Coefficients" << std::endl;
      std::cout << "\t--------------------" << std::endl;
      for (int i = 0; i < 8; i++){
	for (int j = i+1; j < 8; j++){
	  std::stringstream fileNameSS;
	  std::cout << "Non-zero variables: " << cteNames[i] << ", " << cteNames[j] << std::endl;  
	  fileNameSS << cteNames[i] << cteNames[j] << "_comp_gt_12.txt";
	  fileStream.open(fileNameSS.str());
	  fileStream << DoubleVarTest(i, j, true, 0., 1., 1., 2.);
	  fileStream.flush();
	  fileStream.close();
	  fileStream.close();
	  fileNameSS.str("");
	}
      }
    }
    void RunCConstMixedTest(){
      std::ofstream fileStream;
      std::cout << "\t--------------------" << std::endl;
      std::cout << "\t  Real Coefficients " << std::endl;
      std::cout << "\t--------------------" << std::endl;
      for (int i = 0; i < 8; i++){
	for (int j = i+1; j < 8; j++){
	  std::stringstream fileNameSS;
	  std::cout << "Non-zero variables: " << cteNames[i] << ", " << cteNames[j] << std::endl;  
	  fileNameSS << cteNames[i] << cteNames[j] << "_real_mixed_11.txt";
	  fileStream.open(fileNameSS.str());
	  fileStream << DoubleVarTest(i, j, false, 1., 1., 1., 1.);
	  fileStream.flush();
	  fileStream.close();
	  fileNameSS.str("");
	}
      }
      std::cout << "\t--------------------" << std::endl;
      std::cout << "\tComplex Coefficients" << std::endl;
      std::cout << "\t--------------------" << std::endl;
      for (int i = 0; i < 8; i++){
	for (int j = i+1; j < 8; j++){
	  std::stringstream fileNameSS;
	  std::cout << "Non-zero variables: " << cteNames[i] << ", " << cteNames[j] << std::endl;  
	  fileNameSS << cteNames[i] << cteNames[j] << "_comp_mixed_11.txt";
	  fileStream.open(fileNameSS.str());
	  fileStream << DoubleVarTest(i, j, true, 1., 1., 1., 1.);
	  fileStream.flush();
	  fileStream.close();
	  fileNameSS.str("");
	}
      }
    }
    
  }//end of test namespace
}
