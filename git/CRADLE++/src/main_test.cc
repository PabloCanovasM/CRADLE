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

    std::string double_var_test(int c_const1, int c_const2, bool complexValued, double mf, double mgt, double j_in, double j_f) {
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

    inline std::string maximum_inspection_test(double a, double A, double B, double D, double E, int zRes, int phiRes, bool giveMaximum){
      std::stringstream file_content;
      for (int z_e = -zRes; z_e <= zRes; z_e++){
	for (int z_enu = -zRes; z_enu <= zRes; z_enu++){
	  for (int phi = 0; phi < phiRes; phi++){
	    double angCorrFactor = polarisation::CalculateAngularCorrelationFactor(a, 0, A, B, D, E, ((double) z_e)/zRes, ((double) z_enu)/zRes, phi*2*polarisation::PI/phiRes);
	    angCorrFactor = std::abs(angCorrFactor) > 1e-4 ? angCorrFactor : 0;  
	    file_content << std::fixed << std::setprecision(4) << angCorrFactor << '\t';
	  }
	  file_content << '\n'; 
	}
	if (z_e%5 == 0) std::cout << "cos(theta_e) = "<< ((double) z_e)/zRes << std::endl;	
      }
      if (giveMaximum){
	std::cout << "Finding maximum" << std::endl;
	double maxAngCorrFactor = polarisation::MaximumAngCorrFactor(a,0,A,B,D,E);
	std::cout << "Maximum value: " << maxAngCorrFactor << std::endl;
      }	
      return file_content.str();      
    }
  }
}

int main(){
  std::ofstream fileStream;
  
  bool lambda_test = false;
  bool double_cte_gt_test = false;
  bool double_cte_mixed_test = false;
  bool single_var_ang_corr_test = true;
  bool double_var_ang_corr_test = true;
  
  if (lambda_test){
      fileStream.open("lambda_output.txt");
      fileStream << CRADLE::test::LambdaTest(10);
      fileStream.flush();
      fileStream.close();
  }

  std::string cteNames[8] = {"cs","csp","ct","ctp","cv","cvp","ca","cap"};

  if (double_cte_gt_test){
    std::cout << "\t--------------------" << std::endl;
    std::cout << "\t  Real Coefficients " << std::endl;
    std::cout << "\t--------------------" << std::endl;
    for (int i = 0; i < 8; i++){
      for (int j = i+1; j < 8; j++){
	std::stringstream fileNameSS;
	std::cout << "Non-zero variables: " << cteNames[i] << ", " << cteNames[j] << std::endl;  
	fileNameSS << cteNames[i] << cteNames[j] << "_real_gt_12.txt";
	fileStream.open(fileNameSS.str());
	fileStream << CRADLE::test::double_var_test(i, j, false, 0., 1., 1., 2.);
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
	fileStream << CRADLE::test::double_var_test(i, j, true, 0., 1., 1., 2.);
	fileStream.flush();
	fileStream.close();
	fileStream.close();
	fileNameSS.str("");
      }
    }
  }

  if (double_cte_mixed_test){
    std::cout << "\t--------------------" << std::endl;
    std::cout << "\t  Real Coefficients " << std::endl;
    std::cout << "\t--------------------" << std::endl;
    for (int i = 0; i < 8; i++){
      for (int j = i+1; j < 8; j++){
	std::stringstream fileNameSS;
	std::cout << "Non-zero variables: " << cteNames[i] << ", " << cteNames[j] << std::endl;  
	fileNameSS << cteNames[i] << cteNames[j] << "_real_mixed_11.txt";
	fileStream.open(fileNameSS.str());
	fileStream << CRADLE::test::double_var_test(i, j, false, 1., 1., 1., 1.);
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
	fileStream << CRADLE::test::double_var_test(i, j, true, 1., 1., 1., 1.);
	fileStream.flush();
	fileStream.close();
	fileNameSS.str("");
      }
    }
  }

  if (single_var_ang_corr_test){
    double varList[4];
    std::string varNames[4] = {"a","A","B","D"};
    for (int i = 0; i < 4; i++){
	std::stringstream fileNameSS;
	for (int j = 0; j < 4; j++)
	  varList[j] = j == i ? 1 : 0;
	std::cout << "Non-zero variable: " << varNames[i] << std::endl;  
	fileNameSS << varNames[i] << "_simple_pos.txt";
	fileStream.open(fileNameSS.str());
	fileStream << CRADLE::test::maximum_inspection_test(varList[0],varList[1],varList[2],varList[3], 5000, 10, 24, false);
	fileStream.flush();
	fileStream.close();
	fileNameSS.str("");
      }
  }

  if (double_var_ang_corr_test){
    double varList[3];
    std::string varNames[3] = {"a","A","D"};
    std::stringstream fileNameSS;
    for (int i = 0; i < 3; i++){
      for (int j = 0; j < 3; j++)
	  varList[j] = j == i ? 1 : 0;
      int zRes = i == 2 ? 40 : 10;
      std::cout << "Non Zero " << varNames[i] << " and B" << std::endl;  
      std::cout << "Positive " << varNames[i] << ", positive B" << std::endl;
      std::cout << "Low Energy" << std::endl;
      fileNameSS << "pos" << varNames[i] << "_posB_lowE.txt";
      fileStream.open(fileNameSS.str());
      fileStream << CRADLE::test::maximum_inspection_test(varList[0], varList[1], 1, varList[2], 520, zRes, 24, true);
      fileStream.flush();
      fileStream.close();
      fileNameSS.str("");
      std::cout << "High Energy" << std::endl;
      fileNameSS << "pos" << varNames[i] << "_posB_hiE.txt"; 
      fileStream.open(fileNameSS.str());
      fileStream << CRADLE::test::maximum_inspection_test(varList[0], varList[1], 1, varList[2], 5000, zRes, 24, true);
      fileStream.flush();
      fileStream.close();
      fileNameSS.str("");  
      std::cout << "Positive " << varNames[i] << ", negative B" << std::endl;
      std::cout << "Low Energy" << std::endl; 
      fileNameSS << "pos" << varNames[i] << "_negB_lowE.txt";
      fileStream.open(fileNameSS.str());
      fileStream << CRADLE::test::maximum_inspection_test(varList[0], varList[1], -1, varList[2], 520, zRes, 24, true);
      fileStream.flush();
      fileStream.close();
      fileNameSS.str("");
      std::cout << "High Energy" << std::endl;
      fileNameSS << "pos" << varNames[i] << "_negB_hiE.txt"; 
      fileStream.open(fileNameSS.str());
      fileStream << CRADLE::test::maximum_inspection_test(varList[0], varList[1], -1 , varList[2], 5000, zRes, 24, true);
      fileStream.flush();
      fileStream.close();
      fileNameSS.str("");
    }
    
  }
  
  delete CRADLE::test::cConst;
}
