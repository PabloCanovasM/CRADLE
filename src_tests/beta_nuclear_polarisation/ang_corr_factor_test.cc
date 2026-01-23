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
#include <boost/numeric/ublas/vector.hpp>

namespace CRADLE{
  namespace test{
    using namespace boost::numeric::ublas;
    namespace minMax{
      double Q;
      double j_i;
      double j_f;
      int Z;
      double m_f;
      double m_gt;
      std::complex<double> CS;
      std::complex<double> CT;
      std::complex<double> CV;
      std::complex<double> CA;
      std::complex<double> CSP;
      std::complex<double> CTP;
      std::complex<double> CVP;
      std::complex<double> CAP;
      std::string oFile = "output.txt";
    }

    inline double MinimumF(double a, double A, double B, double K, double znu){
      return -std::sqrt(std::pow(a*znu+A,2)+K*K*(1-znu*znu))+B*znu;
    }

    inline double AnalyticalMinimumAngCorrFactor(double a, double b, double c, double A, double B, double D, double E){
      /*Search of the minimum value analitically*/
      double beta = std::sqrt(1-utilities::EMASSC2*utilities::EMASSC2/E/E);
      //scale a, c, A and D by beta. Note that c is multiplied already by the alignment, and A, B and D are multiplied by J as well
      a *= beta;
      c *= beta;
      A *= beta;
      D *= beta;

      double K = std::sqrt(D*D+(a+c/3)*(a+c/3));
      double a_st = a-2.*c/3;

      double A_m = (a_st*a_st-K*K)*(a_st*a_st-K*K-B*B);
      double B_m = 2*a_st*A*(a_st*a_st-K*K-B*B);
      double C_m = a_st*a_st*A*A-A*A*B*B-B*B*K*K;
      double znu_m, znu_m2; //candidates for maximum
      vector<double> F_cand(4);

      //computing the values at the extrema of the interval
      F_cand(0) = MinimumF(a_st, A, B, K, 1);
      F_cand(1) = MinimumF(a_st, A, B, K, -1);
      F_cand(2) = 100;
      F_cand(3) = 100;

      if (A_m == 0) {
		znu_m = -C_m/B_m;
		if ((znu_m > -1) && (znu_m < 1)) {
			F_cand(2) = MinimumF(a_st, A, B, K, znu_m);
		}
		} else if (B_m == 0 && C_m == 0){
			F_cand(2) = MinimumF(a_st, A, B, K, 0); //cover some edge cases like only D non-zero that should never happen in reality
		} else {
		double det = B_m*B_m - 4*A_m*C_m;
			if ((det < 0) && (det > -1e-16)) det = 0; //avoid spurious cases
			if (det >= 0){
				znu_m = (-B_m+std::sqrt(det))/2/A_m;
				if ((znu_m > -1) && (znu_m < 1)) {
					F_cand(2) = MinimumF(a_st, A, B, K, znu_m);
				}
				znu_m2 = (-B_m-std::sqrt(det))/2/A_m;
				if ((znu_m2 > -1) && (znu_m2 < 1)) {
					F_cand(3) = MinimumF(a_st, A, B, K, znu_m2);
				}
			}
		}
		//std::cout << F_cand << std::endl;
		double F_min = *std::min_element(F_cand.begin(),F_cand.end());

		F_min += 1 + b*utilities::EMASSC2/E; //adding the constant terms
		return F_min;
    }
    
    inline std::string MaximumInspectionTest(double a, double c, double A, double B, double D, double E, int zeRes, int znuRes, int phiRes, int giveMaximum){
      	std::stringstream file_content;
      	for (int z_e = -zeRes; z_e <= zeRes; z_e++){
			for (int z_enu = -znuRes; z_enu <= znuRes; z_enu++){
				for (int phi = 0; phi < phiRes; phi++){
				double angCorrFactor = polarisation::CalculateAngularCorrelationFactor(a, 0, c, A, B, D, E, ((double) z_e)/zeRes, ((double) z_enu)/znuRes, phi*2*polarisation::PI/phiRes);
				angCorrFactor = std::abs(angCorrFactor) > 1e-6 ? angCorrFactor : 0;  
				file_content << std::fixed << std::setprecision(5) << angCorrFactor << '\t';
				}
				file_content << '\n'; 
			}	
			if (z_e%5 == 0) std::cout << "cos(theta_e) = "<< ((double) z_e)/zeRes << std::endl;	
		}
		if (giveMaximum == 1){
			std::cout << "Finding maximum" << std::endl;
			double maxAngCorrFactor = polarisation::MaximumAngCorrFactor(a, A, B, D, E);
			std::cout << "Maximum value: " << maxAngCorrFactor << std::endl;
		}	else if (giveMaximum == 2){
			std::cout << "Finding maximum" << std::endl;
			double maxAngCorrFactor = polarisation::AnalyticalMaximumAngCorrFactor(a, 0, c, A, B, D, E);
			ublas::vector<double> pos = polarisation::MaximumAngCorrFactorPos(a, 0, c, A, B, D, E);
			std::cout << "Maximum value: " << maxAngCorrFactor << std::endl;
			std::cout << "at z_e = " << pos(0) << ", z_nu = " << pos(1) << ", phi = " << pos(2) << std::endl;
		}	
      return file_content.str();      
    }

    inline std::string MaximumInspectionTest(double a, double c, double A, double B, double D, double E, int zRes, int phiRes, int giveMaximum){	
      return MaximumInspectionTest(a,c,A,B,D,E,zRes,zRes,phiRes,giveMaximum);      
    }

    

    void RunSingleVarAngCorrTest(){
      std::ofstream fileStream;
      double varList[5];
      std::string varNames[5] = {"a","c","A","B","D"};
      for (int i = 0; i < 5; i++){
		std::stringstream fileNameSS;
		for (int j = 0; j < 5; j++)
		varList[j] = j == i ? 1 : 0;
		std::cout << "Non-zero variable: " << varNames[i] << std::endl;  
		fileNameSS << varNames[i] << "_simple_pos.txt";
		fileStream.open(fileNameSS.str());
		fileStream << MaximumInspectionTest(varList[0],varList[1],varList[2],varList[3],varList[4], 5000, 10, 24, 0);
		fileStream.flush();
		fileStream.close();
		fileNameSS.str("");
      }
      std::cout << "Negative c" << std::endl;  
      fileStream.open("c_simple_neg.txt");
      fileStream << MaximumInspectionTest(0, -1, 0, 0, 0, 5000, 10, 24, 0);
      fileStream.flush();
      fileStream.close();
    }

    void RunDoubleVarBAngCorrTest(int compMax){
		std::ofstream fileStream;
		double varList[4];
		std::string varNames[4] = {"a","c","A","D"};
		std::stringstream fileNameSS;
		for (int i = 0; i < 4; i++){
			for (int j = 0; j < 4; j++)
			varList[j] = j == i ? 1 : 0;
			int znuRes = i == 3 ? 100 : 20;
			int zRes = i == 3 ? 40 : 20;
			int phiRes = 96;
			std::cout << "Non Zero " << varNames[i] << " and B" << std::endl;  
			std::cout << "Positive " << varNames[i] << ", positive B" << std::endl;
			std::cout << "Low Energy" << std::endl;
			fileNameSS << "pos" << varNames[i] << "_posB_lowE.txt";
			fileStream.open(fileNameSS.str());
			fileStream << MaximumInspectionTest(varList[0], varList[1], varList[2], 1, varList[3], 520, zRes, znuRes, phiRes, compMax);
			fileStream.flush();
			fileStream.close();
			fileNameSS.str("");
			std::cout << "Medium Energy" << std::endl;
			fileNameSS << "pos" << varNames[i] << "_posB_medE.txt"; 
			fileStream.open(fileNameSS.str());
			fileStream << MaximumInspectionTest(varList[0], varList[1], varList[2], 1, varList[3], 800, zRes, znuRes, phiRes, compMax);
			fileStream.flush();
			fileStream.close();
			fileNameSS.str("");  
			std::cout << "High Energy" << std::endl;
			fileNameSS << "pos" << varNames[i] << "_posB_hiE.txt"; 
			fileStream.open(fileNameSS.str());
			fileStream << MaximumInspectionTest(varList[0], varList[1], varList[2], 1, varList[3], 5000, zRes, znuRes, phiRes, compMax);
			fileStream.flush();
			fileStream.close();
			fileNameSS.str("");
			std::cout << "High " << varNames[i] << std::endl;
			fileNameSS << "pos" << varNames[i] << "_posB_hi" << varNames[i]  << ".txt";
			varList[i] = 2;
			fileStream.open(fileNameSS.str());
			fileStream << MaximumInspectionTest(varList[0], varList[1], varList[2], 1, varList[3], 5000, zRes, znuRes, phiRes, compMax);
			fileStream.flush();
			fileStream.close();
			fileNameSS.str("");
			varList[i] = 5;
			std::cout << "Very High " << varNames[i] << std::endl;
			fileNameSS << "pos" << varNames[i] << "_posB_vhi" << varNames[i]  << ".txt"; 
			fileStream.open(fileNameSS.str());
			fileStream << MaximumInspectionTest(varList[0], varList[1], varList[2], 1, varList[3], 5000, zRes, znuRes, phiRes, compMax);
			fileStream.flush();
			fileStream.close();
			fileNameSS.str("");     
			varList[i] = 1; //reset for the negative part
			std::cout << "Positive " << varNames[i] << ", negative B" << std::endl;
			std::cout << "Low Energy" << std::endl; 
			fileNameSS << "pos" << varNames[i] << "_negB_lowE.txt";
			fileStream.open(fileNameSS.str());
			fileStream << MaximumInspectionTest(varList[0], varList[1], varList[2], -1, varList[3], 520,  zRes, 24,  compMax);
			fileStream.flush();
			fileStream.close();
			fileNameSS.str("");
			std::cout << "High Energy" << std::endl;
			fileNameSS << "pos" << varNames[i] << "_negB_hiE.txt"; 
			fileStream.open(fileNameSS.str());
			fileStream << MaximumInspectionTest(varList[0], varList[1], varList[2], -1 , varList[3], 5000, zRes, 24, compMax);
			fileStream.flush();
			fileStream.close();
			fileNameSS.str("");     
      	}	
    }
    
	void RunDoubleVarcAngCorrTest(int compMax){
		std::ofstream fileStream;
		double varList[3];
		std::string varNames[3] = {"a","A","D"};
		std::stringstream fileNameSS;
		for (int i = 0; i < 3; i++){
			for (int j = 0; j < 3; j++)
			varList[j] = j == i ? 1 : 0;
			int zeRes = i == 2 ? 100 : 20;
			int zRes = i == 2 ? 40 : 20;
			int E = 100000; //very high, to make beta almost 1 
			std::cout << "Non Zero " << varNames[i] << " and c" << std::endl;
			double c_vals[5] = {0.25,0.5,1,2,4};
			std::string c_names[5] = {"xs","s","eq","l","xl"};
			for (int j = 0; j < 5; j++){
				double c = c_vals[j];
				std::cout << "c = +" << c << std::endl;  
				fileNameSS << "pos" << varNames[i] << "_" << c_names[j]  << "posc.txt";
				fileStream.open(fileNameSS.str());
				fileStream << MaximumInspectionTest(varList[0], c, varList[1], 0 , varList[2], E, zeRes, zRes, 96, compMax);
				fileStream.flush();
				fileStream.close();
				fileNameSS.str("");
				std::cout << "c = -" << c << std::endl;  
				fileNameSS << "pos" << varNames[i] << "_" << c_names[j]  << "negc.txt";
				fileStream.open(fileNameSS.str());
				fileStream << MaximumInspectionTest(varList[0], -c, varList[1], 0, varList[2], E, zeRes, zRes, 96, compMax);
				fileStream.flush();
				fileStream.close();
				fileNameSS.str("");
			}
      	}
    }
    
    void RunDoubleVarAAngCorrTest(int compMax){
		std::ofstream fileStream;
		double varList[2];
		std::string varNames[2] = {"a","D"};
		std::stringstream fileNameSS;
		for (int i = 0; i < 2; i++){
			for (int j = 0; j < 2; j++)
			varList[j] = j == i ? 1 : 0;
			int zeRes = i == 1 ? 100 : 20;
			int zRes = i == 1 ? 40 : 20;
			int E = 100000; //very high, to make beta almost 1 
			std::cout << "Non Zero " << varNames[i] << " and A" << std::endl;
			double A_vals[5] = {0.1,0.5,1,2,5};
			std::string A_names[5] = {"xs","s","eq","l","xl"};
			for (int j = 0; j < 5; j++){
				double A = A_vals[j];
				std::cout << "A = +" << A << std::endl;  
				fileNameSS << "pos" << varNames[i] << "_" << A_names[j]  << "posA.txt";
				fileStream.open(fileNameSS.str());
				fileStream << MaximumInspectionTest(varList[0], 0, A, 0, varList[1], E, zeRes, zRes, 96, compMax);
				fileStream.flush();
				fileStream.close();
				fileNameSS.str("");
				std::cout << "A = -" << A << std::endl;  
				fileNameSS << "pos" << varNames[i] << "_" << A_names[j]  << "negA.txt";
				fileStream.open(fileNameSS.str());
				fileStream << MaximumInspectionTest(varList[0], 0, -A, 0, varList[1], E, zeRes, zRes, 96, compMax);
				fileStream.flush();
				fileStream.close();
				fileNameSS.str("");
			}
    	}
    }
    
    void RunDoubleVarDAngCorrTest(int compMax){
		std::ofstream fileStream;
		std::stringstream fileNameSS;
		std::cout << "Non Zero a and D" << std::endl;
		double D_vals[5] = {0.25,0.5,1,2,4};
		std::string D_names[5] = {"xs","s","eq","l","xl"};
		double E = 100000;
		int zRes = 100;
		for (int j = 0; j < 5; j++){
			double D = D_vals[j];
			std::cout << "D = +" << D << std::endl;  
			fileNameSS << "posa_" << D_names[j]  << "posD.txt";
			fileStream.open(fileNameSS.str());
			fileStream << MaximumInspectionTest(1, 0, 0, 0, D, E, zRes, zRes, 96, compMax);
			fileStream.flush();
			fileStream.close();
			fileNameSS.str("");
			std::cout << "D = -" << D << std::endl;  
			fileNameSS << "posa_" << D_names[j]  << "negD.txt";
			fileStream.open(fileNameSS.str());
			fileStream << MaximumInspectionTest(1, 0, 0, 0, -D, E, zRes, zRes, 96, compMax);
			fileStream.flush();
			fileStream.close();
			fileNameSS.str("");
      	}
    }

    void RunTripleVarAngCorrTest(int compMax){
		std::ofstream fileStream;
		std::stringstream fileNameSS;
		std::cout << "Non Zero a, B, A" << std::endl;
		double AB_vals[5] = {0.25,0.5,0.6,0.8,1};
		double E_vals[5] = {520, 600, 800, 1500, 3000};
		double a = -1.0/3;
		for (double B : AB_vals){
			double A = -B;
			std::cout << "|A| = |B| = " << B << std::endl;
			for (double E : E_vals){
			std::cout << "E = " << E << std::endl;
			fileNameSS << std::setprecision(2) << "AB_" << B << "_E_" << std::setprecision(4) << E <<".txt";
			fileStream.open(fileNameSS.str());
			fileStream << MaximumInspectionTest(a, 0, A, B, 0, E, 100, 12, compMax);
			fileStream.flush();
			fileStream.close();
			fileNameSS.str("");
			}
		}
    }
    
    void RunMaximumAngCorrTest(){
		std::string varNames[5] = {"a","c","A","B","D"};
		double varList[5];
		std::ofstream fileStream;
		std::stringstream fileNameSS;
		double E = 100000;
		for (int i = 0; i < 5; i++){
			for (int j = 0; j < 5; j++) varList[j] = i == j ? 1 : 0; //setting rest of vars to 0 except the first variable of the test
			std::cout << "Non-Zero element: " << varNames[i] << std::endl;
			fileNameSS << varNames[i] << "_maxF_test.txt" ;
			fileStream.open(fileNameSS.str());
			for (int val = -10; val < 11; val += 1){
				varList[i] = val/10.; 
				double maxF = polarisation::AnalyticalMaximumAngCorrFactor(varList[0], 0, varList[1], varList[2], varList[3], varList[4], E);
				ublas::vector<double> maxpos = polarisation::MaximumAngCorrFactorPos(varList[0], 0, varList[1], varList[2], varList[3], varList[4], E);
				fileStream << std::fixed << std::setprecision(4) << varList[i] << "\t" <<  maxF;
				for (double ang : maxpos) fileStream << "\t" << ang;
				fileStream << "\n";
			}
			fileStream.flush();
			fileStream.close();
			fileNameSS.str("");
			varList[i] = 1; //reset first var for 2 var tests. 
			for (int j = i+1; j < 5; j++){
			std::cout << "Non-Zero elements: " << varNames[i] << " and " << varNames[j] << std::endl;
			fileNameSS << varNames[i] << varNames[j] << "_maxF_test.txt" ;
			fileStream.open(fileNameSS.str());
			for (int ratio = -30; ratio < 31; ratio++){
				if (ratio == 0) continue;
				varList[j] = ratio/10.;
				double maxF = polarisation::AnalyticalMaximumAngCorrFactor(varList[0], 0, varList[1], varList[2], varList[3], varList[4], E);
				ublas::vector<double> maxpos = polarisation::MaximumAngCorrFactorPos(varList[0], 0, varList[1], varList[2], varList[3], varList[4], E);
				fileStream << std::fixed << std::setprecision(4) << varList[j] << "\t" <<  maxF;
				for (double ang : maxpos){
					fileStream << "\t" << ang;
				}
				fileStream << "\n";
			}
			fileStream.flush();
			fileStream.close();
			fileNameSS.str("");
			varList[j] = 0; //reset
			}
		}
    }

    void RunMinimumMaximumFromCouplingCTest(){
		int Z = std::abs(minMax::Z);
		int betaType = (Z >= 0) -  (Z < 0); //n always decays beta-
		Z += betaType; //Z to be used in computations is always the final nuclei

		std::complex<double> CS = minMax::CS;
		std::complex<double> CSP = minMax::CSP;
		std::complex<double> CT = minMax::CT;
		std::complex<double> CTP = minMax::CTP;
		std::complex<double> CV = minMax::CV;
		std::complex<double> CVP = minMax::CVP;
		std::complex<double> CA = minMax::CA;
		std::complex<double> CAP = minMax::CAP;

		double m_f = minMax::m_f;
		double m_gt = minMax::m_gt;
		double j_i = minMax::j_i;
		double j_f = minMax::j_f;

		double xi = utilities::CalculateXiBetaDecay(CS, CSP, CT, CTP, CV, CVP, CA, CAP, m_f, m_gt);

		if (xi == 0){
			std::cout << "Impossible decay, terminating..." << std::endl;
			return;
		}

		std::ofstream fileStream;
		fileStream.open(minMax::oFile);

		double b = utilities::CalculateFierz(CS, CSP, CT, CTP, CV, CVP, CA, CAP, m_f, m_gt, std::nan("1"), std::nan("1"), Z, betaType); //energy independent
		for (double kinE = 1.0; kinE < minMax::Q; kinE += 0.5){
			double E = kinE + utilities::EMASSC2;
			double a = utilities::CalculateBetaNeutrinoAsymmetry(CS, CSP, CT, CTP, CV, CVP, CA, CAP, m_f, m_gt, std::nan("1"), std::nan("1"), E, Z, betaType);
			double A = polarisation::CalculateBetaAssymetry(CS, CSP, CT, CTP, CV, CVP, CA, CAP, m_f, m_gt, j_i, j_f, betaType, Z, E);
			double B = polarisation::CalculateNeutrinoAssymetry(CS, CSP, CT, CTP, CV, CVP, CA, CAP, m_f, m_gt, j_i, j_f, betaType, Z, E);
			double D = polarisation::CalculateDTripleCorrelation(CS, CSP, CT, CTP, CV, CVP, CA, CAP, m_f, m_gt, j_i, j_f, betaType, Z, E);
			double c = polarisation::CalculateAlignmentCorrelation(CS, CSP, CT, CTP, CV, CVP, CA, CAP, m_f, m_gt, j_i, j_f, betaType, Z, E);

			c = -c; //maximum alignment

			double max_F = polarisation::AnalyticalMaximumAngCorrFactor(a, b, c, A, B, D, E);
			double min_F = AnalyticalMinimumAngCorrFactor(a, b, c, A, B, D, E);

			fileStream << std::fixed << std::setprecision(4) << E << '\t' << max_F << '\t' << min_F << '\n';
		}

		fileStream.flush();
		fileStream.close();
		
	}
		
}//end of test namespace
}
