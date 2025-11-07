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
    void RunLambdaTest();
    void RunCConstGamovTellerTest();
    void RunCConstMixedTest();
    void RunCConst2FermiTest();
    void RunCConst2GamovTellerTest();
    void RunSingleVarAngCorrTest();
    void RunDoubleVarBAngCorrTest(int);
    void RunDoubleVarAAngCorrTest(int);
    void RunDoubleVarcAngCorrTest(int);
    void RunDoubleVarDAngCorrTest(int);
    void RunTripleVarAngCorrTest(int);
    void RunMaximumAngCorrTest();
    void RunDefaultVNRSamplingTest(double, double, double, double, double, double, double, int, std::string);
    void RunDoubleVarBVNRSamplingTest();
    void RunDoubleVarAVNRSamplingTest();
    void RunDoubleVarDVNRSamplingTest();
    void RunDoubleVarcVNRSamplingTest();
    void RunMinimumMaximumFromCouplingCTest();
    
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
  }
}

int main(int argc, char **argv) {
  std::vector<std::string> args(argv + 1, argv + argc);
  std::string testName;
  std::string genParam = "";

  //von Neumann specific parameters
  double a = 0;
  double b = 0;
  double c = 0;
  double A = 0;
  double B = 0;
  double D = 0;
  double E = 1000;
  int nDecays = 1000;
  std::string oFile = "output.txt";

  // Loop over command-line args
  testName = *args.begin(); //can be -h or --help, it will get detected later
  for (auto i = args.begin(); i != args.end(); ++i) {
    if (*i == "-h" || *i == "--help") {
      std::cout << "Syntax: Test <name-of-test> -p <general-param>" << std::endl;
      std::cout << "Lambda Test -> lambda" << std::endl;
      std::cout << "--------------------------------------------------" << std::endl;
      std::cout << "Coupling Constant Test -> c_const" << std::endl;
      std::cout << "-p gt for Gamov-Teller ABD, -p  m for Mixed ABD" << std::endl;
      std::cout << "-p gt2 for Gamov-Teller abc, -p f2 for Fermi abc" << std::endl;
      std::cout << "--------------------------------------------------" << std::endl;
      std::cout << "Angular Correlation Factor Test -> ang_corr_factor" << std::endl;
      std::cout << "Parameters:\n";
      std::cout << "-p 1 for 1 variable test\n";
      std::cout << "-p 2 for 2 variable test, non-zero B and other;\n";
      std::cout << "-p 3 for 2 variable test, non-zero A and one of a or D;\n";
      std::cout << "-p 4 for 2 variable test, non-zero D and a;\n";
      std::cout << "-p 5 for 2 variable test, non-zero c and one of a, A or D;\n";
      std::cout << "-p 6 for triple aAB test\n";
      std::cout << "-p <20|30|40|60> gives maximum computed with grid search\n";
      std::cout << "-p <21|31|41|51|61> gives maximum computed analytically" << std::endl;
      std::cout << "--------------------------------------------------" << std::endl;
      std::cout << "Maximum Angular Correlation -> max_f" << std::endl;
      std::cout << "--------------------------------------------------" << std::endl;
      std::cout << "Von Neumann Rejection Sampling -> sampling" << std::endl;
      std::cout << "Syntax: Test sampling -a <a> -b <b> -c <c> -A <A> -B <B> -D <D> -E <E> -N <n-decays> -o <output-file> -p <premade-test>" << std::endl;
      std::cout << "Premade Test:\n";
      std::cout << "-p 1 for 2 variable test, non-zero B and other;\n";
      std::cout << "-p 2 for 2 variable test, non-zero c and one of a, A or D;\n";
      std::cout << "-p 3 for 2 variable test, non-zero A and one of a or D;\n";
      std::cout << "-p 4 for 2 variable test, non-zero D and a" << std::endl;
      std::cout << "--------------------------------------------------" << std::endl;
      std::cout << "MinMax Angular Correlation -> minmax_f" << std::endl;
      std::cout << "Syntax Test minmax_f -o <output-file> -Q <Q> -M <M_F> <M_GT> -J <j_i> <j_f> -Z <+-Z>\n";
      std::cout << "-C <CS.Re> <CS.Im> <CSP.Re> <CSP.Im> <CT.Re> <CT.Im> <CTP.Re> <CTP.Im>\n";
      std::cout << "<CV.Re> <CV.Im> <CVP.Re> <CVP.Im> <CA.Re> <CA.Im> <CAP.Re> <CAP.Im>\n";
      return 0;
    }
    if (testName == "sampling"){
       if (*i == "-a")
	 a = std::stod(*(++i));
       else if (*i == "-b")
	 b = std::stod(*(++i));
       else if (*i == "-c")
	 c = std::stod(*(++i));
       else if (*i == "-A")
	 A = std::stod(*(++i));
       else if (*i == "-B")
	 B = std::stod(*(++i));
       else if (*i == "-D")
	 D = std::stod(*(++i));
       else if (*i == "-E")
	 E = std::stod(*(++i));
       else if (*i == "-N")
	 nDecays = std::stod(*(++i));
       else if (*i == "-o")
	 oFile = *(++i);
    } else if (testName == "minmax_f"){
      if (*i == "-Q")
	CRADLE::test::minMax::Q = std::stod(*(++i));
      else if (*i == "-M"){
	CRADLE::test::minMax::m_f = std::stod(*(++i));
	CRADLE::test::minMax::m_gt = std::stod(*(++i));
      }else if (*i == "-J"){
	CRADLE::test::minMax::j_i = std::stod(*(++i));
	CRADLE::test::minMax::j_f = std::stod(*(++i));
      }else if (*i == "-Z")
	CRADLE::test::minMax::Z = std::stoi(*(++i));
      else if (*i == "-C"){
	CRADLE::test::minMax::CS.real(std::stod(*(++i)));
	CRADLE::test::minMax::CS.imag(std::stod(*(++i)));
	CRADLE::test::minMax::CSP.real(std::stod(*(++i)));
	CRADLE::test::minMax::CSP.imag(std::stod(*(++i)));
	CRADLE::test::minMax::CT.real(std::stod(*(++i)));
	CRADLE::test::minMax::CT.imag(std::stod(*(++i)));
	CRADLE::test::minMax::CTP.real(std::stod(*(++i)));
	CRADLE::test::minMax::CTP.imag(std::stod(*(++i)));
	CRADLE::test::minMax::CV.real(std::stod(*(++i)));
	CRADLE::test::minMax::CV.imag(std::stod(*(++i)));
	CRADLE::test::minMax::CVP.real(std::stod(*(++i)));
	CRADLE::test::minMax::CVP.imag(std::stod(*(++i)));
	CRADLE::test::minMax::CA.real(std::stod(*(++i)));
	CRADLE::test::minMax::CA.imag(std::stod(*(++i)));
	CRADLE::test::minMax::CAP.real(std::stod(*(++i)));
	CRADLE::test::minMax::CAP.imag(std::stod(*(++i)));
      } else if (*i == "-o")
	CRADLE::test::minMax::oFile = *(++i);
    }
    if (*i == "-p")
      genParam = *(++i);
  }
  
  bool lambda_test = false;
  bool double_cte_gt_test = false;
  bool double_cte_mixed_test = false;
  bool double_cte2_gt_test = false;
  bool double_cte2_f_test = false;
  bool single_var_ang_corr_test = false;
  bool double_var_B_ang_corr_test = false;
  bool double_var_A_ang_corr_test = false;
  bool double_var_D_ang_corr_test = false;
  bool double_var_c_ang_corr_test = false;
  bool triple_var_aAB_ang_corr_test = false;
  bool maximum_ang_coor_factor_test = false;
  int show_maximum = 0;
  bool vonNeumann_reject_sampling_test = false;
  bool minmax_ang_coor_factor_test = false;

  if (testName == "lambda") lambda_test = true;
  else if (testName == "max_f") maximum_ang_coor_factor_test = true;
  else if (testName == "c_const"){
    if (genParam == "gt") double_cte_gt_test = true;
    else if (genParam == "m") double_cte_mixed_test = true;
    else if (genParam == "f2") double_cte2_f_test = true;
    else if (genParam == "gt2") double_cte2_gt_test = true;
    else{
      std::cout << "Parameter does not match keyword of any test" << std::endl;
      return 0;
    }
  }
  else if (testName == "ang_corr_factor"){
    int param_id = std::stoi(genParam);
    switch(param_id) {
    case 1:
      single_var_ang_corr_test = true;
      break;
    case 2:
      double_var_B_ang_corr_test = true;
      break;
    case 20:
      double_var_B_ang_corr_test = true;
      show_maximum = 1;
      break;
    case 21:
      double_var_B_ang_corr_test = true;
      show_maximum = 2;
      break;
    case 3:
      double_var_A_ang_corr_test = true;
      break;
    case 30:
      double_var_A_ang_corr_test = true;
      show_maximum = 1;
      break;
    case 31:
      double_var_A_ang_corr_test = true;
      show_maximum = 2;
      break;
    case 4:
      double_var_D_ang_corr_test = true;
      break;
    case 40:
      double_var_D_ang_corr_test = true;
      show_maximum = 1;
      break;
    case 41:
      double_var_D_ang_corr_test = true;
      show_maximum = 2;
      break;
    case 5:
      double_var_c_ang_corr_test = true;
      break;
    case 51:
      double_var_c_ang_corr_test = true;
      show_maximum = 2;
      break;
    case 6:
      triple_var_aAB_ang_corr_test = true;
      break;
    case 60:
      triple_var_aAB_ang_corr_test = true;
      show_maximum = 1;
      break;
    case 61:
      triple_var_aAB_ang_corr_test = true;
      show_maximum = 2;
      break;
    default:
      std::cout << "Parameter does not match id of any test" << std::endl;
      return 0;
    }
  }
  else if (testName == "sampling"){
    vonNeumann_reject_sampling_test = true;
  } else if (testName == "minmax_f"){
    minmax_ang_coor_factor_test = true;
  }
  else{
    std::cout << "Input does not match keyword of any test" << std::endl;
    return 0;
  }
  
  if (lambda_test){
    CRADLE::test::RunLambdaTest();
  }
  
  if (double_cte_gt_test){
    CRADLE::test::RunCConstGamovTellerTest();
  }

  if (double_cte_mixed_test){
    CRADLE::test::RunCConstMixedTest();
  }

  if (double_cte2_gt_test){
    CRADLE::test::RunCConst2GamovTellerTest();
  }

  if (double_cte2_f_test){
    CRADLE::test::RunCConst2FermiTest();
  }

  if (single_var_ang_corr_test){
    CRADLE::test::RunSingleVarAngCorrTest();
  }

  if (double_var_B_ang_corr_test){
    CRADLE::test::RunDoubleVarBAngCorrTest(show_maximum);
  }
  
  if (double_var_A_ang_corr_test){
    CRADLE::test::RunDoubleVarAAngCorrTest(show_maximum);
  }
  
  if (double_var_D_ang_corr_test){
    CRADLE::test::RunDoubleVarDAngCorrTest(show_maximum);
  }

  if (triple_var_aAB_ang_corr_test){
    CRADLE::test::RunTripleVarAngCorrTest(show_maximum);
  }

  if (double_var_c_ang_corr_test){
    CRADLE::test::RunDoubleVarcAngCorrTest(show_maximum);
  }

  if (maximum_ang_coor_factor_test){
    CRADLE::test::RunMaximumAngCorrTest();
  }

  if (minmax_ang_coor_factor_test){
    CRADLE::test::MinimumMaximumFromCouplingCTest();
  }

  if (vonNeumann_reject_sampling_test){
    if (genParam.empty()){
      std::cout << "a = " << a << "\n";
      std::cout << "b = " << b << "\n";
      std::cout << "c = " << c << "\n";
      std::cout << "A = " << A << "\n";
      std::cout << "B = " << B << "\n";
      std::cout << "D = " << D << "\n";
      std::cout << "E = " << E << "\n";
      std::cout << "Number of decays = " << nDecays << std::endl;

      CRADLE::test::RunDefaultVNRSamplingTest(a, b, c, A, B, D, E, nDecays, oFile);
      std::cout << "Results saved in " << oFile << std::endl; 
    } else {
      int param_id = std::stoi(genParam);
      switch(param_id) {
      case 1:
	CRADLE::test::RunDoubleVarBVNRSamplingTest();
	break;
      case 2:
	CRADLE::test::RunDoubleVarcVNRSamplingTest();
	break;
       case 3:
	CRADLE::test::RunDoubleVarAVNRSamplingTest();
	break;
      case 4:
	CRADLE::test::RunDoubleVarDVNRSamplingTest();
	break;
      default:
	std::cout << "Parameter does not match id of any test" << std::endl;
	return 0;
      } 
    }
  }
  
  return 0;
}  
