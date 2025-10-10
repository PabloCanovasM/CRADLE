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
    void RunSingleVarAngCorrTest();
    void RunDoubleVarBAngCorrTest();
    void RunDoubleVarAAngCorrTest();
    void RunDoubleVarDAngCorrTest();
    void RunTripleVarAngCorrTest();
    void RunDefaultVNRSampling(double, double, double, double, double, int, std::string);
  }
}

int main(int argc, char **argv) {
  std::vector<std::string> args(argv + 1, argv + argc);
  std::string testName;
  std::string genParam = "";
  // Von Neumann specific params

  std::string ofile;
  double a = 0;
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
      std::cout << "\n--------------------------------------------------\n" << std::endl;
      std::cout << "Coupling Constant Test -> c_const" << std::endl;
      std::cout << "-p gt for Gamov-Teller, -p  m for Mixed" << std::endl;
      std::cout << "\n--------------------------------------------------\n" << std::endl;
      std::cout << "Angular Correlation Factor Test -> ang_corr_factor" << std::endl;
      std::cout << "Parameters:\n";
      std::cout << "-p 1 for 1 variable test\n";
      std::cout << "-p <2|3|4> for 2 variable test with B, A or D;\n";
      std::cout << "-p 5 for triple aAB test" << std::endl;
      std::cout << "\n--------------------------------------------------\n";
      std::cout << "Syntax: Test sampling -a <a> -A <A> -B <B> -D <D> -E <E> -N <n-decays> -o <output-file> -p <premade-test>" << std::endl;
      std::cout << "Von Neumann Rejection Sampling -> sampling" << std::endl;
      return 0;
    }
    if (testName == "sampling"){

       if (*i == "-a")
	 a = std::stod(*(++i));
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
    }
    if (*i == "-p")
      genParam = *(++i);
  }
  
  bool lambda_test = false;
  bool double_cte_gt_test = false;
  bool double_cte_mixed_test = false;
  bool single_var_ang_corr_test = false;
  bool double_var_B_ang_corr_test = false;
  bool double_var_A_ang_corr_test = false;
  bool double_var_D_ang_corr_test = false;
  bool triple_var_aAB_ang_corr_test = false;
  bool vonNeumann_reject_sampling_test = false;

  if (testName == "lambda") lambda_test = true;
  else if (testName == "c_const"){
    if (genParam == "gt") double_cte_gt_test = true;
    else if (genParam == "m") double_cte_mixed_test = true;
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
    case 3:
      double_var_A_ang_corr_test = true;
      break;
    case 4:
      double_var_D_ang_corr_test = true;
      break;
    case 5:
      triple_var_aAB_ang_corr_test = true;
      break;
    default:
      std::cout << "Parameter does not match keyword of any test" << std::endl;
      return 0;
    }
  }
  else if (testName == "sampling"){
    vonNeumann_reject_sampling_test = true;
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

  if (single_var_ang_corr_test){
    CRADLE::test::RunSingleVarAngCorrTest();
  }

  if (double_var_B_ang_corr_test){
    CRADLE::test::RunDoubleVarBAngCorrTest();
  }
  
  if (double_var_A_ang_corr_test){
    CRADLE::test::RunDoubleVarAAngCorrTest();
  }
  
  if (double_var_D_ang_corr_test){
    CRADLE::test::RunDoubleVarDAngCorrTest();
  }

  if (triple_var_aAB_ang_corr_test){
    CRADLE::test::RunTripleVarAngCorrTest();
  }

  if (vonNeumann_reject_sampling_test){
    if (genParam.empty()){
      std::cout << "a = " << a << "\nA = " << A << "\nB = " << B << "\nD = " << D << "\nE = " << E << "\nNumber of decays = " << nDecays << std::endl;
      CRADLE::test::RunDefaultVNRSampling(a, A, B, D, E, nDecays, oFile);
      std::cout << "Results saved in " << oFile << std::endl; 
    }
  }
  
  return 0;
}  
