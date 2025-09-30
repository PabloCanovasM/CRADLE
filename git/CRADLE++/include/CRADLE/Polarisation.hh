#ifndef POLARISATION
#define POLARISATION

#include <math.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <map>
#include <limits>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <complex>
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_sf_result.h"
#include "gsl/gsl_complex_math.h"
#include "gsl/gsl_sf_dilog.h"

#include "CRADLE/Utilities.hh"
/*#include "CRADLE/Screening.hh"
#include "CRADLE/DecayManager.hh"
#include "CRADLE/ConfigParser.hh"
#include "CRADLE/DecayMode.hh"
#include "CRADLE/Particle.hh"*/


namespace CRADLE {

namespace polarisation {
  using namespace boost::numeric::ublas;

  const double PI = 3.14159265359;
  const double C = 299792458;//m/s
  const double EMASSC2 = 510.9989461;//keV
  const double PMASSC2 = 938272.046;//keV
  const double NMASSC2 = 939565.4133;//keV
  const double UMASSC2 = 931494.10242;//keV
  const double ALPHAMASSC2 = 3727379.4066;//keV
  const double FINESTRUCTUREMASSC2 = 3727379.508;//keV
  const double FINESTRUCTURE = 0.0072973525664;
  const double E = 2.718281828459045;
  const double HBAR = 6.58211889e-16;//ev*s
  const double NATURALLENGTH = HBAR*C/EMASSC2/1000.;//m
  const double a_CORR = -1.0;
  const double EULER_MASCHERONI_CONSTANT = 0.577215664901532;  /**< the Euler-Mascheroni constant */

  enum DecayType { FERMI, GAMOW_TELLER, MIXED };

  const std::string atoms[] = {"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt"};


  //Polarisation Variables

  inline double lambdaJiJf_factor(double j_i, double j_f){
    if (j_f == (j_i - 1)){
      return 1;
    }
    else if (j_f == j_i){
      return 1/(j_i+1);
    }
    else if (j_f == (j_i + 1)){
      return -j_i/j_f;
    }
    else
      return 0;
  }

  inline double CalculateBetaAssymetry(std::complex<double> cs, std::complex<double> csp, std::complex<double> ct, std::complex<double> ctp, 
  std::complex<double> cv, std::complex<double> cvp, std::complex<double> ca, std::complex<double> cap, double mf, double mgt, double j_i, double j_f, int betaType, int Z, double energy){
    double coulombCorr = FINESTRUCTURE*Z/std::sqrt(energy*energy/EMASSC2/EMASSC2-1);
    double gamma_ratio = std::sqrt(1. - std::pow(FINESTRUCTURE * Z, 2.))*EMASSC2/energy; //gamma*m_e/E
    double A = mgt*mgt*lambdaJiJf_factor(j_i,j_f)*(2.*betaType*(ct*conj(ctp)-ca*conj(cap)).real()+2.*coulombCorr*(ct*conj(cap)+ctp*conj(ca)).imag()); 
    if (j_i == j_f){
      A += mf*mgt*std::sqrt(j_i/(j_i+1))*(2.*(cs*conj(ctp)+csp*conj(ct)-cv*conj(cap)-cvp*conj(ca)).real()+
                                          2*betaType*coulombCorr*(cs*conj(cap)+csp*conj(ca)-cv*conj(ctp)-cvp*conj(ct)).imag()); 
      }
    A /= utilities::CalculateXiBetaDecay(cs, csp, ct, ctp, cv, cvp, ca, cap, mf, mgt);
    return A;
  }

  inline double CalculateNeutrinoAssymetry(std::complex<double> cs, std::complex<double> csp, std::complex<double> ct, std::complex<double> ctp, 
  std::complex<double> cv, std::complex<double> cvp, std::complex<double> ca, std::complex<double> cap, double mf, double mgt, double j_i, double j_f, int betaType, int Z, double energy){
    double gamma_ratio = std::sqrt(1. - std::pow(FINESTRUCTURE * Z, 2.))*EMASSC2/energy; //gamma*m_e/E
    double B = 2.*mgt*mgt*lambdaJiJf_factor(j_i,j_f)*(gamma_ratio*(ct*conj(cap)+ctp*conj(ca)).real()+betaType*(ct*conj(ctp)+ca*conj(cap)).real());
    if (j_i == j_f){
      B -= mf*mgt*std::sqrt(j_i/(j_i+1))*(2.*(cs*conj(ctp)+csp*conj(ct)+cv*conj(cap)+cvp*conj(ca)).real()+
                                              2.*betaType*gamma_ratio*(cs*conj(cap)+csp*conj(ca)+cv*conj(ctp)+cvp*conj(ct)).real()); 
    }
    B /= utilities::CalculateXiBetaDecay(cs, csp, ct, ctp, cv, cvp, ca, cap, mf, mgt);
    return B;
  }

  inline double CalculateDTripleCorrelation(std::complex<double> cs, std::complex<double> csp, std::complex<double> ct, std::complex<double> ctp, 
  std::complex<double> cv, std::complex<double> cvp, std::complex<double> ca, std::complex<double> cap, double mf, double mgt, double j_i, double j_f, int betaType, int Z, double energy){
    if (j_i == j_f){
    double coulombCorr = FINESTRUCTURE*Z/std::sqrt(energy*energy/EMASSC2/EMASSC2-1);
    double D = mf*mgt*std::sqrt(j_i/(j_i+1))*(-2.*coulombCorr*(cs*conj(ca)+csp*conj(cap)-cv*conj(ct)-cvp*conj(ctp)).real()+
                                              2.*(cs*conj(ct)+csp*conj(ctp)-cv*conj(ca)-cvp*conj(cap)).imag());
    return D/utilities::CalculateXiBetaDecay(cs, csp, ct, ctp, cv, cvp, ca, cap, mf, mgt);
    }
    else
      return 0;
  }

  //Maximum Functions 

  //A: iterate over the energy due to non-zero coulomb Corrections
  inline double* maximumA(std::complex<double> cs, std::complex<double> csp, std::complex<double> ct, std::complex<double> ctp, 
    std::complex<double> cv, std::complex<double> cvp, std::complex<double> ca, std::complex<double> cap, double mf, double mgt, double j_i, double j_f, int betaType, int Z, double E0)  {
    double* max_A_point = new double[2]; // (energy of the maximum, value of the maximum for the abs)
    max_A_point[1] = -1.0; //since we look at the abs, first one will trigger; 
    double kin_en = 0.1; //avoid dividing by zero
    while (kin_en < E0+0.1){
      double en = kin_en + EMASSC2;
      double A = CalculateBetaAssymetry(cs, csp, ct, ctp, cv, cvp, ca, cap, mf, mgt, j_i, j_f, betaType, Z, en);
      if(std::abs(A) > max_A_point[1]){
        //std::cout << "New maximum at " << en << "keV. Previous B: " << max_B << " ; new B: " << B << std::endl;
        max_A_point[0] = en;
        max_A_point[1] = std::abs(A);
        }
      kin_en += 1;
      }
    return max_A_point;
  }

  //B: iterate over the energy
  inline double* maximumB(std::complex<double> cs, std::complex<double> csp, std::complex<double> ct, std::complex<double> ctp, 
  std::complex<double> cv, std::complex<double> cvp, std::complex<double> ca, std::complex<double> cap, double mf, double mgt, double j_i, double j_f, int betaType, int Z, double E0)  {
    double* max_B_point = new double[2]; // (energy of the maximum, value of the maximum for the abs)
    max_B_point[1] = -1.0; //since we look at the abs
    double kin_en = 0.1; //avoid dividing by zero
    while (kin_en < E0+0.1){
      double en = kin_en + EMASSC2;
      double B = CalculateNeutrinoAssymetry(cs, csp, ct, ctp, cv, cvp, ca, cap, mf, mgt, j_i, j_f, betaType, Z, en);
      if(std::abs(B) > max_B_point[1]){
        //std::cout << "New maximum at " << en << "keV. Previous B: " << max_B << " ; new B: " << B << std::endl;
        max_B_point[0] = en;
        max_B_point[1] = std::abs(B);
        }
      kin_en += 1;
      }
    return max_B_point;
  }

  //D: iterate over the energy
  inline double* maximumD(std::complex<double> cs, std::complex<double> csp, std::complex<double> ct, std::complex<double> ctp, 
  std::complex<double> cv, std::complex<double> cvp, std::complex<double> ca, std::complex<double> cap, double mf, double mgt, double j_i, double j_f, int betaType, int Z, double E0)  {
     double* max_D_point = new double[2];  // (energy of the maximum, value of the maximum)
    max_D_point[1] = -1.0;
    double kin_en = 0.1; //avoid dividing by zero
    while (kin_en < (E0 + 0.1)){
      double en = kin_en + EMASSC2;
      double D = CalculateDTripleCorrelation(cs, csp, ct, ctp, cv, cvp, ca, cap, mf, mgt, j_i, j_f, betaType, Z, en);
      if(std::abs(D) > max_D_point[1]){
        //std::cout << "New maximum at " << en << "keV. Previous B: " << max_B << " ; new B: " << B << std::endl;
        max_D_point[0] = en;
        max_D_point[1] = std::abs(D);
        }
      kin_en += 1;
      }
    return max_D_point;
  }

}//closing polarisation namespace
}//closing CRADLE namespace
#endif
