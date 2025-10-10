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

  //Polarisation Variables (from Nuclear Physics 4 (1957) 206-212; by J. D. Jackson)

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
    double coulombCorr = FINESTRUCTURE*Z/std::sqrt(energy*energy/EMASSC2/EMASSC2-1); //alpha*Z*m_e/p_e
    double gamma_ratio = std::sqrt(1. - std::pow(FINESTRUCTURE * Z, 2.))*EMASSC2/energy; //gamma*m_e/E
    double A = mgt*mgt*lambdaJiJf_factor(j_i,j_f)*(2.*betaType*(ct*conj(ctp)-ca*conj(cap)).real()+2.*coulombCorr*(ct*conj(cap)+ctp*conj(ca)).imag()); 
    if (j_i == j_f){ //Mixed decay factors
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
    if (j_i == j_f){ //Mixed decay factors
      B -= mf*mgt*std::sqrt(j_i/(j_i+1))*(2.*(cs*conj(ctp)+csp*conj(ct)+cv*conj(cap)+cvp*conj(ca)).real()+
                                              2.*betaType*gamma_ratio*(cs*conj(cap)+csp*conj(ca)+cv*conj(ctp)+cvp*conj(ct)).real()); 
    }
    B /= utilities::CalculateXiBetaDecay(cs, csp, ct, ctp, cv, cvp, ca, cap, mf, mgt);
    return B;
  }

  inline double CalculateDTripleCorrelation(std::complex<double> cs, std::complex<double> csp, std::complex<double> ct, std::complex<double> ctp, 
  std::complex<double> cv, std::complex<double> cvp, std::complex<double> ca, std::complex<double> cap, double mf, double mgt, double j_i, double j_f, int betaType, int Z, double energy){
    //Only mixed decay factors
    if (j_i == j_f){ 
    double coulombCorr = FINESTRUCTURE*Z/std::sqrt(energy*energy/EMASSC2/EMASSC2-1);
    double D = mf*mgt*std::sqrt(j_i/(j_i+1))*(-2.*coulombCorr*(cs*conj(ca)+csp*conj(cap)-cv*conj(ct)-cvp*conj(ctp)).real()+
                                              2.*(cs*conj(ct)+csp*conj(ctp)-cv*conj(ca)-cvp*conj(cap)).imag());
    return D/utilities::CalculateXiBetaDecay(cs, csp, ct, ctp, cv, cvp, ca, cap, mf, mgt);
    }
    else
      return 0;
  }

  //Maximum Functions for A, B, D

  //A: iterate over the energy due to non-zero coulomb Corrections, dependnt on energy
  inline double* MaximumA(std::complex<double> cs, std::complex<double> csp, std::complex<double> ct, std::complex<double> ctp, 
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
  inline double* MaximumB(std::complex<double> cs, std::complex<double> csp, std::complex<double> ct, std::complex<double> ctp, 
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
  inline double* MaximumD(std::complex<double> cs, std::complex<double> csp, std::complex<double> ct, std::complex<double> ctp, 
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

 
  inline double CalculateAngularCorrelationFactor(double a, double A, double B, double D, double E, double cosTheta_e, double cosTheta_enu, double phi){
  /*Computation of the angular dependent factor (ie proportional to xi) in formula 1 from the Jackson 1957 paper referenced above*/
    ublas::vector<double> elDir (3);
    ublas::vector<double> enuDir (3);
    
    //electron in XZ plane, using axial symmetry in Z direction (direction of J)
    double sinTheta_e = std::sqrt(1-cosTheta_e*cosTheta_e);
    elDir(0) = sinTheta_e;
    elDir(1) = 0; 
    elDir(2) = cosTheta_e;

    double sinTheta_enu = std::sqrt(1-cosTheta_enu*cosTheta_enu);
    enuDir(0) = sinTheta_enu*std::cos(phi);
    enuDir(1) = sinTheta_enu*std::sin(phi);
    enuDir(2) = cosTheta_enu;

    double beta_e = std::sqrt(1-EMASSC2*EMASSC2/E/E); // p_e/E_e; p_nu/E_nu = 1
    
    double angCorrFactor = 1.;
    angCorrFactor += a*beta_e*(elDir(0)*enuDir(0)+elDir(1)*enuDir(1)+elDir(2)*enuDir(2));
    angCorrFactor += A*beta_e*elDir(2);
    angCorrFactor += B*enuDir(2);
    angCorrFactor += D*beta_e*(elDir(0)*enuDir(1)-elDir(1)*enuDir(0));

    return angCorrFactor;
  }

  inline double MaximumAngCorrFactor(double a, double A, double B, double D, double E){
    /*Search of the maximum using grid search*/
    double maxAngCorrFactor =  std::numeric_limits<double>::lowest(); //not sure if value is nonpositive
    for (int z_e = -200; z_e <= 200; z_e++){
      for (int z_enu = -200; z_enu <= 200; z_enu++){
	for (int phi = 0; phi < 720; phi++){
	  double cosTheta_e = ((double) z_e)/200;
	  double cosTheta_enu = ((double) z_enu)/200;
	  double radPhi = phi*PI/360;
	  double angCorrFactor = CalculateAngularCorrelationFactor(a, A, B, D, E, cosTheta_e, cosTheta_enu, radPhi);
	  if (angCorrFactor > maxAngCorrFactor)
	    maxAngCorrFactor = angCorrFactor;
	  }
	}
      if (z_e%10 == 0) std::cout << "cos(theta_e) = "<< ((double) z_e)/200 << std::endl;
    }
    return maxAngCorrFactor;
  }
  
}//closing polarisation namespace
}//closing CRADLE namespace
#endif
