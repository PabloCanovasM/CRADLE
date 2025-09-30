#include "CRADLE/DecayMode.hh"
#include "CRADLE/DecayManager.hh"
#include "CRADLE/Particle.hh"
#include "CRADLE/Utilities.hh"
#include "CRADLE/SpectrumGenerator.hh"
#include "CRADLE/RadiativeCorrections.hh"

#include <string>
#include <sstream>
#include <complex>

namespace CRADLE {

void DecayMode::FourBodyDecay(ublas::vector<double>& velocity, Particle* finalState1, Particle* finalState2, Particle* finalState3, Particle* finalState4) //, ublas::vector<double>& dir1, ublas::vector<double>& dir2, ublas::vector<double>& dirg, double EnergyElectron, double EnergyNeutrino, double EnergyBrPhoton) 
{
  ublas::vector<double> momentum1 (4);
  ublas::vector<double> momentum2 (4);
  ublas::vector<double> momentum3 (4);
  ublas::vector<double> momentumg (4);

  double mass1 = finalState1->GetMass();
  double mass2 = finalState2->GetMass();
  double mass3 = finalState3->GetMass();
  double massg = finalState4->GetMass();

  ublas::vector<double> p1 = finalState1->Get3Momentum();
  ublas::vector<double> p2 = finalState2->Get3Momentum();
  ublas::vector<double> pg = finalState4->Get3Momentum();
  ublas::vector<double> p3 = -(p1+p2+pg) ;

  double p1Norm = utilities::GetNorm(p1);
  double p2Norm = utilities::GetNorm(p2);
  double p3Norm = utilities::GetNorm(p3);
  double pgNorm = utilities::GetNorm(pg);

  momentum1(0) = std::sqrt(mass1*mass1 + std::pow(p1Norm, 2)) ;
  momentum1(1) = p1[0] ;
  momentum1(2) = p1[1] ;
  momentum1(3) = p1[2] ;

  momentum2(0) = std::sqrt(mass2*mass2 + std::pow(p2Norm, 2)) ;
  momentum2(1) = p2[0] ;
  momentum2(2) = p2[1] ;
  momentum2(3) = p2[2] ;

  momentumg(0) = std::sqrt(massg*massg + std::pow(pgNorm, 2)) ;
  momentumg(1) = pg[0] ;
  momentumg(2) = pg[1] ;
  momentumg(3) = pg[2] ;

  momentum3(0) = std::sqrt(mass3*mass3 + std::pow(p3Norm, 2)) ;
  momentum3(1) = p3[0] ;
  momentum3(2) = p3[1] ;
  momentum3(3) = p3[2] ;

  /*
  ublas::vector<double> p1 = EnergyElectron * dir1 ;  
  ublas::vector<double> p2 = EnergyNeutrino * dir2 ;
  ublas::vector<double> pg =  EnergyBrPhoton * dirg ;
  ublas::vector<double> p3 = -(p1+p2+pg) ;
  double p3Norm = utilities::GetNorm(p3);

  momentum1(0) = std::sqrt(mass1*mass1 + EnergyElectron*EnergyElectron) ;
  momentum1(1) = p1[0] ;
  momentum1(2) = p1[1] ;
  momentum1(3) = p1[2] ;

  momentum2(0) = std::sqrt(mass2*mass2 + EnergyNeutrino*EnergyNeutrino) ;
  momentum2(1) = p2[0] ;
  momentum2(2) = p2[1] ;
  momentum2(3) = p2[2] ;

  momentumg(0) = std::sqrt(massg*massg + EnergyBrPhoton*EnergyBrPhoton) ;
  momentumg(1) = pg[0] ;
  momentumg(2) = pg[1] ;
  momentumg(3) = pg[2] ;

  momentum3(0) = std::sqrt(mass3*mass3 + p3Norm*p3Norm) ;
  momentum3(1) = p3[0] ;
  momentum3(2) = p3[1] ;
  momentum3(3) = p3[2] ;
  */

  finalState1->SetMomentum(utilities::LorentzBoost(velocity, momentum1));
  finalState2->SetMomentum(utilities::LorentzBoost(velocity, momentum2));
  finalState3->SetMomentum(utilities::LorentzBoost(velocity, momentum3));
  finalState4->SetMomentum(utilities::LorentzBoost(velocity, momentumg));
}

void DecayMode::ThreeBodyDecay(ublas::vector<double>& velocity, Particle* finalState1, Particle* finalState2, Particle* finalState3, ublas::vector<double>& dir2, double Q) {
  //Perform decay in CoM frame
  ublas::vector<double> momentum1 = finalState1->GetMomentum();
  ublas::vector<double> momentum2 (4);
  ublas::vector<double> momentum3 (4);

  ublas::vector<double> p2 (3);
  double p2Norm = 0.;

  double mass1 = finalState1->GetMass();
  double mass2 = finalState2->GetMass();
  double mass3 = finalState3->GetMass();
  ublas::vector<double> p1 = finalState1->Get3Momentum();

  double a = mass2*mass2;
  double b = mass3*mass3;
  double c = utilities::GetNorm(p1);
  double d = Q + mass1 + mass2 + mass3 - momentum1(0);
  double e = inner_prod(p1, dir2)/c;

  double first = 1./2./(c*c*e*e-d*d);
  double second = a*a*d*d-2*a*b*d*d+4.*a*c*c*d*d*e*e-2.*a*c*c*d*d-2.*a*d*d*d*d+b*b*d*d+2*b*c*c*d*d-2.*b*d*d*d*d+c*c*c*c*d*d-2.*c*c*d*d*d*d+d*d*d*d*d*d;
  double third = a*c*e-b*c*e-c*c*c*e+c*d*d*e;

  p2Norm = first*(-std::sqrt(second)+third);
  p2 = p2Norm*dir2;

  ublas::vector<double> p3 = -(p1+p2);
  double p3Norm = utilities::GetNorm(p3);

  momentum2(0) = std::sqrt(a+p2Norm*p2Norm);
  momentum2(1) = p2(0);
  momentum2(2) = p2(1);
  momentum2(3) = p2(2);

  momentum3(0) = std::sqrt(b+p3Norm*p3Norm);
  momentum3(1) = p3(0);
  momentum3(2) = p3(1);
  momentum3(3) = p3(2);

  //std::cout << "\t" << inner_prod(p1, p2)/p2Norm/c << std::endl;

  // Perform Lorentz boost back to lab frame
  finalState1->SetMomentum(utilities::LorentzBoost(velocity, momentum1));
  finalState2->SetMomentum(utilities::LorentzBoost(velocity, momentum2));
  finalState3->SetMomentum(utilities::LorentzBoost(velocity, momentum3));

}

void DecayMode::TwoBodyDecay(ublas::vector<double>& velocity, Particle* finalState1, Particle* finalState2, double Q) {
  ublas::vector<double> momentum1 (4);
  ublas::vector<double> momentum2 (4);

  ublas::vector<double> dir = utilities::RandomDirection();

  double mass1 = finalState1->GetMass();
  double mass2 = finalState2->GetMass();
  
  double M = Q + mass1 + mass2;

  double p = 1./(2.*M)*std::sqrt((M*M-std::pow(mass1-mass2, 2.))*(M*M-std::pow(mass1+mass2, 2.)));
  
  double p1 = 1./(2.*M) * std::sqrt( std::pow(M*M + mass1*mass1 - mass2*mass2, 2) - 4*M*M*mass1*mass1 );
  double p2 = 1./(2.*M) * std::sqrt( std::pow(M*M + mass2*mass2 - mass1*mass1, 2) - 4*M*M*mass2*mass2 );

  double energy1 = std::sqrt(mass1*mass1+p1*p1);
  double energy2 = std::sqrt(mass2*mass2+p2*p2);

  
  momentum1(0) = energy1;
  momentum1(1) = p1*dir[0];
  momentum1(2) = p1*dir[1];
  momentum1(3) = p1*dir[2];

  momentum2(0) = energy2;
  momentum2(1) = -p2*dir[0];
  momentum2(2) = -p2*dir[1];
  momentum2(3) = -p2*dir[2];
  
  /*std::cout << "mass li : " << mass1 << "\n" ;
  std::cout << "mass gamma : " << mass2 << "\n";
  std::cout << "Q : " << Q << "\n" ;
  std::cout << "M : " << M << "\n";
  std::cout << "p : " << p << "\n";
  std::cout << "p1 : " << p1 << "\n";
  std::cout << "p2 : " << p2 << "\n" ;
  std::cout << "energy 1 : " << energy1 << "\n";
  std::cout << "energy 2 : " << energy2 << "\n";*/
  ;

  // Perform Lorentz boost back to lab frame
  finalState1->SetMomentum(utilities::LorentzBoost(velocity, momentum1));
  finalState2->SetMomentum(utilities::LorentzBoost(velocity, momentum2));
  
}

std::vector<Particle*> BetaMinus::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates;

  // std::cout << "In BetaMinus Decay " << std::endl;
  // std::cout << "Address: " << initState << std::endl;
  std::ostringstream oss;
  oss << initState->GetCharge()+initState->GetNeutrons() << utilities::atoms[initState->GetCharge()];
  //std::cout << oss.str() << std::endl;
  Particle* recoil = DecayManager::GetInstance().GetNewParticle(oss.str(), initState->GetCharge()+1, initState->GetCharge()+initState->GetNeutrons());
  recoil->SetExcitationEnergy(daughterExEn);
    
  Particle* e = DecayManager::GetInstance().GetNewParticle("e-");
  Particle* enu = DecayManager::GetInstance().GetNewParticle("enubar");

  // std::cout << "Recoil " << recoil->GetCharge() << " " << recoil->GetNeutrons() << " " << recoil << std::endl;

  oss.str("");
  oss.clear();
  oss << "BetaMinus:Z" << recoil->GetCharge() << "A" << recoil->GetCharge() + recoil->GetNeutrons() << "Q" << Q;
  //Work in the COM frame
  ublas::vector<double> elFourMomentum (4);

  DecayManager& dm = DecayManager::GetInstance();
  double CS = dm.configOptions.couplingConstants.CS.real();
  double CSP = dm.configOptions.couplingConstants.CSP.real();
  double CV = dm.configOptions.couplingConstants.CV.real();
  double CVP = dm.configOptions.couplingConstants.CVP.real();
  double CA = dm.configOptions.couplingConstants.CA.real();
  double CAP = dm.configOptions.couplingConstants.CAP.real();
  double CT = dm.configOptions.couplingConstants.CT.real();
  double CTP = dm.configOptions.couplingConstants.CTP.real();
  double a_conf = dm.configOptions.couplingConstants.a;
  double b_conf = dm.configOptions.couplingConstants.b;

  double mf = 0.;
  double mgt = 0.;
  
  if (dm.configOptions.betaDecay.Default == "Fermi") {
    mf = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Gamow-Teller") {
    mgt = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Mixed") {
    mf = 1. ;
    mgt = 1. ;
  }
  else if (dm.configOptions.betaDecay.Default == "Auto") {
    std::string Type;
    try{
      Type = DecayManager::GetInstance().GetBetaType(oss.str());
    } catch (const std::invalid_argument& e) {
      Type = utilities::FindBetaType(initState, recoil);
      DecayManager::GetInstance().RegisterBetaType(oss.str(), Type);
    }
    if (Type == "Fermi") {
      mf = 1.;
    } else if (Type == "Gamow-Teller") {
      mgt = 1.;
    } else {
      mgt = 1. ;
      mf = 1. ;
    }
  } 
  //mgt = 0. ;
  //mf = 1. ;
  double a = utilities::CalculateBetaNeutrinoAsymmetry(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf);
  double fierz = utilities::CalculateFierz(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf);

  //std::cout << "fierz " << fierz << " a " << a << std::endl;
  
  std::vector<std::vector<double> >* dist;
  try {
    dist = DecayManager::GetInstance().GetDistribution(oss.str());
  } catch (const std::invalid_argument& e) {
    bool advancedFermi = false;
    if (dm.configOptions.betaDecay.FermiFunction == "Advanced") { //// réécriture de la condition par SL 10/05/2023
      advancedFermi = true;
    }
    dist = spectrumGen->GenerateSpectrum(initState, recoil, Q); //// changement de Q en E0 par SL 10/05/2023
    double gamma = std::sqrt(1-std::pow(utilities::FINESTRUCTURE*recoil->GetCharge(), 2.));
    for ( int i = 0; i<dist->size(); i++) {                                  ////// changement de boucle et element par ((*dist)[i]) par SL 10/05/2023
      double E = ((*dist)[i])[0]+utilities::EMASSC2;
      double SH = ((*dist)[i])[1];
      ((*dist)[i])[1] = SH*(1+gamma*fierz*utilities::EMASSC2/E);
    }
    DecayManager::GetInstance().RegisterDistribution(oss.str(), dist);
    }

  //std::cout << "Found distribution" << std::endl;

  std::vector<double> p;

  double elEnergy = utilities::RandomFromDistribution(*dist)+utilities::EMASSC2;
  double elMomentum = std::sqrt(elEnergy*elEnergy-std::pow(utilities::EMASSC2, 2.));
  ublas::vector<double> enuDir = utilities::RandomDirection();

  p.push_back(1.);
  p.push_back(a*elMomentum/elEnergy);
  ublas::vector<double> eDir = utilities::GetParticleDirection(enuDir, p);

  elFourMomentum(0) = elEnergy;
  elFourMomentum(1) = elMomentum*eDir[0];
  elFourMomentum(2) = elMomentum*eDir[1];
  elFourMomentum(3) = elMomentum*eDir[2];

  e->SetMomentum(elFourMomentum);

  ublas::vector<double> velocity = -initState->GetVelocity();
  ThreeBodyDecay(velocity, e, enu, recoil, enuDir, Q);

  finalStates.push_back(recoil);
  finalStates.push_back(e);
  finalStates.push_back(enu);

  return finalStates;
}

std::vector<Particle*> BetaPlus::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates;

  DecayManager& dm = DecayManager::GetInstance();

  double E0 = Q-2*utilities::EMASSC2;

  //std::cout << "In BetaPlus Decay " << std::endl;
  //std::cout << "Address: " << initState << std::endl;
  std::ostringstream oss;
  oss << initState->GetCharge()+initState->GetNeutrons() << utilities::atoms[initState->GetCharge()-2];
  Particle* recoil = DecayManager::GetInstance().GetNewParticle(oss.str(), initState->GetCharge()-1, initState->GetCharge()+initState->GetNeutrons());
  recoil->SetExcitationEnergy(daughterExEn);
  Particle* pos = DecayManager::GetInstance().GetNewParticle("e+");
  Particle* enubar = DecayManager::GetInstance().GetNewParticle("enu");

  oss.str("");
  oss.clear();

  oss << "BetaPlus:Z" << recoil->GetCharge() << "A" << recoil->GetCharge() + recoil->GetNeutrons() << "Q" << Q;
  //Work in the COM frame
  ublas::vector<double> enubarDir = utilities::RandomDirection();

  ublas::vector<double> posFourMomentum (4);

  double mf = 0.;
  double mgt = 0.;

  /////////ajout de SL 12/05/2023//////////////

  if (dm.configOptions.betaDecay.Default == "Fermi") {
    mf = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Gamow-Teller") {
    mgt = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Mixed") {
    mf = 1. ;
    mgt = 1. ;
  }
  else if (dm.configOptions.betaDecay.Default == "Auto") {
    std::string Type;
    try{
      Type = DecayManager::GetInstance().GetBetaType(oss.str());
    } catch (const std::invalid_argument& e) {
      Type = utilities::FindBetaType(initState, recoil);
      DecayManager::GetInstance().RegisterBetaType(oss.str(), Type);
    }
    if (Type == "Fermi") {
      mf = 1.;
    } else if (Type == "Gamow-Teller") {
      mgt = 1.;
    } else {
      mgt = 1. ;
      mf = 1. ;
    }
  }
  //std::cout << "mf : " << mf << "\n";
  //std::cout << "mgt : " << mgt << "\n";
  //mf = 1.;
  //mgt = 0.;
  ////////////////////////////////////////////////

  double CS = dm.configOptions.couplingConstants.CS.real();
  double CSP = dm.configOptions.couplingConstants.CSP.real();
  double CV = dm.configOptions.couplingConstants.CV.real();
  double CVP = dm.configOptions.couplingConstants.CVP.real();
  double CA = dm.configOptions.couplingConstants.CA.real();
  double CAP = dm.configOptions.couplingConstants.CAP.real();
  double CT = dm.configOptions.couplingConstants.CT.real();
  double CTP = dm.configOptions.couplingConstants.CTP.real();
  double a_conf = dm.configOptions.couplingConstants.a;
  double b_conf = dm.configOptions.couplingConstants.b;

  double a = utilities::CalculateBetaNeutrinoAsymmetry(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf);
  double fierz = utilities::CalculateFierz(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf);

  //std::cout <<" b = " << fierz <<"\t a = " << a << std::endl;
  //std::cout << "a : " << a << "\n";

  std::vector<std::vector<double> >* dist;
  try {
    dist = DecayManager::GetInstance().GetDistribution(oss.str());
  } catch (const std::invalid_argument& e) {
    bool advancedFermi = false;
    if (dm.configOptions.betaDecay.FermiFunction == "Advanced") { //// réécriture de la condition par SL 10/05/2023
      advancedFermi = true;
    }
    dist = spectrumGen->GenerateSpectrum(initState, recoil, E0); //// changement de Q en E0 par SL 10/05/2023
    double gamma = std::sqrt(1-std::pow(utilities::FINESTRUCTURE*recoil->GetCharge(), 2.));
    int i=0;
    for ( int i = 0; i<dist->size(); i++) {                                  ////// changement de boucle et element par ((*dist)[i]) par SL 10/05/2023
      double E = ((*dist)[i])[0]+utilities::EMASSC2;
      double SH = ((*dist)[i])[1];
      ((*dist)[i])[1] = SH*(1+gamma*fierz*utilities::EMASSC2/E);
    }
    DecayManager::GetInstance().RegisterDistribution(oss.str(), dist);
    }


  // ublas::vector<std::vector<double> >* dist;
  double posEnergy = utilities::RandomFromDistribution(*dist) + utilities::EMASSC2;
  double posMomentum = std::sqrt(posEnergy*posEnergy-std::pow(utilities::EMASSC2, 2.));

  std::vector<double> p;
  p.push_back(1.);
  p.push_back(a*posMomentum/posEnergy);
  ublas::vector<double> posDir = utilities::GetParticleDirection(enubarDir, p);
  posFourMomentum(0) = posEnergy;
  posFourMomentum(1) = posMomentum*posDir[0];
  posFourMomentum(2) = posMomentum*posDir[1];
  posFourMomentum(3) = posMomentum*posDir[2];

  pos->SetMomentum(posFourMomentum);

  ublas::vector<double> velocity = -initState->GetVelocity();
  ThreeBodyDecay(velocity, pos, enubar, recoil, enubarDir, E0);


  finalStates.push_back(recoil);
  finalStates.push_back(pos);
  finalStates.push_back(enubar);

  return finalStates;
}

////////////////////////////////////////////////////////////
std::vector<Particle*> BetaMinusRadiative::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates;

  // std::cout << "In BetaMinusRadiativeCorrection Decay " << std::endl;
  // std::cout << "Address: " << initState << std::endl;
  std::ostringstream oss;
  oss << initState->GetCharge()+initState->GetNeutrons() << utilities::atoms[initState->GetCharge()];
  //std::cout << oss.str() << std::endl;
  Particle* recoil = DecayManager::GetInstance().GetNewParticle(oss.str()+"RC", initState->GetCharge()+1, initState->GetCharge()+initState->GetNeutrons());
  recoil->SetExcitationEnergy(daughterExEn);
  
  //std::string parentName = oss.str();
  
  Particle* e = DecayManager::GetInstance().GetNewParticle("e-RC");
  Particle* enubar = DecayManager::GetInstance().GetNewParticle("enubarRC");
  Particle* gamma = DecayManager::GetInstance().GetNewParticle("gammaRC");


  // std::cout << "Recoil " << recoil->GetCharge() << " " << recoil->GetNeutrons() << " " << recoil << std::endl;

  oss.str("");
  oss.clear();
  oss << "BetaMinusRadiative:Z" << recoil->GetCharge() << "A" << recoil->GetCharge() + recoil->GetNeutrons() << "Q" << Q;
  //Work in the COM frame
  ublas::vector<double> elFourMomentum (4);

  DecayManager& dm = DecayManager::GetInstance();
  double CS = dm.configOptions.couplingConstants.CS.real();
  double CSP = dm.configOptions.couplingConstants.CSP.real();
  double CV = dm.configOptions.couplingConstants.CV.real();
  double CVP = dm.configOptions.couplingConstants.CVP.real();
  double CA = dm.configOptions.couplingConstants.CA.real();
  double CAP = dm.configOptions.couplingConstants.CAP.real();
  double CT = dm.configOptions.couplingConstants.CT.real();
  double CTP = dm.configOptions.couplingConstants.CTP.real();
  double a_conf = dm.configOptions.couplingConstants.a;
  double b_conf = dm.configOptions.couplingConstants.b;

  ublas::vector<double> eFourMomentum (4);
  ublas::vector<double> enubarFourMomentum (4);
  ublas::vector<double> gammaFourMomentum (4);

  double mf = 0.;
  double mgt = 0.;

  if (dm.configOptions.betaDecay.Default == "Fermi") {
    mf = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Gamow-Teller") {
    mgt = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Mixed") {
    mf = 1. ;
    mgt = 1. ;
  }
  else if (dm.configOptions.betaDecay.Default == "Auto") {
    std::string Type;
    try{
      Type = DecayManager::GetInstance().GetBetaType(oss.str());
    } catch (const std::invalid_argument& e) {
      Type = utilities::FindBetaType(initState, recoil);
      DecayManager::GetInstance().RegisterBetaType(oss.str(), Type);
    }
    if (Type == "Fermi") {
      mf = 1.;
    } else if (Type == "Gamow-Teller") {
      mgt = 1.;
    } else {
      mgt = 1. ;
      mf = 1. ;
    }
  } 
  //mgt = 0. ;
  //mf = 1. ;

  double a = utilities::CalculateBetaNeutrinoAsymmetry(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf);
  double fierz = utilities::CalculateFierz(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf);

  //std::cout << "fierz " << fierz << " a " << a << std::endl;
 
  int nH = 0;
  double mass_i = initState->GetMass();
  double mass_f = mass_i - Q ; //recoil->GetMass() ;
  double Cs = dm.configOptions.betaDecay.OmegaValue ;
  int A = recoil-> GetCharge() + recoil->GetNeutrons() ; 
  double R = utilities::ApproximateRadius(A) ;
  int Z = (recoil->GetCharge() - initState->GetCharge())*recoil->GetCharge();
  int betaType = (int)((Z > 0) - (Z < 0));
  Z = std::abs(Z) ;

  std::vector<double> W;
  try {
    W = DecayManager::GetInstance().GetParameterMC(oss.str()) ;
  } catch (const std::invalid_argument& e) {
    double WHmax = radiativecorrections::WH_max(1000000, Cs, a, mass_i, mass_f, Z, R, betaType, "BetaMinus");
    W = {WHmax, 0} ;
    DecayManager::GetInstance().RegisterParameterMC(oss.str(), W) ;
  }

  std::random_device rd;
  std::mt19937 generator(rd()); 
  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  std::uniform_real_distribution<double> distribution_wHmax(0.0, W[0]); 

  int n_verif = nH ;
  while (n_verif == nH) {
    double wH_NR = distribution_wHmax(generator) ;

    double U[8] = {distribution(generator), distribution(generator), distribution(generator), distribution(generator), distribution(generator), distribution(generator), distribution(generator), distribution(generator)};
    double E2 = 1. + (radiativecorrections::delta(mass_i, mass_f, "BetaMinus") - 1.) * U[0] ;
    double E10 = radiativecorrections::delta(mass_i, mass_f, "BetaMinus") - E2 ;
    double omega = Cs * E10 ;
    double K = omega * exp(-U[1] * log(Cs)) ;
    double E1 = E10 - K ;

    double BETA = std::sqrt(1. - 1./std::pow(E2, 2) ) ;
    double N = 0.5 * log( (1. + BETA)/(1. - BETA) ) ; 

    double COS_GAMMA = ( 1. - (1. + BETA) * exp(-2. * N * U[2]) )/BETA ;
    double COS_NEUTRINO = 2. * U[3] - 1. ;
    double COS_ELECTRON = 2. * U[4] - 1. ;

    double PHI_GAMMA = 2. * utilities::PI * U[5] ;
    double PHI_NEUTRINO = 2. * utilities::PI * U[6] ;
    double PHI_ELECTRON = 2. * utilities::PI * U[7] ;


    double SIN_GAMMA = std::sqrt( (1. - std::pow(COS_GAMMA, 2)) ) ;
    double SIN_NEUTRINO = std::sqrt( (1. - std::pow(COS_NEUTRINO, 2)) ) ;
    double SIN_ELECTRON = std::sqrt( (1. - std::pow(COS_ELECTRON, 2)) ) ;
                
    ublas::vector<double> n_ELECTRON (3); 
    n_ELECTRON[0]=SIN_ELECTRON * cos(PHI_ELECTRON); n_ELECTRON[1]=SIN_ELECTRON * sin(PHI_ELECTRON) ; n_ELECTRON[2]=COS_ELECTRON; 
    double n_ELECTRON_PRIME[3] = {-sin(PHI_ELECTRON), cos(PHI_ELECTRON), 0}; 
    double n_ELECTRON_SECOND[3] = {-COS_ELECTRON * cos(PHI_ELECTRON), -COS_ELECTRON * sin(PHI_ELECTRON), SIN_ELECTRON};

    double n_PERPENDICULAIRE_GAMMA[3] ; 
    ublas::vector<double> n_GAMMA (3) ;
    ublas::vector<double> n_NEUTRINO (3);
    n_NEUTRINO[0]=SIN_NEUTRINO*cos(PHI_NEUTRINO) ; n_NEUTRINO[1]=SIN_NEUTRINO * sin(PHI_NEUTRINO) ; n_NEUTRINO[2]=COS_NEUTRINO;
    for (int j = 0; j < 3; j++) {
      n_PERPENDICULAIRE_GAMMA[j] = n_ELECTRON_PRIME[j] * cos(PHI_GAMMA) + n_ELECTRON_SECOND[j] * sin(PHI_GAMMA) ;
      n_GAMMA[j] = n_ELECTRON[j] * COS_GAMMA + n_PERPENDICULAIRE_GAMMA[j] * SIN_GAMMA ;
    }

    double N1_N2 = n_NEUTRINO[0]*n_ELECTRON[0] + n_NEUTRINO[1]*n_ELECTRON[1] + n_NEUTRINO[2]*n_ELECTRON[2] ;
    double N1_K = n_NEUTRINO[0]*n_GAMMA[0] + n_NEUTRINO[1]*n_GAMMA[1] + n_NEUTRINO[2]*n_GAMMA[2] ;

    if (wH_NR < radiativecorrections::WH(E2, K, COS_GAMMA, N1_K, N1_N2, a, mass_i, mass_f, Z, R, betaType, "BetaMinus")) {
      nH += 1 ;

      ublas::vector<double> velocity = -initState->GetVelocity();
      double eMomentum = std::sqrt(std::pow(E2*utilities::EMASSC2, 2) -std::pow(utilities::EMASSC2, 2.));
      double enubarMomentum = E1*utilities::EMASSC2 ;
      double gammaMomentum = K*utilities::EMASSC2 ;
      
      eFourMomentum(0) = E2*utilities::EMASSC2;
      eFourMomentum(1) = eMomentum*n_ELECTRON[0];
      eFourMomentum(2) = eMomentum*n_ELECTRON[1];
      eFourMomentum(3) = eMomentum*n_ELECTRON[2];

      enubarFourMomentum(0) = enubarMomentum ; 
      enubarFourMomentum(1) = enubarMomentum*n_NEUTRINO[0];
      enubarFourMomentum(2) = enubarMomentum*n_NEUTRINO[1];
      enubarFourMomentum(3) = enubarMomentum*n_NEUTRINO[2];

      gammaFourMomentum(0) = gammaMomentum ;
      gammaFourMomentum(1) = gammaMomentum*n_GAMMA[0];
      gammaFourMomentum(2) = gammaMomentum*n_GAMMA[1];
      gammaFourMomentum(3) = gammaMomentum*n_GAMMA[2];

      e->SetMomentum(eFourMomentum);
      enubar->SetMomentum(enubarFourMomentum);
      gamma->SetMomentum(gammaFourMomentum);
      FourBodyDecay(velocity, e, enubar, recoil, gamma);

      finalStates.push_back(e);
      finalStates.push_back(enubar);
      finalStates.push_back(recoil);
      finalStates.push_back(gamma) ;
    } 
  }
  return finalStates;
}

std::vector<Particle*> BetaMinusVirtualSoft::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates;

  // std::cout << "In BetaMinusRadiativeCorrection Decay " << std::endl;
  // std::cout << "Address: " << initState << std::endl;
  std::ostringstream oss;
  oss << initState->GetCharge()+initState->GetNeutrons() << utilities::atoms[initState->GetCharge()];
  //std::cout << oss.str() << std::endl;
  Particle* recoil = DecayManager::GetInstance().GetNewParticle(oss.str(), initState->GetCharge()+1, initState->GetCharge()+initState->GetNeutrons());
  recoil->SetExcitationEnergy(daughterExEn);
  
  //std::string parentName = oss.str();
  
  Particle* e = DecayManager::GetInstance().GetNewParticle("e-");
  Particle* enubar = DecayManager::GetInstance().GetNewParticle("enubar");

  // std::cout << "Recoil " << recoil->GetCharge() << " " << recoil->GetNeutrons() << " " << recoil << std::endl;
    
  oss.str("");
  oss.clear();
  oss << "BetaMinusVirtualSoft:Z" << recoil->GetCharge() << "A" << recoil->GetCharge() + recoil->GetNeutrons() << "Q" << Q;
  //Work in the COM frame
  ublas::vector<double> elFourMomentum (4);

  DecayManager& dm = DecayManager::GetInstance();
  double CS = dm.configOptions.couplingConstants.CS.real();
  double CSP = dm.configOptions.couplingConstants.CSP.real();
  double CV = dm.configOptions.couplingConstants.CV.real();
  double CVP = dm.configOptions.couplingConstants.CVP.real();
  double CA = dm.configOptions.couplingConstants.CA.real();
  double CAP = dm.configOptions.couplingConstants.CAP.real();
  double CT = dm.configOptions.couplingConstants.CT.real();
  double CTP = dm.configOptions.couplingConstants.CTP.real();
  double a_conf = dm.configOptions.couplingConstants.a;
  double b_conf = dm.configOptions.couplingConstants.b;

  ublas::vector<double> eFourMomentum (4);
  ublas::vector<double> enubarFourMomentum (4);

  double mf = 0.;
  double mgt = 0.;

  if (dm.configOptions.betaDecay.Default == "Fermi") {
    mf = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Gamow-Teller") {
    mgt = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Mixed") {
    mf = 1. ;
    mgt = 1. ;
  }
  else if (dm.configOptions.betaDecay.Default == "Auto") {
    std::string Type;
    try{
      Type = DecayManager::GetInstance().GetBetaType(oss.str());
    } catch (const std::invalid_argument& e) {
      Type = utilities::FindBetaType(initState, recoil);
      DecayManager::GetInstance().RegisterBetaType(oss.str(), Type);
    }
    if (Type == "Fermi") {
      mf = 1.;
    } else if (Type == "Gamow-Teller") {
      mgt = 1.;
    } else {
      mgt = 1. ;
      mf = 1. ;
    }
  } 
  //mgt = 0. ;
  //mf = 1. ;
  
  double a = utilities::CalculateBetaNeutrinoAsymmetry(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf);
  double fierz = utilities::CalculateFierz(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf);
  
  //std::cout << "a : " << a << "\n";
  //std::cout << "fierz " << fierz << " a " << a << std::endl;
 
  int n0VS = 0;
  double mass_i = initState->GetMass();
  double mass_f = mass_i - Q ; //recoil->GetMass() ;
  double Cs = dm.configOptions.betaDecay.OmegaValue ;
  int A = recoil-> GetCharge() + recoil->GetNeutrons() ; 
  double R = utilities::ApproximateRadius(A) ;
  int Z = (recoil->GetCharge() - initState->GetCharge())*recoil->GetCharge();
  int betaType = (int)((Z > 0) - (Z < 0));
  Z = std::abs(Z) ;

  std::vector<double> W;
  try {
    W = DecayManager::GetInstance().GetParameterMC(oss.str()) ;
  } catch (const std::invalid_argument& e) {
    double W0VSmax = radiativecorrections::W0VS_max(Cs, a, mass_i, mass_f, Z, R, betaType, "BetaMinus");
    W = {W0VSmax, 0} ;
    DecayManager::GetInstance().RegisterParameterMC(oss.str(), W) ;
  }
  //std::cout << "w0vs : " << W[0] << "\n";
  std::random_device rd;
  std::mt19937 generator(rd()); 
  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  std::uniform_real_distribution<double> distribution_w0VSmax(0.0, W[0]); 

  int n_verif = n0VS ;
  while (n_verif == n0VS) {
    double w0VS_NR = distribution_w0VSmax(generator) ;

    double U[5] = {distribution(generator), distribution(generator), distribution(generator), distribution(generator), distribution(generator)};
    double E2 = 1. + (radiativecorrections::delta(mass_i, mass_f, "BetaMinus") - 1.) * U[0] ;
    double COS_NEUTRINO = 2. * U[1] - 1. ;
    
    if (w0VS_NR < radiativecorrections::W0VS(E2, COS_NEUTRINO, Cs, a, mass_i, mass_f, Z, R, betaType, "BetaMinus")) {
      n0VS += 1 ;
      double E10 = radiativecorrections::delta(mass_i, mass_f, "BetaMinus") - E2 ;
      double BETA = std::sqrt(1. - 1./std::pow(E2, 2) ) ;
      
      double COS_ELECTRON = 2. * U[2] - 1. ;
      double PHI_NEUTRINO = 2. * utilities::PI * U[3] ;
      double PHI_ELECTRON = 2. * utilities::PI * U[4] ;
      double SIN_NEUTRINO = std::sqrt( (1. - std::pow(COS_NEUTRINO, 2)) ) ;
      double SIN_ELECTRON = std::sqrt( (1. - std::pow(COS_ELECTRON, 2)) ) ;

      ublas::vector<double> n_ELECTRON (3); 
      n_ELECTRON[0]=SIN_ELECTRON * cos(PHI_ELECTRON); n_ELECTRON[1]=SIN_ELECTRON * sin(PHI_ELECTRON) ; n_ELECTRON[2]=COS_ELECTRON; 
      double n_ELECTRON_PRIME[3] = {-sin(PHI_ELECTRON), cos(PHI_ELECTRON), 0}; 
      double n_ELECTRON_SECOND[3] = {-COS_ELECTRON * cos(PHI_ELECTRON), -COS_ELECTRON * sin(PHI_ELECTRON), SIN_ELECTRON};

      double n_PERPENDICULAIRE_NEUTRINO[3] ; 
      ublas::vector<double> n_NEUTRINO (3);
      for (int j = 0; j < 3; j++) {
        n_PERPENDICULAIRE_NEUTRINO[j] = n_ELECTRON_PRIME[j] * cos(PHI_NEUTRINO) + n_ELECTRON_SECOND[j] * sin(PHI_NEUTRINO) ;
        n_NEUTRINO[j] = n_ELECTRON[j] * COS_NEUTRINO + n_PERPENDICULAIRE_NEUTRINO[j] * SIN_NEUTRINO ;
      }

      ublas::vector<double> velocity = -initState->GetVelocity();
      double eMomentum = std::sqrt(std::pow(E2*utilities::EMASSC2, 2) -std::pow(utilities::EMASSC2, 2.));
      double enubarMomentum = E10*utilities::EMASSC2 ;
      
      eFourMomentum(0) = E2*utilities::EMASSC2;
      eFourMomentum(1) = eMomentum*n_ELECTRON[0];
      eFourMomentum(2) = eMomentum*n_ELECTRON[1];
      eFourMomentum(3) = eMomentum*n_ELECTRON[2];

      /*enubarFourMomentum(0) = enubarMomentum ; 
      enubarFourMomentum(1) = enubarMomentum*n_NEUTRINO[0];
      enubarFourMomentum(2) = enubarMomentum*n_NEUTRINO[1];
      enubarFourMomentum(3) = enubarMomentum*n_NEUTRINO[2];
      enubar->SetMomentum(enubarFourMomentum); */
      
      e->SetMomentum(eFourMomentum);
      ThreeBodyDecay(velocity, e, enubar, recoil, n_NEUTRINO, Q);

      finalStates.push_back(e);
      finalStates.push_back(enubar);
      finalStates.push_back(recoil);
    } 
  }
  return finalStates;
}

std::vector<Particle*> BetaMinusPolarised::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates;

  // std::cout << "In BetaMinus Polarised Decay " << std::endl;
  // std::cout << "Address: " << initState << std::endl;
  std::ostringstream oss;
  oss << initState->GetCharge()+initState->GetNeutrons() << utilities::atoms[initState->GetCharge()];
  //std::cout << oss.str() << std::endl;
  Particle* recoil = DecayManager::GetInstance().GetNewParticle(oss.str(), initState->GetCharge()+1, initState->GetCharge()+initState->GetNeutrons());
  recoil->SetExcitationEnergy(daughterExEn);
    
  Particle* e = DecayManager::GetInstance().GetNewParticle("e-");
  Particle* enu = DecayManager::GetInstance().GetNewParticle("enubar");

  // std::cout << "Recoil " << recoil->GetCharge() << " " << recoil->GetNeutrons() << " " << recoil << std::endl;

  oss.str("");
  oss.clear();
  oss << "BetaMinusPolarised:Z" << recoil->GetCharge() << "A" << recoil->GetCharge() + recoil->GetNeutrons() << "Q" << Q;
  //Work in the COM frame
  ublas::vector<double> elFourMomentum (4);

  DecayManager& dm = DecayManager::GetInstance();
  double CS = dm.configOptions.couplingConstants.CS.real();
  double CSP = dm.configOptions.couplingConstants.CSP.real();
  double CV = dm.configOptions.couplingConstants.CV.real();
  double CVP = dm.configOptions.couplingConstants.CVP.real();
  double CA = dm.configOptions.couplingConstants.CA.real();
  double CAP = dm.configOptions.couplingConstants.CAP.real();
  double CT = dm.configOptions.couplingConstants.CT.real();
  double CTP = dm.configOptions.couplingConstants.CTP.real();
  double a_conf = dm.configOptions.couplingConstants.a;
  double b_conf = dm.configOptions.couplingConstants.b;

  double mf = 0.;
  double mgt = 0.;
  
  if (dm.configOptions.betaDecay.Default == "Fermi") {
    mf = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Gamow-Teller") {
    mgt = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Mixed") {
    mf = 1. ;
    mgt = 1. ;
  }
  else if (dm.configOptions.betaDecay.Default == "Auto") {
    std::string Type;
    try{
      Type = DecayManager::GetInstance().GetBetaType(oss.str());
    } catch (const std::invalid_argument& e) {
      Type = utilities::FindBetaType(initState, recoil);
      DecayManager::GetInstance().RegisterBetaType(oss.str(), Type);
    }
    if (Type == "Fermi") {
      mf = 1.;
    } else if (Type == "Gamow-Teller") {
      mgt = 1.;
    } else {  
      mgt = 1. ;
      mf = 1. ;
    }
  } 
  //mgt = 0. ;
  //mf = 1. ;
  double a = utilities::CalculateBetaNeutrinoAsymmetry(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf);
  double fierz = utilities::CalculateFierz(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf);

  //std::cout << "fierz " << fierz << " a " << a << std::endl;
  
  std::vector<std::vector<double> >* dist;
  try {
    dist = DecayManager::GetInstance().GetDistribution(oss.str());
  } catch (const std::invalid_argument& e) {
    bool advancedFermi = false;
    if (dm.configOptions.betaDecay.FermiFunction == "Advanced") { //// réécriture de la condition par SL 10/05/2023
      advancedFermi = true;
    }
    dist = spectrumGen->GenerateSpectrum(initState, recoil, Q); //// changement de Q en E0 par SL 10/05/2023
    double gamma = std::sqrt(1-std::pow(utilities::FINESTRUCTURE*recoil->GetCharge(), 2.));
    int i=0;
    for ( int i = 0; i<dist->size(); i++) {                                  ////// changement de boucle et element par ((*dist)[i]) par SL 10/05/2023
      double E = ((*dist)[i])[0]+utilities::EMASSC2;
      double SH = ((*dist)[i])[1];
      ((*dist)[i])[1] = SH*(1+gamma*fierz*utilities::EMASSC2/E);
    }
    DecayManager::GetInstance().RegisterDistribution(oss.str(), dist);
    }

  std::vector<double> p;

  double elEnergy = utilities::RandomFromDistribution(*dist)+utilities::EMASSC2;
  double elMomentum = std::sqrt(elEnergy*elEnergy-std::pow(utilities::EMASSC2, 2.));
  ublas::vector<double> enuDir = utilities::RandomDirection();

  p.push_back(1.);
  p.push_back(a*elMomentum/elEnergy);
  ublas::vector<double> eDir = utilities::GetParticleDirection(enuDir, p);

  elFourMomentum(0) = elEnergy;
  elFourMomentum(1) = elMomentum*eDir[0];
  elFourMomentum(2) = elMomentum*eDir[1];
  elFourMomentum(3) = elMomentum*eDir[2];

  e->SetMomentum(elFourMomentum);

  ublas::vector<double> velocity = -initState->GetVelocity();
  ThreeBodyDecay(velocity, e, enu, recoil, enuDir, Q);

  finalStates.push_back(recoil);
  finalStates.push_back(e);
  finalStates.push_back(enu);

  return finalStates;
}


////////////////////////////////////////////////////////////////

std::vector<Particle*> BetaPlusRadiative::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates;

  DecayManager& dm = DecayManager::GetInstance();

  double E0 = Q-2*utilities::EMASSC2;
  //std::cout << "E0 : " << E0 << "\n";
  

  //std::cout << "In BetaPlus Decay " << std::endl;
  //std::cout << "Address: " << initState << std::endl;
  std::ostringstream oss;
  oss << initState->GetCharge()+initState->GetNeutrons() << utilities::atoms[initState->GetCharge()-2];
  Particle* recoil = DecayManager::GetInstance().GetNewParticle(oss.str()+"RC", initState->GetCharge()-1, initState->GetCharge()+initState->GetNeutrons());
  recoil->SetExcitationEnergy(daughterExEn);
  Particle* pos = DecayManager::GetInstance().GetNewParticle("e+RC");
  Particle* enu = DecayManager::GetInstance().GetNewParticle("enuRC");
  Particle* gamma = DecayManager::GetInstance().GetNewParticle("gammaRC");

  oss.str("");
  oss.clear();

  oss << "BetaPlusRadiative:Z" << recoil->GetCharge() << "A" << recoil->GetCharge() + recoil->GetNeutrons() << "Q" << Q;
  //Work in the COM frame

  double mf = 0.;
  double mgt = 0.;

  ublas::vector<double> posFourMomentum (4);
  ublas::vector<double> enuFourMomentum (4);
  ublas::vector<double> gammaFourMomentum (4);

  /////////ajout de SL 12/05/2023//////////////   

  if (dm.configOptions.betaDecay.Default == "Fermi") {
    mf = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Gamow-Teller") {
    mgt = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Mixed") {
    mf = 1. ;
    mgt = 1. ;
  }
  else if (dm.configOptions.betaDecay.Default == "Auto") {
    std::string Type;
    try{
      Type = DecayManager::GetInstance().GetBetaType(oss.str());
    } catch (const std::invalid_argument& e) {
      Type = utilities::FindBetaType(initState, recoil);
      DecayManager::GetInstance().RegisterBetaType(oss.str(), Type);
    }
    if (Type == "Fermi") {
      mf = 1.;
    } else if (Type == "Gamow-Teller") {
      mgt = 1.;
    } else {
      mgt = 1. ;
      mf = 1. ;
    }
  } 
  //mf = 1. ;
  //mgt = 0. ;
  ////////////////////////////////////////////////

  double CS = dm.configOptions.couplingConstants.CS.real();
  double CSP = dm.configOptions.couplingConstants.CSP.real();
  double CV = dm.configOptions.couplingConstants.CV.real();
  double CVP = dm.configOptions.couplingConstants.CVP.real();
  double CA = dm.configOptions.couplingConstants.CA.real();
  double CAP = dm.configOptions.couplingConstants.CAP.real();
  double CT = dm.configOptions.couplingConstants.CT.real();
  double CTP = dm.configOptions.couplingConstants.CTP.real();
  double a_conf = dm.configOptions.couplingConstants.a;
  double b_conf = dm.configOptions.couplingConstants.b;

  double a = utilities::CalculateBetaNeutrinoAsymmetry(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf);
  double fierz = utilities::CalculateFierz(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf);

  //std::cout <<" b = " << fierz <<"\t a = " << a << std::endl;

  int nH = 0;
  double mass_i = initState->GetMass();
  double mass_f = mass_i - Q ; //recoil->GetMass() ; 
  /*std::cout << "mass f : " << mass_f << "\n";
  double m1 = mass_f + 2*utilities::EMASSC2 ;
  std::cout << "m1 : " << m1 << "\n";
  double m2 = mass_f + 2*(radiativecorrections::EMASSC2/radiativecorrections::UMASSC2)*std::pow(10, 6) ;
  std::cout << "m2 : " << m2 << "\n"; */
  //std::cout << "delta2 : " << mass_i - mass_f << "\n"; 
  double Cs = dm.configOptions.betaDecay.OmegaValue ;
  int A = recoil-> GetCharge() + recoil->GetNeutrons() ; 
  double R = utilities::ApproximateRadius(A) ;
  int Z = (recoil->GetCharge() - initState->GetCharge())*recoil->GetCharge();
  int betaType = (int)((Z > 0) - (Z < 0));
  Z = std::abs(Z) ;

  std::vector<double> W;
  try {
    W = DecayManager::GetInstance().GetParameterMC(oss.str()) ;
  } catch (const std::invalid_argument& e) {
    double WHmax = radiativecorrections::WH_max(1000000, Cs, a, mass_i, mass_f, Z, R, betaType, "BetaPlus");
    W = {WHmax, 0} ;
    DecayManager::GetInstance().RegisterParameterMC(oss.str(), W) ;
  }

  std::random_device rd;
  std::mt19937 generator(rd()); 
  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  std::uniform_real_distribution<double> distribution_wHmax(0.0, W[0]); 

  int n_verif = nH ;
  while (n_verif == nH) {
    double wH_NR = distribution_wHmax(generator) ;

    double U[8] = {distribution(generator), distribution(generator), distribution(generator), distribution(generator), distribution(generator), distribution(generator), distribution(generator), distribution(generator)};
    double E2 = 1. + (radiativecorrections::delta(mass_i, mass_f, "BetaPlus") - 1.) * U[0] ;
    double E10 = radiativecorrections::delta(mass_i, mass_f, "BetaPlus") - E2 ;
    double omega = Cs * E10 ;
    double K = omega * exp(-U[1] * log(Cs)) ;
    double E1 = E10 - K ;

    double BETA = std::sqrt(1. - 1./std::pow(E2, 2) ) ;
    double N = 0.5 * log( (1. + BETA)/(1. - BETA) ) ; 

    double COS_GAMMA = ( 1. - (1. + BETA) * exp(-2. * N * U[2]) )/BETA ;
    double COS_NEUTRINO = 2. * U[3] - 1. ;
    double COS_ELECTRON = 2. * U[4] - 1. ;

    double PHI_GAMMA = 2. * utilities::PI * U[5] ;
    double PHI_NEUTRINO = 2. * utilities::PI * U[6] ;
    double PHI_ELECTRON = 2. * utilities::PI * U[7] ; 

    double SIN_GAMMA = std::sqrt( (1. - std::pow(COS_GAMMA, 2)) ) ;
    double SIN_NEUTRINO = std::sqrt( (1. - std::pow(COS_NEUTRINO, 2)) ) ;
    double SIN_ELECTRON = std::sqrt( (1. - std::pow(COS_ELECTRON, 2)) ) ;
                
    ublas::vector<double> n_ELECTRON (3); 
    n_ELECTRON[0]=SIN_ELECTRON * cos(PHI_ELECTRON); n_ELECTRON[1]=SIN_ELECTRON * sin(PHI_ELECTRON) ; n_ELECTRON[2]=COS_ELECTRON; 
    double n_ELECTRON_PRIME[3] = {-sin(PHI_ELECTRON), cos(PHI_ELECTRON), 0}; 
    double n_ELECTRON_SECOND[3] = {-COS_ELECTRON * cos(PHI_ELECTRON), -COS_ELECTRON * sin(PHI_ELECTRON), SIN_ELECTRON};

    double n_PERPENDICULAIRE_GAMMA[3] ; 
    ublas::vector<double> n_GAMMA (3) ;
    ublas::vector<double> n_NEUTRINO (3);
    n_NEUTRINO[0]=SIN_NEUTRINO*cos(PHI_NEUTRINO) ; n_NEUTRINO[1]=SIN_NEUTRINO * sin(PHI_NEUTRINO) ; n_NEUTRINO[2]=COS_NEUTRINO;
    for (int j = 0; j < 3; j++) {
      n_PERPENDICULAIRE_GAMMA[j] = n_ELECTRON_PRIME[j] * cos(PHI_GAMMA) + n_ELECTRON_SECOND[j] * sin(PHI_GAMMA) ;
      n_GAMMA[j] = n_ELECTRON[j] * COS_GAMMA + n_PERPENDICULAIRE_GAMMA[j] * SIN_GAMMA ;
    }

    double N1_N2 = n_NEUTRINO[0]*n_ELECTRON[0] + n_NEUTRINO[1]*n_ELECTRON[1] + n_NEUTRINO[2]*n_ELECTRON[2] ;
    double N1_K = n_NEUTRINO[0]*n_GAMMA[0] + n_NEUTRINO[1]*n_GAMMA[1] + n_NEUTRINO[2]*n_GAMMA[2] ;
    if (wH_NR <= radiativecorrections::WH(E2, K, COS_GAMMA, N1_K, N1_N2, a, mass_i, mass_f, Z, R, betaType, "BetaPlus")) {
      nH += 1 ;

      /*double p2_k = E2 * K - BETA * E2 * K * COS_GAMMA ;
      double E1 = radiativecorrections::delta(mass_i, mass_f, "BetaPlus") - K - E2 ;
      double P_2 = 1/std::pow(K, 2) + 1/std::pow(p2_k, 2) - (2 * E2)/(K * p2_k) ;
      double p1_p2 = BETA * E1 * E2 * N1_N2 ;
      double p1_k = E1 * K *  N1_K ; 

      double h0 = E1 * ( -(E2 + K)*P_2 + K/(p2_k) ) ;
      double h1 = p1_p2 * ( -P_2 + 1/p2_k ) + p1_k * ( (E2 + K)/(K * p2_k) - 1/std::pow(p2_k, 2) ) ;

      double EPSa = mf - std::pow(radiativecorrections::LAMBDA, 2)*mgt/3 ;
      double EPS = mf + std::pow(radiativecorrections::LAMBDA, 2)*mgt ;
      double mbr = 16 * std::pow(radiativecorrections::FERMICONSTANT, 2) * std::pow( (mass_i/radiativecorrections::EMASSC2) , 2 ) * std::pow(radiativecorrections::e, 2) * (EPS * h0 + EPSa * h1) ;

      double g = (BETA * E2)/(2 * N * p2_k) ; 

      double whhh = ( K * BETA * E1 * E2 * mbr )/( std::pow(2, 13) * std::pow(utilities::PI, 8) * g * std::pow( (mass_i/radiativecorrections::EMASSC2) , 2 )) ;
        

      std::cout << "WH = nan" << "\n"; 
      std::cout << "E1 : " << E1 << "\n";
      std::cout << "E2 : " << E2 << "\n";
      std::cout << "K : " << K << "\n";
      std::cout << "delta : " << radiativecorrections::delta(mass_i, mass_f, "BetaPlus") << "\n" ;

      std::cout << "p2_k : " << p2_k << "\n";
      std::cout << "P_2 : " << P_2 << "\n";
      std::cout << "p1_p2 : " << p1_p2 << "\n";
      std::cout << "p1_k : " << p1_k << "\n";

      std::cout << "h0 : " << h0 << "\n";
      std::cout << "h1 : " << h1 << "\n";
      std::cout << "mbr : " << mbr << "\n";
      std::cout << "whhh : " << whhh << "\n";*/

      ublas::vector<double> velocity = -initState->GetVelocity();
      double posMomentum = std::sqrt(std::pow(E2*utilities::EMASSC2, 2) -std::pow(utilities::EMASSC2, 2.));
      double enuMomentum = E1*utilities::EMASSC2 ;
      double gammaMomentum = K*utilities::EMASSC2 ;
      
      posFourMomentum(0) = E2*utilities::EMASSC2;
      posFourMomentum(1) = posMomentum*n_ELECTRON[0];
      posFourMomentum(2) = posMomentum*n_ELECTRON[1];
      posFourMomentum(3) = posMomentum*n_ELECTRON[2];

      enuFourMomentum(0) = enuMomentum ; 
      enuFourMomentum(1) = enuMomentum*n_NEUTRINO[0];
      enuFourMomentum(2) = enuMomentum*n_NEUTRINO[1];
      enuFourMomentum(3) = enuMomentum*n_NEUTRINO[2];

      gammaFourMomentum(0) = gammaMomentum ;
      gammaFourMomentum(1) = gammaMomentum*n_GAMMA[0];
      gammaFourMomentum(2) = gammaMomentum*n_GAMMA[1];
      gammaFourMomentum(3) = gammaMomentum*n_GAMMA[2];

      pos->SetMomentum(posFourMomentum);
      enu->SetMomentum(enuFourMomentum);
      gamma->SetMomentum(gammaFourMomentum);
      FourBodyDecay(velocity, pos, enu, recoil, gamma);

      finalStates.push_back(pos);
      finalStates.push_back(enu);
      finalStates.push_back(recoil);
      finalStates.push_back(gamma) ;
    } 
  }
  return finalStates;
}


std::vector<Particle*> BetaPlusVirtualSoft::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates;

  DecayManager& dm = DecayManager::GetInstance();

  double E0 = Q-2*utilities::EMASSC2;
  //std::cout << "E0 : " << E0 << "\n";
  //std::cout << "daughter ex : " << daughterExEn << "\n";

  //std::cout << "In BetaPlus Decay " << std::endl;
  //std::cout << "Address: " << initState << std::endl;
  std::ostringstream oss;
  oss << initState->GetCharge()+initState->GetNeutrons() << utilities::atoms[initState->GetCharge()-2];
  Particle* recoil = DecayManager::GetInstance().GetNewParticle(oss.str()+"RC", initState->GetCharge()-1, initState->GetCharge()+initState->GetNeutrons());
  recoil->SetExcitationEnergy(daughterExEn);
  Particle* pos = DecayManager::GetInstance().GetNewParticle("e+RC");
  Particle* enu = DecayManager::GetInstance().GetNewParticle("enuRC");

  oss.str("");
  oss.clear();

  oss << "BetaPlusVirtualSoft:Z" << recoil->GetCharge() << "A" << recoil->GetCharge() + recoil->GetNeutrons() << "Q" << Q;
  //Work in the COM frame

  double mf = 0.;
  double mgt = 0.;

  ublas::vector<double> posFourMomentum (4);
  ublas::vector<double> enuFourMomentum (4);

  /////////ajout de SL 12/05/2023//////////////   

  if (dm.configOptions.betaDecay.Default == "Fermi") {
    mf = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Gamow-Teller") {
    mgt = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Mixed") {
    mf = 1. ;
    mgt = 1. ;
  }
  else if (dm.configOptions.betaDecay.Default == "Auto") {
    std::string Type;
    try{
      Type = DecayManager::GetInstance().GetBetaType(oss.str());
    } catch (const std::invalid_argument& e) {
      Type = utilities::FindBetaType(initState, recoil);
      DecayManager::GetInstance().RegisterBetaType(oss.str(), Type);
    }
    if (Type == "Fermi") {
      mf = 1.;
    } else if (Type == "Gamow-Teller") {
      mgt = 1.;
    } else {
      mgt = 1. ;
      mf = 1. ;
    }
  } 
  //std::cout << "mf : " << mf << "\n";
  //std::cout << "mgt : " << mgt << "\n";
  
  //mf = 1.;
  //mgt = 0.;
  ////////////////////////////////////////////////

  double CS = dm.configOptions.couplingConstants.CS.real();
  double CSP = dm.configOptions.couplingConstants.CSP.real();
  double CV = dm.configOptions.couplingConstants.CV.real();
  double CVP = dm.configOptions.couplingConstants.CVP.real();
  double CA = dm.configOptions.couplingConstants.CA.real();
  double CAP = dm.configOptions.couplingConstants.CAP.real();
  double CT = dm.configOptions.couplingConstants.CT.real();
  double CTP = dm.configOptions.couplingConstants.CTP.real();
  double a_conf = dm.configOptions.couplingConstants.a;
  double b_conf = dm.configOptions.couplingConstants.b;

  double a = utilities::CalculateBetaNeutrinoAsymmetry(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf);
  double fierz = utilities::CalculateFierz(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf);

  //std::cout <<" b = " << fierz <<"\t a = " << a << std::endl;
  //std::cout << "a : " << a << "\n";
  
  int n0VS = 0;
  double mass_i = initState->GetMass();
  double mass_f = mass_i - Q ; //recoil->GetMass() ; 
  /*std::cout << "mass f : " << mass_f << "\n";
  double m1 = mass_f + 2*utilities::EMASSC2 ;
  std::cout << "m1 : " << m1 << "\n";
  double m2 = mass_f + 2*(radiativecorrections::EMASSC2/radiativecorrections::UMASSC2)*std::pow(10, 6) ;
  std::cout << "m2 : " << m2 << "\n"; 
  std::cout << "e transi : " << mass_i - mass_f << "\n"; */
  double Cs = dm.configOptions.betaDecay.OmegaValue ;
  int A = recoil-> GetCharge() + recoil->GetNeutrons() ; 
  double R = utilities::ApproximateRadius(A) ;
  int Z = (recoil->GetCharge() - initState->GetCharge())*recoil->GetCharge();
  int betaType = (int)((Z > 0) - (Z < 0));
  Z = std::abs(Z) ;
  
  std::vector<double> W;
  try {
    W = DecayManager::GetInstance().GetParameterMC(oss.str()) ;
  } catch (const std::invalid_argument& e) {
    double W0VSmax = radiativecorrections::W0VS_max(Cs, a, mass_i, mass_f, Z, R, betaType, "BetaPlus");
    W = {W0VSmax, 0} ;
    DecayManager::GetInstance().RegisterParameterMC(oss.str(), W) ;
  }

  std::random_device rd;
  std::mt19937 generator(rd()); 
  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  std::uniform_real_distribution<double> distribution_w0VSmax(0.0, W[0]); 

  int n_verif = n0VS ;
  while (n_verif == n0VS) {
    double w0VS_NR = distribution_w0VSmax(generator) ;

    double U[5] = {distribution(generator), distribution(generator), distribution(generator), distribution(generator), distribution(generator)};
    double E2 = 1. + (radiativecorrections::delta(mass_i, mass_f, "BetaPlus") - 1 ) * U[0] ;
    double COS_NEUTRINO = 2. * U[1] - 1. ;
    
    /*std::cout << "E0 : " << E0 << "\n";
    std::cout << "mi - mf : " << mass_i - mass_f << "\n";
    std::cout << "Q : " << Q << "\n";
    std::cout << "delta : " << radiativecorrections::delta(mass_i, mass_f, "BetaPlus") << "\n";*/
    

    if (w0VS_NR < radiativecorrections::W0VS(E2, COS_NEUTRINO, Cs, a, mass_i, mass_f, Z, R, betaType, "BetaPlus")) {
      n0VS += 1 ;
      double E10 = radiativecorrections::delta(mass_i, mass_f, "BetaPlus") - E2 ;
      double BETA = std::sqrt(1. - 1./std::pow(E2, 2) ) ;
      double COS_ELECTRON = 2. * U[2] - 1. ;

      double PHI_NEUTRINO = 2. * utilities::PI * U[3] ;
      double PHI_ELECTRON = 2. * utilities::PI * U[4] ;
      double SIN_NEUTRINO = std::sqrt( (1. - std::pow(COS_NEUTRINO, 2)) ) ;
      double SIN_ELECTRON = std::sqrt( (1. - std::pow(COS_ELECTRON, 2)) ) ;

      ublas::vector<double> n_ELECTRON (3); 
      n_ELECTRON[0]=SIN_ELECTRON * cos(PHI_ELECTRON); n_ELECTRON[1]=SIN_ELECTRON * sin(PHI_ELECTRON) ; n_ELECTRON[2]=COS_ELECTRON; 
      double n_ELECTRON_PRIME[3] = {-sin(PHI_ELECTRON), cos(PHI_ELECTRON), 0}; 
      double n_ELECTRON_SECOND[3] = {-COS_ELECTRON * cos(PHI_ELECTRON), -COS_ELECTRON * sin(PHI_ELECTRON), SIN_ELECTRON};

      double n_PERPENDICULAIRE_NEUTRINO[3] ; 
      ublas::vector<double> n_NEUTRINO (3);
      for (int j = 0; j < 3; j++) {
        n_PERPENDICULAIRE_NEUTRINO[j] = n_ELECTRON_PRIME[j] * cos(PHI_NEUTRINO) + n_ELECTRON_SECOND[j] * sin(PHI_NEUTRINO) ;
        n_NEUTRINO[j] = n_ELECTRON[j] * COS_NEUTRINO + n_PERPENDICULAIRE_NEUTRINO[j] * SIN_NEUTRINO ;
      }

      ublas::vector<double> velocity = -initState->GetVelocity();
      double posMomentum = std::sqrt(std::pow(E2*utilities::EMASSC2, 2) -std::pow(utilities::EMASSC2, 2.));
      double enuMomentum = E10*utilities::EMASSC2 ;
      //std::cout << "pos : " << posMomentum << "\n";
      //std::cout << "enu : " << enuMomentum << "\n";
      /*if (E2*utilities::EMASSC2 > 1850) {
        std::cout << "E2 : " <<  E2*utilities::EMASSC2 << "\n";
        std::cout << "E10 : " <<  E10*utilities::EMASSC2 << "\n";
      }*/
      posFourMomentum(0) = E2*utilities::EMASSC2;
      posFourMomentum(1) = posMomentum*n_ELECTRON[0];
      posFourMomentum(2) = posMomentum*n_ELECTRON[1];
      posFourMomentum(3) = posMomentum*n_ELECTRON[2];

      enuFourMomentum(0) = enuMomentum ; 
      enuFourMomentum(1) = enuMomentum*n_NEUTRINO[0];
      enuFourMomentum(2) = enuMomentum*n_NEUTRINO[1];
      enuFourMomentum(3) = enuMomentum*n_NEUTRINO[2];

      pos->SetMomentum(posFourMomentum);
      enu->SetMomentum(enuFourMomentum);
      ThreeBodyDecay(velocity, pos, enu, recoil, n_NEUTRINO, E0);

      finalStates.push_back(pos);
      finalStates.push_back(enu);
      finalStates.push_back(recoil);
    } 
  }
  return finalStates;
}


std::vector<Particle*> BetaPlusPolarised::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates;

  DecayManager& dm = DecayManager::GetInstance();

  double E0 = Q-2*utilities::EMASSC2;

  //std::cout << "In BetaPlusPolarised Decay " << std::endl;
  //std::cout << "Address: " << initState << std::endl;
  std::ostringstream oss;
  oss << initState->GetCharge()+initState->GetNeutrons() << utilities::atoms[initState->GetCharge()-2];
  Particle* recoil = DecayManager::GetInstance().GetNewParticle(oss.str(), initState->GetCharge()-1, initState->GetCharge()+initState->GetNeutrons());
  recoil->SetExcitationEnergy(daughterExEn);
  Particle* pos = DecayManager::GetInstance().GetNewParticle("e+");
  Particle* enubar = DecayManager::GetInstance().GetNewParticle("enu");

  oss.str("");
  oss.clear();

  oss << "BetaPlusPolarised:Z" << recoil->GetCharge() << "A" << recoil->GetCharge() + recoil->GetNeutrons() << "Q" << Q;
  //Work in the COM frame
  ublas::vector<double> enubarDir = utilities::RandomDirection();

  ublas::vector<double> posFourMomentum (4);

  double mf = 0.;
  double mgt = 0.;

  /////////ajout de SL 12/05/2023//////////////

  if (dm.configOptions.betaDecay.Default == "Fermi") {
    mf = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Gamow-Teller") {
    mgt = 1.;
  } 
  else if (dm.configOptions.betaDecay.Default == "Mixed") {
    mf = 1. ;
    mgt = 1. ;
  }
  else if (dm.configOptions.betaDecay.Default == "Auto") {
    std::string Type;
    try{
      Type = DecayManager::GetInstance().GetBetaType(oss.str());
    } catch (const std::invalid_argument& e) {
      Type = utilities::FindBetaType(initState, recoil);
      DecayManager::GetInstance().RegisterBetaType(oss.str(), Type);
    }
    if (Type == "Fermi") {
      mf = 1.;
    } else if (Type == "Gamow-Teller") {
      mgt = 1.;
    } else {
      mgt = 1. ;
      mf = 1. ;
    }
  }
  //std::cout << "mf : " << mf << "\n";
  //std::cout << "mgt : " << mgt << "\n";
  //mf = 1.;
  //mgt = 0.;
  ////////////////////////////////////////////////

  double CS = dm.configOptions.couplingConstants.CS.real();
  double CSP = dm.configOptions.couplingConstants.CSP.real();
  double CV = dm.configOptions.couplingConstants.CV.real();
  double CVP = dm.configOptions.couplingConstants.CVP.real();
  double CA = dm.configOptions.couplingConstants.CA.real();
  double CAP = dm.configOptions.couplingConstants.CAP.real();
  double CT = dm.configOptions.couplingConstants.CT.real();
  double CTP = dm.configOptions.couplingConstants.CTP.real();
  double a_conf = dm.configOptions.couplingConstants.a;
  double b_conf = dm.configOptions.couplingConstants.b;

  double a = utilities::CalculateBetaNeutrinoAsymmetry(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf);
  double fierz = utilities::CalculateFierz(CS, CSP, CT, CTP, CV, CVP, CA, CAP, mf, mgt, a_conf, b_conf);

  //std::cout <<" b = " << fierz <<"\t a = " << a << std::endl;
  //std::cout << "a : " << a << "\n";

  std::vector<std::vector<double> >* dist;
  try {
    dist = DecayManager::GetInstance().GetDistribution(oss.str());
  } catch (const std::invalid_argument& e) {
    bool advancedFermi = false;
    if (dm.configOptions.betaDecay.FermiFunction == "Advanced") { //// réécriture de la condition par SL 10/05/2023
      advancedFermi = true;
    }
    dist = spectrumGen->GenerateSpectrum(initState, recoil, E0); //// changement de Q en E0 par SL 10/05/2023
    double gamma = std::sqrt(1-std::pow(utilities::FINESTRUCTURE*recoil->GetCharge(), 2.));
    int i=0;
    for ( int i = 0; i<dist->size(); i++) {                                  ////// changement de boucle et element par ((*dist)[i]) par SL 10/05/2023
      double E = ((*dist)[i])[0]+utilities::EMASSC2;
      double SH = ((*dist)[i])[1];
      ((*dist)[i])[1] = SH*(1+gamma*fierz*utilities::EMASSC2/E);
    }
    DecayManager::GetInstance().RegisterDistribution(oss.str(), dist);
    }


  // ublas::vector<std::vector<double> >* dist;
  double posEnergy = utilities::RandomFromDistribution(*dist) + utilities::EMASSC2;
  double posMomentum = std::sqrt(posEnergy*posEnergy-std::pow(utilities::EMASSC2, 2.));

  std::vector<double> p;
  p.push_back(1.);
  p.push_back(a*posMomentum/posEnergy);
  ublas::vector<double> posDir = utilities::GetParticleDirection(enubarDir, p);
  posFourMomentum(0) = posEnergy;
  posFourMomentum(1) = posMomentum*posDir[0];
  posFourMomentum(2) = posMomentum*posDir[1];
  posFourMomentum(3) = posMomentum*posDir[2];

  pos->SetMomentum(posFourMomentum);

  ublas::vector<double> velocity = -initState->GetVelocity();
  ThreeBodyDecay(velocity, pos, enubar, recoil, enubarDir, E0);


  finalStates.push_back(recoil);
  finalStates.push_back(pos);
  finalStates.push_back(enubar);

  return finalStates;
}

/////////////////////////////////////////////////////

std::vector<Particle*> ShellEC::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates ;
  ublas::vector<double> velocity = -initState->GetVelocity();
  
  std::ostringstream oss ;
  oss << initState->GetCharge()+initState->GetNeutrons() << utilities::atoms[initState->GetCharge()-2];
  Particle* recoil = DecayManager::GetInstance().GetNewParticle(oss.str()+"EC", initState->GetCharge()-1, initState->GetCharge()+initState->GetNeutrons());
  Particle* enu = DecayManager::GetInstance().GetNewParticle("enuEC");  
  recoil->SetExcitationEnergy(daughterExEn);
  
  //double mass_f = recoil -> GetMass() ;
  //std::cout << "mass f ec : " << mass_f << "\n";
  //double Excitation = decayChannel::GetParentExcitationEnergy() ;
  //std::cout << "excitation : " << excitation << "\n";
  
  TwoBodyDecay(velocity, enu, recoil, Q);
  
  finalStates.push_back(recoil) ;
  finalStates.push_back(enu) ;
  
  return finalStates;
}

/////////////////////////////////////////////////////

std::vector<Particle*> ConversionElectron::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates;

  ublas::vector<double> velocity = -initState->GetVelocity();
  Particle* recoil = DecayManager::GetInstance().GetNewParticle(initState->GetRawName());
  Particle* e = DecayManager::GetInstance().GetNewParticle("e+");
  recoil->SetExcitationEnergy(daughterExEn);

  TwoBodyDecay(velocity, recoil, e, Q);

  finalStates.push_back(recoil);
  finalStates.push_back(e);

  return finalStates;
}

std::vector<Particle*> Proton::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates;

  std::ostringstream oss;
  oss << initState->GetCharge()+initState->GetNeutrons() - 1 << utilities::atoms[initState->GetCharge()-2];
  Particle* recoil = DecayManager::GetInstance().GetNewParticle(oss.str(), initState->GetCharge()-1, initState->GetCharge()+initState->GetNeutrons()-1);
  Particle* p = DecayManager::GetInstance().GetNewParticle("p");
  recoil->SetExcitationEnergy(daughterExEn);

  ublas::vector<double> velocity = -initState->GetVelocity();
  TwoBodyDecay(velocity, recoil, p, Q);

  finalStates.push_back(recoil);
  finalStates.push_back(p);

  return finalStates;
}

std::vector<Particle*> Alpha::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates;

  std::ostringstream oss;
  oss << initState->GetCharge()+initState->GetNeutrons() - 4 << utilities::atoms[initState->GetCharge()-3];
  Particle* recoil = DecayManager::GetInstance().GetNewParticle(oss.str(), initState->GetCharge()-2, initState->GetCharge()+initState->GetNeutrons()-4);
  Particle* alpha = DecayManager::GetInstance().GetNewParticle("alpha");
  recoil->SetExcitationEnergy(daughterExEn);

  ublas::vector<double> velocity = -initState->GetVelocity();
  TwoBodyDecay(velocity, recoil, alpha, Q);

  finalStates.push_back(recoil);
  finalStates.push_back(alpha);

  return finalStates;
}

std::vector<Particle*> Gamma::Decay(Particle* initState, double Q, double daughterExEn) {
  std::vector<Particle*> finalStates;

  ublas::vector<double> velocity = -initState->GetVelocity();
  //std::cout << velocity[0] << "\n" ;
  //std::cout << velocity[1] << "\n" ;
  //std::cout << velocity[2] << "\n" ;
  Particle* recoil = DecayManager::GetInstance().GetNewParticle(initState->GetRawName());
  Particle* gamma = DecayManager::GetInstance().GetNewParticle("gamma");
  recoil->SetExcitationEnergy(daughterExEn);
  
  /*ublas::vector<double> momentum1 (4);
  ublas::vector<double> momentum2 (4);

  ublas::vector<double> dir = utilities::RandomDirection();

  double mass1 = recoil->GetMass();
  double mass2 = gamma->GetMass();
  
  double M = Q + mass1 + mass2;

  double p = 1./(2.*M)*std::sqrt((M*M-std::pow(mass1-mass2, 2.))*(M*M-std::pow(mass1+mass2, 2.)));
  
  double p1 = 1./(2.*M) * std::sqrt( std::pow(M*M + mass1*mass1 - mass2*mass2, 2) - 4*M*M*mass1*mass1 );
  double p2 = 1./(2.*M) * std::sqrt( std::pow(M*M + mass2*mass2 - mass1*mass1, 2) - 4*M*M*mass2*mass2 );

  double energy1 = std::sqrt(mass1*mass1+p1*p1);
  double energy2 = std::sqrt(mass2*mass2+p2*p2);
  
  momentum1(0) = energy1;
  momentum1(1) = p1*dir[0];
  momentum1(2) = p1*dir[1];
  momentum1(3) = p1*dir[2];
  
  momentum2(0) = energy2;
  momentum2(1) = -p2*dir[0];
  momentum2(2) = -p2*dir[1];
  momentum2(3) = -p2*dir[2];
  
  recoil->SetMomentum(momentum1) ;
  gamma->SetMomentum(utilities::LorentzBoost(velocity, momentum2)) ;*/
  //double mass_i = initState -> GetMass() ;
  //double mass_f = recoil -> GetMass() ;
  //std::cout << "mass i : " << mass_i << "\n";
  //std::cout << "mass f : " << mass_f << "\n" ;
  //std::cout << "diff mass : " << mass_i - mass_f - Q << "\n";
  //Q = daughterExEn ;
  //std::cout << "Q : " << Q << "\n";
  //std::cout << "daughter exc : " << daughterExEn << "\n";

  TwoBodyDecay(velocity, recoil, gamma, Q);

  finalStates.push_back(recoil);
  finalStates.push_back(gamma);

  return finalStates;
}

DecayMode::DecayMode() { }

DecayMode::~DecayMode() { }

void DecayMode::SetSpectrumGenerator(SpectrumGenerator* sg) {
  spectrumGen = sg;
}

BetaMinus::BetaMinus() { }

BetaPlus::BetaPlus() { }

BetaMinusRadiative::BetaMinusRadiative() { }

BetaPlusRadiative::BetaPlusRadiative() { }

BetaMinusVirtualSoft::BetaMinusVirtualSoft() { }

BetaPlusVirtualSoft::BetaPlusVirtualSoft() { }

BetaMinusPolarised::BetaMinusPolarised() { }

BetaPlusPolarised::BetaPlusPolarised() { }

ShellEC::ShellEC () { } 

ConversionElectron::ConversionElectron() { }

Proton::Proton () { }

Alpha::Alpha () { }

Gamma::Gamma () { }

}//End of CRADLE namespace
