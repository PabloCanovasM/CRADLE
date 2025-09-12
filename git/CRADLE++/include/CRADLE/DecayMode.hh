#ifndef DECAYMODE
#define DECAYMODE

#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include <string>

namespace CRADLE {

class Particle;
class SpectrumGenerator;

namespace ublas = boost::numeric::ublas;

class DecayMode {
  public:
    virtual std::vector<Particle*> Decay(Particle*, double, double) = 0;
    DecayMode();
    virtual ~DecayMode();

    void SetSpectrumGenerator(SpectrumGenerator*);

  protected:
    
    static void FourBodyDecay(ublas::vector<double>&, Particle*, Particle*, Particle*, Particle*) ;
    static void ThreeBodyDecay(ublas::vector<double>&, Particle*, Particle*, Particle*, ublas::vector<double>&, double);
    static void TwoBodyDecay(ublas::vector<double>&, Particle*, Particle*, double);
    SpectrumGenerator* spectrumGen;
};

class BetaMinus: public DecayMode {
  public:
    static BetaMinus& GetInstance() {
      static BetaMinus instance;
      return instance;
    }
    std::vector<Particle*> Decay(Particle*, double, double);

  protected:
    BetaMinus();
    BetaMinus(BetaMinus const& copy);
    BetaMinus& operator=(BetaMinus const& copy);
};

class BetaPlus: public DecayMode {
  public:
    static BetaPlus& GetInstance() {
      static BetaPlus instance;
      return instance;
    }
    std::vector<Particle*> Decay(Particle*, double, double);

  protected:
    BetaPlus();
    BetaPlus(BetaPlus const& copy);
    BetaPlus& operator=(BetaPlus const& copy);
};

///////////////////////////////////////////////////////////////
class BetaMinusRadiative: public DecayMode {
  public:
    static BetaMinusRadiative& GetInstance() {
      static BetaMinusRadiative instance;
      return instance;
    }
    std::vector<Particle*> Decay(Particle*, double, double);

  protected:
    BetaMinusRadiative();
    BetaMinusRadiative(BetaMinusRadiative const& copy);
    BetaMinusRadiative& operator=(BetaMinusRadiative const& copy);
};

class BetaPlusRadiative: public DecayMode {
  public:
    static BetaPlusRadiative& GetInstance() {
      static BetaPlusRadiative instance;
      return instance;
    }
    std::vector<Particle*> Decay(Particle*, double, double);

  protected:
    BetaPlusRadiative();
    BetaPlusRadiative(BetaPlusRadiative const& copy);
    BetaPlusRadiative& operator=(BetaPlusRadiative const& copy);
};

///////////////////////////////////////////////////////////////
class BetaMinusVirtualSoft: public DecayMode {
  public:
    static BetaMinusVirtualSoft& GetInstance() {
      static BetaMinusVirtualSoft instance;
      return instance;
    }
    std::vector<Particle*> Decay(Particle*, double, double);

  protected:
    BetaMinusVirtualSoft();
    BetaMinusVirtualSoft(BetaMinusVirtualSoft const& copy);
    BetaMinusVirtualSoft& operator=(BetaMinusVirtualSoft const& copy);
};

class BetaPlusVirtualSoft: public DecayMode {
  public:
    static BetaPlusVirtualSoft& GetInstance() {
      static BetaPlusVirtualSoft instance;
      return instance;
    }
    std::vector<Particle*> Decay(Particle*, double, double);

  protected:
    BetaPlusVirtualSoft();
    BetaPlusVirtualSoft(BetaPlusVirtualSoft const& copy);
    BetaPlusVirtualSoft& operator=(BetaPlusVirtualSoft const& copy);
};
//////////////////////////////////////////////////////////////////////

class ShellEC: public DecayMode {
  public:
    static ShellEC& GetInstance() {
      static ShellEC instance;
      return instance;
    }
    std::vector<Particle*> Decay(Particle*, double, double);

  protected:
    ShellEC();
    ShellEC(ShellEC const& copy);
    ShellEC& operator=(ShellEC const& copy);
};

//////////////////////////////////////////////////////////////////////

class ConversionElectron: public DecayMode {
    public:
    static ConversionElectron& GetInstance() {
      static ConversionElectron instance;
      return instance;
    }
    std::vector<Particle*> Decay(Particle*, double, double);

  protected:
    ConversionElectron();
    ConversionElectron(ConversionElectron const& copy);
    ConversionElectron& operator=(ConversionElectron const& copy);
};

class Proton: public DecayMode {
    public:
    static Proton& GetInstance() {
      static Proton instance;
      return instance;
    }
    std::vector<Particle*> Decay(Particle*, double, double);

  protected:
    Proton();
    Proton(Proton const& copy);
    Proton& operator=(Proton const& copy);
};

class Alpha: public DecayMode {
    public:
    static Alpha& GetInstance() {
      static Alpha instance;
      return instance;
    }
    std::vector<Particle*> Decay(Particle*, double, double);

  protected:
    Alpha();
    Alpha(Alpha const& copy);
    Alpha& operator=(Alpha const& copy);
};

class Gamma: public DecayMode {
    public:
    static Gamma& GetInstance() {
      static Gamma instance;
      return instance;
    }
    std::vector<Particle*> Decay(Particle*, double, double);

  protected:
    Gamma();
    Gamma(Gamma const& copy);
    Gamma& operator=(Gamma const& copy);
};

}//End of CRADLE namespace
#endif
