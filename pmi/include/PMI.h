// use C++ mode in Emacs: -*- C++ -*-

#ifndef PMI_HEADER
#define PMI_HEADER

static const char PMI_SccsId[] =
  "$Id: //tcad/gemini/D-2010.03/src/dessis/src/PMI.h#1 $, "
  "Copyright Synopsys, Inc. 2010-2010";

// definitions for the physical model interface in Sentaurus Device

#include "PMIFloat.h"
#include "PMIModels.h"



class PMI_Vertex_Base : public PMI_Vertex_Common_Base {

public:
  PMI_Vertex_Base (const PMI_Environment& env);

  virtual ~PMI_Vertex_Base ();

};



class PMI_MSC_Vertex_Base : public PMI_Vertex_Common_Base,
                            public PMI_MSC_Common_Base {

public:
  PMI_MSC_Vertex_Base (const PMI_Environment& env,
                       const std::string& msconfig_name,
                       int model_index,
                       const std::string& model_string);

  virtual ~PMI_MSC_Vertex_Base ();

};



class PMI_Vertex_Input_Base {

private:
  const PMI_Vertex_Common_Base* vertex_base;
  int vertex;

public:
  PMI_Vertex_Input_Base (const PMI_Vertex_Common_Base* vert_base, const int vert);
  virtual ~PMI_Vertex_Input_Base ();

  // coordinates of vertex in micrometers
  void ReadCoordinate (double& x, double& y, double& z) const;

  // distance in micrometers of vertex from closest insulator-semiconductor
  // interface if vertex is in semiconductor, or minus that distance otherwise
  double ReadDistanceFromSemiconductorInsulatorInterface () const;

  // distance in micrometers of vertex from closest high-k insulator
  double ReadDistanceFromHighkInsulator () const;

  // time in seconds during transient simulations
  pmi_float ReadTime () const;

  // step size in seconds during transient simulations
  pmi_float ReadTransientStepSize () const;

  // step type during transient simulations
  PMI_StepType ReadTransientStepType () const;

  // mole fractions
  pmi_float ReadxMoleFraction () const;
  pmi_float ReadyMoleFraction () const;

  // doping concentrations in cm^(-3)
  pmi_float ReadDoping (PMI_DopingSpecies species) const;

  // doping concentrations in cm^(-3)
  pmi_float ReadDoping (const char* SpeciesName) const;

  // was a user defined field loaded from a file?
  int IsUserFieldDefined (PMI_UserFieldIndex index) const;

  // read the value of a user defined field
  pmi_float ReadUserField (PMI_UserFieldIndex index) const;

  // write value to a user defined field
  void WriteUserField (PMI_UserFieldIndex index, pmi_float value) const;

  // stress value in Pascal
  pmi_float ReadStress (PMI_StressIndex index) const;

  // read MSC values for the actual vertex (if possible and available)
  // the function returns false if the MSC does not exist at
  // the actual vertex or if it is not available by internal reasons
  bool ReadMSCOccupations (const std::string& msc_name, pmi_float* values) const;

  pmi_float ReadeSHEDistribution (double energy) const;
  pmi_float ReadhSHEDistribution (double energy) const;
  pmi_float ReadeSHETotalDOS (double energy) const;
  pmi_float ReadhSHETotalDOS (double energy) const;
  pmi_float ReadeSHETotalGSV (double energy) const;
  pmi_float ReadhSHETotalGSV (double energy) const;


  // Synopsys internal use only
  void set_vertex (const int vert);

};



typedef std::vector<pmi_float>  sdevice_pmi_float_vector;
typedef std::vector<pmi_float*> sdevice_array_vector;



class sdevice_data {
  friend class PMIDeviceSupport;

public:
  typedef enum { vertex, edge, element, rivertex } sdevice_location;

  // element-edge coefficients: [element index][local element-edge index]
  const double*const* ReadCoefficient ();

  // element-vertex measures: [element index][local element-vertex index]
  const double*const* ReadMeasure ();

  // interface-vertex surface measures: [interface index][local interface-vertex index]
  const double*const* ReadSurfaceMeasure ();

  // return scalar data: [location index]
  const pmi_float* ReadScalar (sdevice_location location, std::string name);

  // return vector data: [mesh dimension index][location index]
  const pmi_float*const* ReadVector (sdevice_location location, std::string name);

  // set scalar data to values obtained from newvalue: [location index]
  void WriteScalar (sdevice_location location, std::string name, const pmi_float* newvalue);

  // return gradient of variable 'name'
  // NOTE: Actual implementation works for vertex-based datasets only
  const pmi_float*const* ReadGradient (sdevice_location location, std::string name);

  // return surface integral of the gradient of variable 'name' taken over the
  // boundary box divided by the box volume
  // NOTE: Actual implementation works for vertex-based datasets only
  const pmi_float* ReadFlux (sdevice_location location, std::string name);

  // get number of states of the given MSC
  size_t NumberOfMSCStates (const std::string& msc_name) const;

  // get name of the state of the MSC 'msc_name' with given state index
  bool ReadMSCStateName (const std::string& msc_name, size_t state_index,
                         std::string& state_name) const;

  // read occupations of given MSC and (bulk) region
  // for all region vertices and MSC states
  // (memory 'values' has to be allocated by the user: values[iv][is])
  // returns false if MSC is not defined or not computed in the region
  bool ReadMSCOccupations (const std::string& msc_name, const des_region* region,
                           pmi_float*const* values) const;

  pmi_float ReadeSHEDistribution (des_bulk* r, des_vertex* v, double energy) const;
  pmi_float ReadhSHEDistribution (des_bulk* r, des_vertex* v, double energy) const;
  pmi_float ReadeSHETotalDOS (des_bulk* r, double energy) const;
  pmi_float ReadhSHETotalDOS (des_bulk* r, double energy) const;
  pmi_float ReadeSHETotalGSV (des_bulk* r, double energy) const;
  pmi_float ReadhSHETotalGSV (des_bulk* r, double energy) const;


  // Synopsys internal use only
  class sdevice_data_identifier {
  public:
    sdevice_location w_location;
    std::string w_name;
    sdevice_data_identifier ();  // only used in std::map
    sdevice_data_identifier (sdevice_location location, std::string name);

    bool operator () (const sdevice_data_identifier id1, const sdevice_data_identifier id2) const;
  };

private:
  sdevice_data (const PMI_Device_Environment& env);
  ~sdevice_data ();

  typedef std::map <sdevice_data_identifier, pmi_float*, sdevice_data_identifier> sdevice_scalar_map;
  typedef std::map <sdevice_data_identifier, pmi_float**, sdevice_data_identifier> sdevice_vector_map;

  void delete_maps ();

  PMI_Device_Environment* w_env;

  double** w_coefficient;
  double** w_measure;
  double** w_surface_measure;

  sdevice_scalar_map w_scalar_map;
  sdevice_vector_map w_vector_map;

  sdevice_vector_map w_grad_map;
  sdevice_scalar_map w_flux_map;
};



class PMI_Device_Base : public PMI_Device_Common_Base {

public:
  PMI_Device_Base (const PMI_Device_Environment& env);

  virtual ~PMI_Device_Base ();

  // run time support

  // device data
  sdevice_data* Data () const;
};



class PMI_Device_Input_Base {

private:
  const PMI_Device_Base* device_base;

public:
  PMI_Device_Input_Base (const PMI_Device_Base* dev_base);
  virtual ~PMI_Device_Input_Base ();

  // time in seconds during transient simulations
  pmi_float ReadTime () const;

  // step size in seconds during transient simulations
  pmi_float ReadTransientStepSize () const;

  // step type during transient simulations
  PMI_StepType ReadTransientStepType () const;

};



class PMI_Avalanche_Base : public PMI_Vertex_Base {

private:
  const PMI_AvalancheDrivingForce drivingForce;

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_Avalanche_Base* avalanche_base, const int vertex);

    pmi_float F;                // driving force
    pmi_float t;                // lattice temperature
    pmi_float bg;               // bandgap
    pmi_float ct;               // carrier temperature
    pmi_float currentWoMob[3];  // current density (without mobility)
  };

  class Output {
  public:
    pmi_float alpha;  // ionization coefficient
  };

  PMI_Avalanche_Base (const PMI_Environment& env,
                      const PMI_AvalancheDrivingForce force);
  virtual ~PMI_Avalanche_Base ();

  PMI_AvalancheDrivingForce AvalancheDrivingForce () const { return drivingForce; }

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_Recombination_Base : public PMI_Vertex_Base {

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_Recombination_Base* recombination_base, const int vertex);

    pmi_float t;    // lattice temperature
    pmi_float n;    // electron density
    pmi_float p;    // hole density
    pmi_float nie;  // effective intrinsic density
    pmi_float f;    // absolute value of electric field
  };

  class Output {
  public:
    pmi_float r;  // recombination rate
  };

  PMI_Recombination_Base (const PMI_Environment& env);
  virtual ~PMI_Recombination_Base ();

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_DopingDepMobility_Base : public PMI_Vertex_Base {

private:
  const PMI_AnisotropyType anisoType;

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_DopingDepMobility_Base* dopingdepmobility_base, const int vertex);

    pmi_float n;         // electron density
    pmi_float p;         // hole density
    pmi_float t;         // lattice temperature
    pmi_float acceptor;  // total acceptor concentration (for random dopant fluctuations)
    pmi_float donor;     // total donor concentration (for random dopant fluctuations)
  };

  class Output {
  public:
    pmi_float m;  // doping dependent mobility
  };

  PMI_DopingDepMobility_Base (const PMI_Environment& env,
                              const PMI_AnisotropyType anisotype);
  virtual ~PMI_DopingDepMobility_Base ();

  PMI_AnisotropyType AnisotropyType () const { return anisoType; }

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_EnormalMobility_Base : public PMI_Vertex_Base {

private:
  const PMI_EnormalType enormalType;
  const PMI_AnisotropyType anisoType;

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_EnormalMobility_Base* enormalmobility_base, const int vertex);

    pmi_float dist;      // distance to nearest interface
    pmi_float pot;       // electrostatic potential
    pmi_float enorm;     // normal electric field
    pmi_float n;         // electron density
    pmi_float p;         // hole density
    pmi_float t;         // lattice temperature
    pmi_float ct;        // carrier temperature
    pmi_float acceptor;  // total acceptor concentration (for random dopant fluctuations)
    pmi_float donor;     // total donor concentration (for random dopant fluctuations)
  };

  class Output {
  public:
    pmi_float muinv;  // inverse of mobility degradation
  };

  PMI_EnormalMobility_Base (const PMI_Environment& env,
                            const PMI_EnormalType type,
                            const PMI_AnisotropyType anisotype);
  virtual ~PMI_EnormalMobility_Base ();

  PMI_EnormalType EnormalType () const { return enormalType; }
  PMI_AnisotropyType AnisotropyType () const { return anisoType; }

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_HighFieldMobility_Base : public PMI_Vertex_Base {

private:
  const PMI_HighFieldDrivingForce drivingForce;
  const PMI_AnisotropyType anisoType;

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_HighFieldMobility_Base* highfieldmobility_base, const int vertex);

    pmi_float pot;       // electrostatic potential
    pmi_float n;         // electron density
    pmi_float p;         // hole density
    pmi_float t;         // lattice temperature
    pmi_float ct;        // carrier temperature
    pmi_float mulow;     // low field mobility
    pmi_float F;         // driving force
    pmi_float acceptor;  // total acceptor concentration (for random dopant fluctuations)
    pmi_float donor;     // total donor concentration (for random dopant fluctuations)
  };

  class Output {
  public:
    pmi_float mu;  // mobility
  };

  PMI_HighFieldMobility_Base (const PMI_Environment& env,
                              const PMI_HighFieldDrivingForce force,
                              const PMI_AnisotropyType anisotype);
  virtual ~PMI_HighFieldMobility_Base ();

  PMI_HighFieldDrivingForce HighFieldDrivingForce () const { return drivingForce; }
  PMI_AnisotropyType AnisotropyType () const { return anisoType; }

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_HighFieldMobility2_Base : public PMI_Vertex_Base {

private:
  const int mindex;
  const std::string mstring;
  const PMI_AnisotropyType anisoType;

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_HighFieldMobility2_Base* highfieldmobility2_base, const int vertex);

    pmi_float mulow;     // low field mobility
    pmi_float n;         // electron density
    pmi_float p;         // hole density
    pmi_float T;         // lattice temperature
    pmi_float cT;        // carrier temperature
    pmi_float Epar;      // parallel electric field
    pmi_float gradQF;    // gradient of quasi-Fermi potential
    pmi_float EprodQF;   // product gradient QF and electric field
    pmi_float Na0;       // acceptor concentration
    pmi_float Nd0;       // donor concentration
  };

  class Output {
  public:
    pmi_float val;  // mobility
  };

  PMI_HighFieldMobility2_Base (const PMI_Environment& env,
                               const int model_index,
                               const std::string& model_string,
                               const PMI_AnisotropyType anisotype);
  virtual ~PMI_HighFieldMobility2_Base ();

  int model_index () const { return mindex; }
  const std::string& model_string () const { return mstring; }
  PMI_AnisotropyType AnisotropyType () const { return anisoType; }

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_BandGap_Base : public PMI_Vertex_Base {

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_BandGap_Base* bandgap_base, const int vertex);

    pmi_float t;  // lattice temperature
  };

  class Output {
  public:
    pmi_float bg;  // band bap
  };

  PMI_BandGap_Base (const PMI_Environment& env);
  virtual ~PMI_BandGap_Base ();

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_BandGapNarrowing_Base : public PMI_Vertex_Base {

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_BandGapNarrowing_Base* bandgapnarrowing_base, const int vertex);

    pmi_float acceptor;  // total acceptor concentration (for random dopant fluctuations)
    pmi_float donor;     // total donor concentration (for random dopant fluctuations)
  };

  class Output {
  public:
    pmi_float bgn;  // band bap narrowing
  };

  PMI_BandGapNarrowing_Base (const PMI_Environment& env);
  virtual ~PMI_BandGapNarrowing_Base ();

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_Affinity_Base : public PMI_Vertex_Base {

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_Affinity_Base* affinity_base, const int vertex);

    pmi_float t;  // lattice temperature
  };

  class Output {
  public:
    pmi_float affinity;  // electron affinity
  };

  PMI_Affinity_Base (const PMI_Environment& env);
  virtual ~PMI_Affinity_Base ();

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_EffectiveMass_Base : public PMI_Vertex_Base {

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_EffectiveMass_Base* effectivemass_base, const int vertex);

    pmi_float t;   // lattice temperature
    pmi_float bg;  // band gap
  };

  class Output {
  public:
    pmi_float m;  // effective mass
  };

  PMI_EffectiveMass_Base (const PMI_Environment& env);
  virtual ~PMI_EffectiveMass_Base ();

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_EnergyRelaxationTime_Base : public PMI_Vertex_Base {

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_EnergyRelaxationTime_Base* energyrelaxationtime_base, const int vertex);

    pmi_float ct;  // carrier temperature
  };

  class Output {
  public:
    pmi_float tau;  // energy relaxation time
  };

  PMI_EnergyRelaxationTime_Base (const PMI_Environment& env);
  virtual ~PMI_EnergyRelaxationTime_Base ();

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_Lifetime_Base : public PMI_Vertex_Base {

private:
  const PMI_LifetimeModel lifetimeModel;

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_Lifetime_Base* lifetime_base, const int vertex);

    pmi_float t;  // lattice temperature
  };

  class Output {
  public:
    pmi_float tau;  // lifetime
  };

  PMI_Lifetime_Base (const PMI_Environment& env,
                     const PMI_LifetimeModel model);
  virtual ~PMI_Lifetime_Base ();

  PMI_LifetimeModel LifetimeModel () const { return lifetimeModel; }

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_ThermalConductivity_Base : public PMI_Vertex_Base {

private:
  const PMI_AnisotropyType anisoType;

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_ThermalConductivity_Base* thermalconductivity_base, const int vertex);

    pmi_float t;  // lattice temperature
  };

  class Output {
  public:
    pmi_float kappa;  // thermal conductivity
  };

  PMI_ThermalConductivity_Base (const PMI_Environment& env,
                                const PMI_AnisotropyType anisotype);

  virtual ~PMI_ThermalConductivity_Base ();

  PMI_AnisotropyType AnisotropyType () const { return anisoType; }

  virtual void compute (const Input& input, Output& output) = 0;

};


class PMI_MetalResistivity_Base : public PMI_Vertex_Base {

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_MetalResistivity_Base* MetalResistivity_base, const int vertex);

    pmi_float t;  // lattice temperature
  };

  class Output {
  public:
    pmi_float Resist;  // metal resistivity
  };

  PMI_MetalResistivity_Base (const PMI_Environment& env);

  virtual ~PMI_MetalResistivity_Base ();

  virtual void compute (const Input& input, Output& output) = 0;

};




class PMI_HeatCapacity_Base : public PMI_Vertex_Base {

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_HeatCapacity_Base* heatcapacity_base, const int vertex);

    pmi_float t;  // lattice temperature
  };

  class Output {
  public:
    pmi_float c;  // heat capacity
  };

  PMI_HeatCapacity_Base (const PMI_Environment& env);

  virtual ~PMI_HeatCapacity_Base ();

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_MSC_HeatCapacity_Base : public PMI_MSC_Vertex_Base {

public:
  class Input : public PMI_Vertex_Input_Base {
  private:
    const NS_PMI_MSC::idata* w;
  public:
    Input (const PMI_MSC_HeatCapacity_Base* msc_heatcapacity_base,
           const NS_PMI_MSC::idata* idata);

    const pmi_float& n () const;   // electron density
    const pmi_float& p () const;   // hole density
    const pmi_float& T () const;   // lattice temperature
    const pmi_float& eT () const;  // electron temperature
    const pmi_float& hT () const;  // hole temperature
    const pmi_float& s (size_t ind) const;  // phase fraction
  };

  class Output {
  private:
    NS_PMI_MSC::odata* w;
  public:
    Output (NS_PMI_MSC::odata* odata);

    pmi_float& val ();  // heat capacity
  };

  PMI_MSC_HeatCapacity_Base (const PMI_Environment& env,
                             const std::string& msconfig_name,
                             const int model_index,
                             const std::string& model_string);

  virtual ~PMI_MSC_HeatCapacity_Base ();

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_MSC_ThermalConductivity_Base : public PMI_MSC_Vertex_Base {

private:
  const PMI_AnisotropyType anisoType;

public:
  class Input : public PMI_Vertex_Input_Base {
  private:
    const NS_PMI_MSC::idata* w;
  public:
    Input (const PMI_MSC_ThermalConductivity_Base* msc_thermalconductivity_base,
           const NS_PMI_MSC::idata* idata);

    const pmi_float& n () const;   // electron density
    const pmi_float& p () const;   // hole density
    const pmi_float& T () const;   // lattice temperature
    const pmi_float& eT () const;  // electron temperature
    const pmi_float& hT () const;  // hole temperature
    const pmi_float& s (size_t ind) const;  // phase fraction
  };

  class Output {
  private:
    NS_PMI_MSC::odata* w;
  public:
    Output (NS_PMI_MSC::odata* odata);

    pmi_float& val ();  // thermal conductivity
  };

  PMI_MSC_ThermalConductivity_Base (const PMI_Environment& env,
                                    const std::string& msconfig_name,
                                    const int model_index,
                                    const std::string& model_string,
                                    const PMI_AnisotropyType anisotype);

  virtual ~PMI_MSC_ThermalConductivity_Base ();

  PMI_AnisotropyType AnisotropyType () const { return anisoType; }

  virtual void compute (const Input& input, Output& output) = 0;

};


class PMI_MSC_Mobility_Base : public PMI_MSC_Vertex_Base {

private:
  const PMI_AnisotropyType anisoType;

public:
  class Input : public PMI_Vertex_Input_Base {
  private:
    const NS_PMI_MSC::idata* w;
  public:
    Input (const PMI_MSC_Mobility_Base* msc_mobility_base,
           const NS_PMI_MSC::idata* idata);

    const pmi_float& n () const;   // electron density
    const pmi_float& p () const;   // hole density
    const pmi_float& T () const;   // lattice temperature
    const pmi_float& eT () const;  // electron temperature
    const pmi_float& hT () const;  // hole temperature
    const pmi_float& s (size_t ind) const;  // phase fraction
  };

  class Output {
  private:
    NS_PMI_MSC::odata* w;
  public:
    Output (NS_PMI_MSC::odata* odata);

    pmi_float& val ();  // mobility
  };

  PMI_MSC_Mobility_Base (const PMI_Environment& env,
                         const std::string& msconfig_name,
                         const int model_index,
                         const std::string& model_string,
                         const PMI_AnisotropyType anisotype);

  virtual ~PMI_MSC_Mobility_Base ();

  PMI_AnisotropyType AnisotropyType () const { return anisoType; }

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_MSC_ApparentBandEdgeShift_Base : public PMI_MSC_Vertex_Base {

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_MSC_ApparentBandEdgeShift_Base* msc_apparentbandedgeshift_base,
           const int vertex);

    pmi_float n;   // electron density
    pmi_float p;   // hole density
    pmi_float T;   // lattice temperature
    pmi_float eT;  // electron temperature
    pmi_float hT;  // hole temperature
    std::vector<pmi_float> s;  // phase fraction
  };

  class Output {
  public:
    pmi_float val;  // apparent band-edge shift
  };

  PMI_MSC_ApparentBandEdgeShift_Base (const PMI_Environment& env,
                                      const std::string& msconfig_name,
                                      const int model_index);

  virtual ~PMI_MSC_ApparentBandEdgeShift_Base ();

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_TrapCaptureEmission_Base : public PMI_Vertex_Base {

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_TrapCaptureEmission_Base* trapcaptureemission_base, const int vertex);

    pmi_float n;   // electron density
    pmi_float p;   // hole density
    pmi_float t;   // lattice temperature
    pmi_float tn;  // electron temperature
    pmi_float tp;  // hole temperature
    pmi_float f;   // absolute value of electric field
    pmi_float nc;  // lattice effective state density for electrons
    pmi_float nv;  // lattice effective state density for holes
    pmi_float egeff; // effective band gap
  };

  class Output {
  public:
    pmi_float capture;   // capture rate
    pmi_float emission;  // emission rate
  };

  PMI_TrapCaptureEmission_Base (const PMI_Environment& env);
  virtual ~PMI_TrapCaptureEmission_Base ();

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_TrapEnergyShift_Base : public PMI_Vertex_Base {

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_TrapEnergyShift_Base* trapenergyshift_base, const int vertex);

    pmi_float f[3];  // electric field vector
    pmi_float t;     // lattice temperature
  };

  class Output {
  public:
    pmi_float shift;  // trap energy shift
  };

  PMI_TrapEnergyShift_Base (const PMI_Environment& env);
  virtual ~PMI_TrapEnergyShift_Base ();

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_ApparentBandEdgeShift_Base : public PMI_Vertex_Base {

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_ApparentBandEdgeShift_Base* apparentbandedgeshift_base, const int vertex);

    pmi_float n;  // electron density
    pmi_float p;  // hole density
    pmi_float t;  // lattice temperature
    pmi_float f;  // absolute value of electric field
  };

  class Output {
  public:
    pmi_float shift;  // apparent band edge shift
  };

  PMI_ApparentBandEdgeShift_Base (const PMI_Environment& env);
  virtual ~PMI_ApparentBandEdgeShift_Base ();

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_PiezoresistanceFactor_Base : public PMI_Vertex_Base {

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_PiezoresistanceFactor_Base* piezoresistancefactor_base, const int vertex);

    pmi_float Enormal;  // normal to interface electric field
  };

  class Output {
  public:
    pmi_float p11;  // piezoresistive prefactor
    pmi_float p12;  // piezoresistive prefactor
    pmi_float p44;  // piezoresistive prefactor
  };

  PMI_PiezoresistanceFactor_Base (const PMI_Environment& env);
  virtual ~PMI_PiezoresistanceFactor_Base ();

  virtual bool IsPrefactor () = 0;

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_DistributionFunction_Base : public PMI_Vertex_Base {

private:
  const PMI_SpeciesType speciesType;
  const char* speciesName;

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_DistributionFunction_Base* distributionfunction_base, const int vertex);

    pmi_float t;  // lattice temperature
  };

  class Output {
  public:
    pmi_float g;  // ionization factor
  };

  PMI_DistributionFunction_Base (const PMI_Environment& env,
                                 const char* name,
                                 const PMI_SpeciesType type = PMI_acceptor);
  virtual ~PMI_DistributionFunction_Base ();

  PMI_SpeciesType SpeciesType () const { return speciesType; }
  const char* SpeciesName () const { return speciesName; }

  // read parameter from Sentaurus Device parameter file
  // (override for PMI_Vertex_Base::ReadParameter)
  const PMIBaseParam* ReadParameter (const char* name) const;

  // initialize parameter from Sentaurus Device parameter file or from default value
  // (override for PMI_Vertex_Base::InitParameter)
  double InitParameter (const char* name, double defaultvalue) const;

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_SpatialDistributionFunction_Base : public PMI_Vertex_Base {

private:
  const char* HeavyIonType;

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_SpatialDistributionFunction_Base* spatialdistributionfunction_base, const int vertex);

    pmi_float w;  // radius (perpendicular distance from track)
    pmi_float l;  // coordinate along track
    pmi_float E;  // energy of heavy ion
  };

  class Output {
  public:
    pmi_float R;  // spatial distribution
  };

  PMI_SpatialDistributionFunction_Base (const PMI_Environment& env, const char* name);
  virtual ~PMI_SpatialDistributionFunction_Base ();

  const char* GetHeavyIonType () const { return HeavyIonType; }

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_StimEmissionCoeff_Base : public PMI_Vertex_Base {

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_StimEmissionCoeff_Base* stimemissioncoeff_base, const int vertex);

    pmi_float E;   // transition energy
    pmi_float n;   // electron density
    pmi_float p;   // hole density
    pmi_float et;  // electron temperature if hydro, lattice temperature otherwise
    pmi_float ht;  // hole temperature if hydro, lattice temperature otherwise
  };

  class Output {
  public:
    pmi_float rstim;  // stimulated emission coefficient
  };

  PMI_StimEmissionCoeff_Base (const PMI_Environment& env);
  virtual ~PMI_StimEmissionCoeff_Base ();

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_PhotonPhaseCoeff_Base : public PMI_Vertex_Base {

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_PhotonPhaseCoeff_Base* photonphasecoeff_base, const int vertex);

    pmi_float E;   // transition energy
    pmi_float n;   // electron density
    pmi_float p;   // hole density
    pmi_float et;  // electron temperature if hydro, lattice temperature otherwise
    pmi_float ht;  // hole temperature if hydro, lattice temperature otherwise
  };

  class Output {
  public:
    pmi_float pphase;  // photon phase coefficient
  };

  PMI_PhotonPhaseCoeff_Base (const PMI_Environment& env);
  virtual ~PMI_PhotonPhaseCoeff_Base ();

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_SponEmissionCoeff_Base : public PMI_Vertex_Base {

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_SponEmissionCoeff_Base* sponemissioncoeff_base, const int vertex);

    pmi_float E;   // transition energy
    pmi_float n;   // electron density
    pmi_float p;   // hole density
    pmi_float et;  // electron temperature if hydro, lattice temperature otherwise
    pmi_float ht;  // hole temperature if hydro, lattice temperature otherwise
  };

  class Output {
  public:
    pmi_float rspon;  // spontaneous emission coefficient
  };

  PMI_SponEmissionCoeff_Base (const PMI_Environment& env);
  virtual ~PMI_SponEmissionCoeff_Base ();

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_Stress_Base : public PMI_Vertex_Base {

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_Stress_Base* stress_base, const int vertex);
  };

  class Output {
  public:
    pmi_float stress_xx;  // xx component of stress
    pmi_float stress_yy;  // yy component of stress
    pmi_float stress_zz;  // zz component of stress
    pmi_float stress_yz;  // yz component of stress
    pmi_float stress_xz;  // xz component of stress
    pmi_float stress_xy;  // xy component of stress
  };

  PMI_Stress_Base (const PMI_Environment& env);
  virtual ~PMI_Stress_Base ();

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_SpaceFactor_Base : public PMI_Vertex_Base {

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_SpaceFactor_Base* spacefactor_base, const int vertex);
  };

  class Output {
  public:
    pmi_float spacefactor;  // trap space factor
  };

  PMI_SpaceFactor_Base (const PMI_Environment& env);
  virtual ~PMI_SpaceFactor_Base ();

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_Polarization_Base : public PMI_Vertex_Base {

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_Polarization_Base* polarization_base, const int vertex);
  };

  class Output {
  public:
    pmi_float pol [3];  // piezoelectric polarization
  };

  PMI_Polarization_Base (const PMI_Environment& env);
  virtual ~PMI_Polarization_Base ();

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_Absorption_Base : public PMI_Vertex_Base {

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_Absorption_Base* absorption_base, const int vertex);

    pmi_float energy;  // optical wave energy
    pmi_float t;       // device or lattice temperature
  };

  class Output {
  public:
    pmi_float alpha;  // optical absorption coefficient
  };

  PMI_Absorption_Base (const PMI_Environment& env);
  virtual ~PMI_Absorption_Base ();

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_RefractiveIndex_Base : public PMI_Vertex_Base {

public:
  class Input : public PMI_Vertex_Input_Base {
  public:
    Input (const PMI_RefractiveIndex_Base* refractiveindex_base, const int vertex);

    pmi_float energy;  // optical wave energy
    pmi_float t;       // device or lattice temperature
  };

  class Output {
  public:
    pmi_float Refract;  // optical refractive index
  };

  PMI_RefractiveIndex_Base (const PMI_Environment& env);
  virtual ~PMI_RefractiveIndex_Base ();

  virtual void compute (const Input& input, Output& output) = 0;

};



class PMI_CurrentPlot_Base : public PMI_Device_Base {

public:
  class Input : public PMI_Device_Input_Base {
  public:
    Input (const PMI_CurrentPlot_Base* currentplot_base);
  };

  class Output_Header {
  public:
    des_string_vector dataset;   // array of dataset names
    des_string_vector function;  // array of function names
  };

  class Output_Body {
  public:
    sdevice_pmi_float_vector value;  // array of plot values
  };

  PMI_CurrentPlot_Base (const PMI_Device_Environment& env);
  virtual ~PMI_CurrentPlot_Base ();

  virtual void compute_header (Output_Header& output) = 0;
  virtual void compute_body (const Input& input, Output_Body& output) = 0;

};



class PMI_HotCarrierInjection_Base : public PMI_Device_Base {

public:
  class Input : public PMI_Device_Input_Base {
  public:
    Input (const PMI_HotCarrierInjection_Base* hotcarrierinjection_base);

    des_regioninterface_vector regioninterfaces;  // region interfaces
                                                  // associated with model name
  };

  class Output {
  public:
    sdevice_array_vector gCurr;  // gate injection current in each vertex of
                                 // specified region interfaces
  };

  PMI_HotCarrierInjection_Base (const PMI_Device_Environment& env);
  virtual ~PMI_HotCarrierInjection_Base ();

  virtual void compute (const Input& input, Output& output) = 0;
     
};



class PMI_PostProcess_Base : public PMI_Device_Base {

public:
  class Input : public PMI_Device_Input_Base {
  public:
    Input (const PMI_PostProcess_Base* postprocess_base);
  };

  class Output {
  public:
  };

  PMI_PostProcess_Base (const PMI_Device_Environment& env);
  virtual ~PMI_PostProcess_Base ();

  virtual void compute (const Input& input, Output& output) = 0;
     
};



extern "C" {

typedef PMI_Avalanche_Base* new_PMI_Avalanche_Base_func
  (const PMI_Environment& env, const PMI_AvalancheDrivingForce force);

typedef PMI_Recombination_Base* new_PMI_Recombination_Base_func
  (const PMI_Environment& env);

typedef PMI_DopingDepMobility_Base* new_PMI_DopingDepMobility_Base_func
  (const PMI_Environment& env, const PMI_AnisotropyType anisotype);

typedef PMI_EnormalMobility_Base* new_PMI_EnormalMobility_Base_func
  (const PMI_Environment& env, const PMI_EnormalType type,
   const PMI_AnisotropyType anisotype);

typedef PMI_HighFieldMobility_Base* new_PMI_HighFieldMobility_Base_func
  (const PMI_Environment& env, const PMI_HighFieldDrivingForce force,
   const PMI_AnisotropyType anisotype);

typedef PMI_HighFieldMobility2_Base* new_PMI_HighFieldMobility2_Base_func
  (const PMI_Environment& env, const int model_index,
   const std::string& model_string, const PMI_AnisotropyType anisotype);

typedef PMI_BandGap_Base* new_PMI_BandGap_Base_func
  (const PMI_Environment& env);

typedef PMI_BandGapNarrowing_Base* new_PMI_BandGapNarrowing_Base_func
  (const PMI_Environment& env);

typedef PMI_Affinity_Base* new_PMI_Affinity_Base_func
  (const PMI_Environment& env);

typedef PMI_EffectiveMass_Base* new_PMI_EffectiveMass_Base_func
  (const PMI_Environment& env);

typedef PMI_EnergyRelaxationTime_Base* new_PMI_EnergyRelaxationTime_Base_func
  (const PMI_Environment& env);

typedef PMI_Lifetime_Base* new_PMI_Lifetime_Base_func
  (const PMI_Environment& env, const PMI_LifetimeModel model);

typedef PMI_ThermalConductivity_Base* new_PMI_ThermalConductivity_Base_func
  (const PMI_Environment& env, const PMI_AnisotropyType anisotype);

typedef PMI_MetalResistivity_Base* new_PMI_MetalResistivity_Base_func
  (const PMI_Environment& env);

typedef PMI_HeatCapacity_Base* new_PMI_HeatCapacity_Base_func
  (const PMI_Environment& env);

typedef PMI_MSC_HeatCapacity_Base* new_PMI_MSC_HeatCapacity_Base_func
  (const PMI_Environment& env, const std::string& msconfig_name,
   const int model_index, const std::string& model_string);

typedef PMI_MSC_ThermalConductivity_Base* new_PMI_MSC_ThermalConductivity_Base_func
  (const PMI_Environment& env, const std::string& msconfig_name,
   const int model_index, const std::string& model_string,
   const PMI_AnisotropyType anisotype);

typedef PMI_MSC_Mobility_Base* new_PMI_MSC_Mobility_Base_func
  (const PMI_Environment& env, const std::string& msconfig_name,
   const int model_index, const std::string& model_string,
   const PMI_AnisotropyType anisotype);

typedef PMI_MSC_ApparentBandEdgeShift_Base* new_PMI_MSC_ApparentBandEdgeShift_Base_func
  (const PMI_Environment& env, const std::string& msconfig_name,
   const int model_index);

typedef PMI_TrapCaptureEmission_Base* new_PMI_TrapCaptureEmission_Base_func
  (const PMI_Environment& env, int id);

typedef PMI_TrapEnergyShift_Base* new_PMI_TrapEnergyShift_Base_func
  (const PMI_Environment& env, int id);

typedef PMI_ApparentBandEdgeShift_Base* new_PMI_ApparentBandEdgeShift_Base_func
  (const PMI_Environment& env);

typedef PMI_PiezoresistanceFactor_Base* new_PMI_PiezoresistanceFactor_Base_func
  (const PMI_Environment& env);

typedef PMI_DistributionFunction_Base* new_PMI_DistributionFunction_Base_func
  (const PMI_Environment& env, const char* name, const PMI_SpeciesType type);

typedef PMI_SpatialDistributionFunction_Base* new_PMI_SpatialDistributionFunction_Base_func
  (const PMI_Environment& env, const char* HeavyIonName);

typedef PMI_StimEmissionCoeff_Base* new_PMI_StimEmissionCoeff_Base_func
  (const PMI_Environment& env);

typedef PMI_PhotonPhaseCoeff_Base* new_PMI_PhotonPhaseCoeff_Base_func
  (const PMI_Environment& env);

typedef PMI_SponEmissionCoeff_Base* new_PMI_SponEmissionCoeff_Base_func
  (const PMI_Environment& env);

typedef PMI_Stress_Base* new_PMI_Stress_Base_func
  (const PMI_Environment& env);

typedef PMI_SpaceFactor_Base* new_PMI_SpaceFactor_Base_func
  (const PMI_Environment& env);

typedef PMI_Polarization_Base* new_PMI_Polarization_Base_func
  (const PMI_Environment& env);

typedef PMI_Absorption_Base* new_PMI_Absorption_Base_func
  (const PMI_Environment& env);

typedef PMI_RefractiveIndex_Base* new_PMI_RefractiveIndex_Base_func
  (const PMI_Environment& env);

typedef PMI_CurrentPlot_Base* new_PMI_CurrentPlot_Base_func
  (const PMI_Device_Environment& env);

typedef PMI_HotCarrierInjection_Base* new_PMI_HotCarrierInjection_Base_func
  (const PMI_Device_Environment& env);

typedef PMI_PostProcess_Base* new_PMI_PostProcess_Base_func
  (const PMI_Device_Environment& env);


#ifndef PMI_DESSIS_SRC

// functions to be implemented by user

// virtual constructor
new_PMI_Avalanche_Base_func new_PMI_e_Avalanche_Base;
new_PMI_Avalanche_Base_func new_PMI_h_Avalanche_Base;

// virtual constructor
new_PMI_Recombination_Base_func new_PMI_Recombination_Base;

// virtual constructor
new_PMI_DopingDepMobility_Base_func new_PMI_DopingDep_e_Mobility_Base;
new_PMI_DopingDepMobility_Base_func new_PMI_DopingDep_h_Mobility_Base;

// virtual constructor
new_PMI_EnormalMobility_Base_func new_PMI_Enormal_e_Mobility_Base;
new_PMI_EnormalMobility_Base_func new_PMI_Enormal_h_Mobility_Base;

// virtual constructor
new_PMI_HighFieldMobility_Base_func new_PMI_HighField_e_Mobility_Base;
new_PMI_HighFieldMobility_Base_func new_PMI_HighField_h_Mobility_Base;

// virtual constructor
new_PMI_HighFieldMobility2_Base_func new_PMI_HighField_Mobility2_Base;

// virtual constructor
new_PMI_BandGap_Base_func new_PMI_BandGap_Base;

// virtual constructor
new_PMI_BandGapNarrowing_Base_func new_PMI_BandGapNarrowing_Base;

// virtual constructor
new_PMI_Affinity_Base_func new_PMI_Affinity_Base;

// virtual constructor
new_PMI_EffectiveMass_Base_func new_PMI_e_EffectiveMass_Base;
new_PMI_EffectiveMass_Base_func new_PMI_h_EffectiveMass_Base;

// virtual constructor
new_PMI_EnergyRelaxationTime_Base_func new_PMI_e_EnergyRelaxationTime_Base;
new_PMI_EnergyRelaxationTime_Base_func new_PMI_h_EnergyRelaxationTime_Base;

// virtual constructor
new_PMI_Lifetime_Base_func new_PMI_e_Lifetime_Base;
new_PMI_Lifetime_Base_func new_PMI_h_Lifetime_Base;

// virtual constructor
new_PMI_ThermalConductivity_Base_func new_PMI_ThermalConductivity_Base;

// virtual constructor
new_PMI_MetalResistivity_Base_func new_PMI_MetalResistivity_Base;

// virtual constructor
new_PMI_HeatCapacity_Base_func new_PMI_HeatCapacity_Base;

// virtual constructor
new_PMI_MSC_HeatCapacity_Base_func new_PMI_MSC_HeatCapacity_Base;

// virtual constructor
new_PMI_MSC_ThermalConductivity_Base_func new_PMI_MSC_ThermalConductivity_Base;

// virtual constructor
new_PMI_MSC_Mobility_Base_func new_PMI_MSC_Mobility_Base;

// virtual constructor
new_PMI_MSC_ApparentBandEdgeShift_Base_func new_PMI_MSC_ApparentBandEdgeShift_Base;

// virtual constructor
new_PMI_TrapCaptureEmission_Base_func new_PMI_TrapCaptureEmission_Base;

// virtual constructor
new_PMI_TrapEnergyShift_Base_func new_PMI_TrapEnergyShift_Base;

// virtual constructor
new_PMI_ApparentBandEdgeShift_Base_func new_PMI_ApparentBandEdgeShift_Base;

// virtual constructor
new_PMI_PiezoresistanceFactor_Base_func new_PMI_ePiezoresistanceFactor_Base;
new_PMI_PiezoresistanceFactor_Base_func new_PMI_hPiezoresistanceFactor_Base;

// virtual constructor
new_PMI_DistributionFunction_Base_func new_PMI_DistributionFunction_Base;

// virtual constructor
new_PMI_SpatialDistributionFunction_Base_func new_PMI_SpatialDistributionFunction_Base;

// virtual constructor
new_PMI_StimEmissionCoeff_Base_func new_PMI_StimEmissionCoeff_Base;

// virtual constructor
new_PMI_PhotonPhaseCoeff_Base_func new_PMI_PhotonPhaseCoeff_Base;

// virtual constructor
new_PMI_SponEmissionCoeff_Base_func new_PMI_SponEmissionCoeff_Base;

// virtual constructor
new_PMI_Stress_Base_func new_PMI_Stress_Base;

// virtual constructor
new_PMI_SpaceFactor_Base_func new_PMI_SpaceFactor_Base;

// virtual constructor
new_PMI_Polarization_Base_func new_PMI_Polarization_Base;

// virtual constructor
new_PMI_Absorption_Base_func new_PMI_Absorption_Base;

// virtual constructor
new_PMI_RefractiveIndex_Base_func new_PMI_RefractiveIndex_Base;

// virtual constructor
new_PMI_CurrentPlot_Base_func new_PMI_CurrentPlot_Base;

// virtual constructor
new_PMI_HotCarrierInjection_Base_func new_PMI_e_HotCarrierInjection_Base;
new_PMI_HotCarrierInjection_Base_func new_PMI_h_HotCarrierInjection_Base;

// virtual constructor
new_PMI_PostProcess_Base_func new_PMI_PostProcess_Base;

#endif

} // extern "C"

#endif
