// use C++ mode in Emacs: -*- C++ -*-

#ifndef PMIModels_HEADER
#define PMIModels_HEADER

#define PMI_VERSION 170


static const char PMIModels_SccsId[] =
  "$Id: //tcad/gemini/D-2010.03/src/dessis/src/PMIModels.h#1 $, "
  "Copyright Synopsys, Inc. 2000-2010";

// definitions for the physical model interface in Sentaurus Device

#include <map>
#include <vector>
#include <string>
#include <iostream>


#if defined(_MSC_VER)
#ifdef PMI_DESSIS_SRC
#define PMI_EXTERNAL __declspec (dllexport)
#else
#define PMI_EXTERNAL __declspec (dllimport)
#endif
#else
#define PMI_EXTERNAL
#endif

// backward compatibility only
#define PMI_Dessis_Interface PMI_Vertex_Interface
#define PMI_Dessis_Device_Interface PMI_Device_Interface
#define PMI_Boron PMI_BoronActive
#define PMI_Indium PMI_IndiumActive
#define PMI_PDopant PMI_PDopantActive
#define PMI_Phosphorus PMI_PhosphorusActive
#define PMI_Arsenic PMI_ArsenicActive
#define PMI_Antimony PMI_AntimonyActive
#define PMI_Nitrogen PMI_NitrogenActive
#define PMI_NDopant PMI_NDopantActive


class PMI_EXTERNAL PMIBaseParam {
public:

  PMIBaseParam ();
  virtual ~PMIBaseParam ();

  // name of parameter
  virtual const char* Name () const = 0;

  // conversion operators
  virtual operator double () const = 0;

  // assignment operators
  virtual PMIBaseParam& operator = (const double val) = 0;
  virtual PMIBaseParam& operator = (const PMIBaseParam& par) = 0;

  // print the parameter
  virtual void Print (std::ostream& stream) const;

};



PMI_EXTERNAL std::ostream& operator << (std:: ostream& stream, const PMIBaseParam* par);



class PMIDeviceSupport;
class PMIModel;
class PMI_EXTERNAL PMI_Environment;
class PMI_EXTERNAL PMI_Device_Environment;



PMI_EXTERNAL enum PMI_DopingSpecies {
  // Acceptors
  PMI_BoronActive,      // active Boron concentration
  PMI_BoronChemical,    // chemical Boron concentration
  PMI_IndiumActive,     // active Indium concentration
  PMI_IndiumChemical,   // chemical Indium concentration
  PMI_PDopantActive,    // active PDopant concentration
  PMI_PDopantChemical,  // chemical PDopant concentration
  PMI_Acceptor,         // total acceptor concentration

  // incomplete ionization entries
  PMI_AcceptorMinus,  // total incomplete ionization acceptor concentration
  PMI_AcceptorMinusPer_hDensity,
  PMI_AcceptorMinusPerT,

  // Donors
  PMI_PhosphorusActive,    // active Phosphorus concentration
  PMI_PhosphorusChemical,  // chemical Phosphorus concentration
  PMI_ArsenicActive,       // active Arsenic concentration
  PMI_ArsenicChemical,     // chemical Arsenic concentration
  PMI_AntimonyActive,      // active Antimony concentration
  PMI_AntimonyChemical,    // chemical Antimony concentration
  PMI_NitrogenActive,      // active Nitrogen concentration
  PMI_NitrogenChemical,    // chemical Nitrogen concentration
  PMI_NDopantActive,       // active NDopant concentration
  PMI_NDopantChemical,     // chemical NDopant concentration
  PMI_Donor,               // total donor concentration

  // incomplete ionization entries
  PMI_DonorPlus,  // total incomplete ionization donor concentration
  PMI_DonorPlusPer_eDensity,
  PMI_DonorPlusPerT
};



PMI_EXTERNAL enum PMI_UserFieldIndex {
  PMI_UserField0,
  PMI_UserField1,
  PMI_UserField2,
  PMI_UserField3,
  PMI_UserField4,
  PMI_UserField5,
  PMI_UserField6,
  PMI_UserField7,
  PMI_UserField8,
  PMI_UserField9,
  PMI_UserField10,
  PMI_UserField11,
  PMI_UserField12,
  PMI_UserField13,
  PMI_UserField14,
  PMI_UserField15,
  PMI_UserField16,
  PMI_UserField17,
  PMI_UserField18,
  PMI_UserField19,
  PMI_UserField20,
  PMI_UserField21,
  PMI_UserField22,
  PMI_UserField23,
  PMI_UserField24,
  PMI_UserField25,
  PMI_UserField26,
  PMI_UserField27,
  PMI_UserField28,
  PMI_UserField29,
  PMI_UserField30,
  PMI_UserField31,
  PMI_UserField32,
  PMI_UserField33,
  PMI_UserField34,
  PMI_UserField35,
  PMI_UserField36,
  PMI_UserField37,
  PMI_UserField38,
  PMI_UserField39,
  PMI_UserField40,
  PMI_UserField41,
  PMI_UserField42,
  PMI_UserField43,
  PMI_UserField44,
  PMI_UserField45,
  PMI_UserField46,
  PMI_UserField47,
  PMI_UserField48,
  PMI_UserField49,
  PMI_UserField50,
  PMI_UserField51,
  PMI_UserField52,
  PMI_UserField53,
  PMI_UserField54,
  PMI_UserField55,
  PMI_UserField56,
  PMI_UserField57,
  PMI_UserField58,
  PMI_UserField59,
  PMI_UserField60,
  PMI_UserField61,
  PMI_UserField62,
  PMI_UserField63,
  PMI_UserField64,
  PMI_UserField65,
  PMI_UserField66,
  PMI_UserField67,
  PMI_UserField68,
  PMI_UserField69,
  PMI_UserField70,
  PMI_UserField71,
  PMI_UserField72,
  PMI_UserField73,
  PMI_UserField74,
  PMI_UserField75,
  PMI_UserField76,
  PMI_UserField77,
  PMI_UserField78,
  PMI_UserField79,
  PMI_UserField80,
  PMI_UserField81,
  PMI_UserField82,
  PMI_UserField83,
  PMI_UserField84,
  PMI_UserField85,
  PMI_UserField86,
  PMI_UserField87,
  PMI_UserField88,
  PMI_UserField89,
  PMI_UserField90,
  PMI_UserField91,
  PMI_UserField92,
  PMI_UserField93,
  PMI_UserField94,
  PMI_UserField95,
  PMI_UserField96,
  PMI_UserField97,
  PMI_UserField98,
  PMI_UserField99
};



PMI_EXTERNAL enum PMI_StressIndex {
  PMI_StressXX,
  PMI_StressYY,
  PMI_StressZZ,
  PMI_StressYZ,
  PMI_StressXZ,
  PMI_StressXY
};



PMI_EXTERNAL enum PMI_StepType {
  PMI_UndefStepType = 0,
  PMI_TR = 1,
  PMI_BDF = 2,
  PMI_BE = 3
};



class PMI_EXTERNAL PMI_Vertex_Common_Base {
  friend class PMI_Vertex_Input_Base;

protected:
  PMI_Environment* pmi_env;
  PMIModel* pmi_model;

public:
  PMI_Vertex_Common_Base (const PMI_Environment& env);

  virtual ~PMI_Vertex_Common_Base ();

  // run time support

  // return name of this PMI model
  const char* Name () const;

  // region name
  const char* ReadRegionName () const;

  // region material
  const char* ReadRegionMaterial () const;

  // device name
  const char* ReadDeviceName () const;

  // read parameter from Sentaurus Device parameter file
  const PMIBaseParam* ReadParameter (const char* name) const;

  // initialize parameter from Sentaurus Device parameter file or from default value
  double InitParameter (const char* name, double defaultvalue) const;

  // dimension of mesh (1, 2, or 3)
  int ReadDimension () const;

  // read reference coordinate system:
  //
  //     [ 1  0  0 ]
  //     [ 0  1  0 ] is unified coordinate system
  //     [ 0  0  1 ]
  //     
  //     otherwise dfise coordinate system
  //
  void ReadReferenceCoordinates (double ref [3][3]) const;

  // number of worker threads in parallel assembly
  static int ReadNumberOfThreads ();

  // thread id of this thread
  // 0 <= thread id < number of threads
  static int ReadThreadId ();

  // read number of states of MSConfig
  size_t NumberOfMSConfigStates (const std::string& msconfig_name) const;

  // read state name of MSConfig state
  const std::string& MSConfigStateName (const std::string& msconfig_name,
                                        size_t state_index) const;

};



class PMI_EXTERNAL PMI_Vertex_Interface : public PMI_Vertex_Common_Base {

public:
  PMI_Vertex_Interface (const PMI_Environment& env);

  virtual ~PMI_Vertex_Interface ();

  // run time support

  // coordinates of vertex in micrometers (available in Compute functions)
  void ReadCoordinate (double& x, double& y, double& z) const;

  // distance in micrometers of vertex from closest insulator-semiconductor
  // interface if vertex is in semiconductor, or minus that distance
  // otherwise (available in Compute functions)
  double ReadDistanceFromSemiconductorInsulatorInterface() const;

  // distance in micrometers of vertex from closest high-k insulator
  // (not available in Compute functions yet)
  double ReadDistanceFromHighkInsulator() const;

  // time in seconds during transient simulations
  double ReadTime () const;

  // step size in seconds during transient simulations
  double ReadTransientStepSize () const;

  // step type during transient simulations
  PMI_StepType ReadTransientStepType () const;

  // mole fractions
  double ReadxMoleFraction () const;
  double ReadyMoleFraction () const;

  // doping concentrations in cm^(-3)
  double ReadDoping (PMI_DopingSpecies species) const;

  // doping concentrations in cm^(-3)
  double ReadDoping (const char* SpeciesName) const;

  // was a user defined field loaded from a file?
  int IsUserFieldDefined (PMI_UserFieldIndex index) const;

  // read the value of a user defined field
  double ReadUserField (PMI_UserFieldIndex index) const;

  // write value to a user defined field
  void WriteUserField (PMI_UserFieldIndex index, double value) const;

  // stress value in Pascal
  double ReadStress (PMI_StressIndex index) const;

  // read MSC values for the actual vertex (if possible and available)
  // the function returns false if the MSC does not exist at
  // the actual vertex or if it is not available by internal reasons
  bool ReadMSCOccupations (
    const std::string& msc_name, double* values ) const;

  double ReadeSHEDistribution (double energy) const;
  double ReadhSHEDistribution (double energy) const;
  double ReadeSHETotalDOS (double energy) const;
  double ReadhSHETotalDOS (double energy) const;
  double ReadeSHETotalGSV (double energy) const;
  double ReadhSHETotalGSV (double energy) const;

  // the following functions are not documented, use at your own risk!

  // print a list of recognized Sentaurus Device field names
  void PrintFieldNames () const;

  // convert Sentaurus Device field name into index (returns -1 on failure)
  int FieldName2Index (const char* fieldname) const;

  // read the value of an Sentaurus Device field
  double ReadField (int index) const;

};

//----------------------------------------------------------------------
namespace NS_PMI_MSC {
    class Vertex_Interface;
    class idata;
    class odata;
}



class PMI_EXTERNAL PMI_MSC_Common_Base {

private:
  NS_PMI_MSC::Vertex_Interface* w_impl;

public:
  // not to be used by user models
  NS_PMI_MSC::Vertex_Interface* impl ();

  PMI_MSC_Common_Base (const PMI_Environment& env,
                       const std::string& msconfig_name,
                       int model_index,
                       const std::string& model_string);

  virtual ~PMI_MSC_Common_Base ();

  // the name of the reference MSC
  const std::string& msconfig_name () const;

  // the number of states of the reference MSC
  size_t nb_states () const;

  // the name of the state (of the reference MSC) for the given index
  const std::string& state (size_t index) const;

  // the model index for free use by the PMI model
  int model_index () const;

  // the model string for free use by the PMI model
  const std::string& model_string () const;

  // initialize the internal parameters
  // (called by the simulator at least if the parameter values
  // have changed e.g. due to quasitationary ramping)
  //
  // the user should overload this function
  // and update all relevant parameters by a call to InitParameter
  virtual void init_parameter () {}

};



// the base class for the multistate configuration (MSC) models
class PMI_EXTERNAL PMI_MSC_Vertex_Interface : public PMI_Vertex_Interface,
                                              public PMI_MSC_Common_Base {

public:
  PMI_MSC_Vertex_Interface (const PMI_Environment& env,
                            const std::string& msconfig_name,
                            int model_index,
                            const std::string& model_string);

  virtual ~PMI_MSC_Vertex_Interface ();

};

//----------------------------------------------------------------------



class des_vertex;
class des_edge;
class des_element;
class des_region;
class des_bulk;
class des_contact;
class des_regioninterface;

typedef std::vector<std::string>          des_string_vector;
typedef std::vector<double>               des_double_vector;
typedef std::vector<des_vertex*>          des_vertex_vector;
typedef std::vector<des_edge*>            des_edge_vector;
typedef std::vector<des_element*>         des_element_vector;
typedef std::vector<des_region*>          des_region_vector;
typedef std::vector<des_regioninterface*> des_regioninterface_vector;
typedef std::vector<double*>              des_array_vector;



class PMI_EXTERNAL des_vertex {
  friend class PMIDeviceSupport;

public:
  // index for vertex data
  size_t index () const { return w_index; }

  // vertex coordinates
  const double* coord () const { return w_coord; }

  // has this vertex the same coordinates as vertex v?
  bool equal_coord (des_vertex* v) const;

  // number of edges connected to vertex
  size_t size_edge () const { return w_edge.size (); }

  // return edge i connected to vertex
  des_edge* edge (size_t i) const { return w_edge.at (i); }

  // number of elements connected to vertex
  size_t size_element () const { return w_element.size (); }

  // return element i connected to vertex
  des_element* element (size_t i) const { return w_element.at (i); }

  // number of regions containing vertex
  size_t size_region () const { return w_region.size (); }

  // return region i containing vertex
  des_region* region (size_t i) const { return w_region.at (i); }

  // number of region interfaces containing vertex
  size_t size_regioninterface () const { return w_regioninterface.size (); }

  // return regioninterface i containing vertex
  des_regioninterface* regioninterface (size_t i) const { return w_regioninterface.at (i); }

private:
  des_vertex (size_t index, double* coord);

  void add_edge (des_edge* e) { w_edge.push_back (e); }
  void add_element (des_element* e) { w_element.push_back (e); }
  void add_region (des_region* r) { w_region.push_back (r); }
  void add_regioninterface (des_regioninterface* ri) { w_regioninterface.push_back (ri); }

  size_t w_index;
  double w_coord [3];
  des_edge_vector w_edge;
  des_element_vector w_element;
  des_region_vector w_region;
  des_regioninterface_vector w_regioninterface;
};



class PMI_EXTERNAL des_edge {
  friend class PMIDeviceSupport;

public:
  // index for edge data
  size_t index () const { return w_index; }

  // return first vertex of edge
  des_vertex* start () const { return w_start; }

  // return second vertex of edge
  des_vertex* end () const { return w_end; }

  // number of elements connected to edge
  size_t size_element () const { return w_element.size (); }

  // return element i connected to edge
  des_element* element (size_t i) const { return w_element.at (i); }

  // number of regions containing edge
  size_t size_region () const { return w_region.size (); }

  // return region i containing edge
  des_region* region (size_t i) const { return w_region.at (i); }

private:
  des_edge (size_t index, des_vertex* start, des_vertex* end);

  void add_element (des_element* e) { w_element.push_back (e); }
  void add_region (des_region* r) { w_region.push_back (r); }

  size_t w_index;
  des_vertex* w_start;
  des_vertex* w_end;
  des_element_vector w_element;
  des_region_vector w_region;
};



class PMI_EXTERNAL des_element {
  friend class PMIDeviceSupport;

public:
  typedef enum { point, line, triangle, rectangle, tetrahedron,
                 pyramid, prism, cuboid, tetrabrick } des_type;

  // index for element data
  size_t index () const { return w_index; }

  // type of element
  des_type type () const { return w_type; }

  // number of vertices in element
  size_t size_vertex () const { return w_vertex.size (); }

  // return vertex i in element
  des_vertex* vertex (size_t i) const { return w_vertex.at (i); }

  // number of edges in element
  size_t size_edge () const { return w_edge.size (); }

  // return edge i in element
  des_edge* edge (size_t i) const { return w_edge.at (i); }

  // return bulk region containing element
  des_bulk* bulk () const { return w_bulk; }

private:
  des_element (size_t index, des_type type,
               const des_vertex_vector& vertex,
               const des_edge_vector& edge,
               des_bulk* bulk);

  size_t w_index;
  des_type w_type;
  des_vertex_vector w_vertex;
  des_edge_vector w_edge;
  des_bulk* w_bulk;
};



class PMI_EXTERNAL des_region {
  friend class PMIDeviceSupport;

public:
  typedef enum { bulk, contact } des_type;

  // index for region used for accessing it from the des_mesh class
  size_t index () const { return w_index; }

  // type of region
  virtual des_type type () const = 0;

  // name of region
  std::string name () const { return w_name; }

  // number of vertices in region
  size_t size_vertex () const { return w_vertex.size (); }

  // return vertex i in region
  des_vertex* vertex (size_t i) const { return w_vertex.at (i); }

  // number of edges in region
  size_t size_edge () const { return w_edge.size (); }

  // return edge i in region
  des_edge* edge (size_t i) const { return w_edge.at (i); }

  virtual ~des_region ();

protected:
  des_region (size_t index, std::string name);

private:
  void add_vertex (des_vertex* v) { w_vertex.push_back (v); }
  void add_edge (des_edge* e) { w_edge.push_back (e); }

  size_t w_index;
  std::string w_name;
  des_vertex_vector w_vertex;
  des_edge_vector w_edge;
};



class PMI_EXTERNAL des_bulk : public des_region {
  friend class PMIDeviceSupport;

public:
  // bulk region
  des_type type () const { return bulk; }

  // material of bulk region
  std::string material () const { return w_material; }

  // number of elements in region
  size_t size_element () const { return w_element.size (); }

  // return element i in region
  des_element* element (size_t i) const { return w_element.at (i); }

  // number of regioninterfaces for the region
  size_t size_regioninterface () const { return w_regioninterface.size (); }

  // return regioninterface i in region
  des_regioninterface* regioninterface (size_t i) const { return w_regioninterface.at (i); }

private:
  des_bulk (size_t index, std::string name, std::string material);
  void add_element (des_element* e) { w_element.push_back (e); }
  void add_regioninterface (des_regioninterface* ri) { w_regioninterface.push_back (ri); }

  std::string w_material;
  des_element_vector w_element;
  des_regioninterface_vector w_regioninterface;
};



class PMI_EXTERNAL des_contact : public des_region {
  friend class PMIDeviceSupport;

public:
  // contact region
  des_type type () const { return contact; }

private:
  des_contact (size_t index, std::string name);
};



class PMI_EXTERNAL des_regioninterface {
  friend class PMIDeviceSupport;

public:
  // index for region interface
  size_t index () const { return w_index; }

  // return first bulk region connected to region interface
  des_bulk* bulk1 () const { return w_bulk1; }

  // return second bulk region connected to region interface
  des_bulk* bulk2 () const { return w_bulk2; }

  // is this region interface a hetero interface?
  bool is_heterointerface () const { return w_ishetero; }

  // number of vertices in region interface
  size_t size_vertex () const { return w_vertex.size (); }

  // return vertex i in region interface
  des_vertex* vertex (size_t i) const { return w_vertex.at (i); }

  // index for region interface data
  size_t index (size_t local_vertex_index) const;

private:
  des_regioninterface (size_t index,
                       des_bulk* b1, des_bulk* b2,
                       bool ishetero, size_t offset);

  void add_vertex (des_vertex* v) { w_vertex.push_back (v); }

  size_t w_index;
  des_bulk*const w_bulk1;
  des_bulk*const w_bulk2;
  const bool w_ishetero;
  const size_t w_offset;
  des_vertex_vector w_vertex;
};



class PMI_EXTERNAL des_mesh {
  friend class PMIDeviceSupport;

public:
  // dimension of mesh (1, 2, or 3)
  int dim () const { return w_dim; }

  // return reference coordinate system:
  //
  //     [ 1  0  0 ]
  //     [ 0  1  0 ] is unified coordinate system
  //     [ 0  0  1 ]
  //     
  //     otherwise dfise coordinate system
  //
  void ref_coordinates (double ref [3][3]) const;

  // number of vertices in mesh
  size_t size_vertex () const { return w_vertex.size (); }

  // return vertex i in mesh
  des_vertex* vertex (size_t i) const { return w_vertex.at (i); }

  // number of edges in mesh
  size_t size_edge () const { return w_edge.size (); }

  // return edge i in mesh
  des_edge* edge (size_t i) const { return w_edge.at (i); }

  // number of elements in mesh
  size_t size_element () const { return w_element.size (); }

  // return element i in mesh
  des_element* element (size_t i) const { return w_element.at (i); }

  // number of regions in mesh
  size_t size_region () const { return w_region.size (); }

  // return region i in mesh
  des_region* region (size_t i) const { return w_region.at (i); }

  // number of region interfaces in mesh
  size_t size_regioninterface () const { return w_regioninterface.size (); }

  // return region interface i in mesh
  des_regioninterface* regioninterface (size_t i) const { return w_regioninterface.at (i); }

private:
  des_mesh (int dim);
  ~des_mesh ();

  void add_vertex (des_vertex* v) { w_vertex.push_back (v); }
  void add_edge (des_edge* e) { w_edge.push_back (e); }
  void add_element (des_element* e) { w_element.push_back (e); }
  void add_region (des_region* r) { w_region.push_back (r); }
  void add_regioninterface (des_regioninterface* ri) { w_regioninterface.push_back (ri); }

  int w_dim;
  double ref_coords [3][3];
  des_vertex_vector w_vertex;
  des_edge_vector w_edge;
  des_element_vector w_element;
  des_region_vector w_region;
  des_regioninterface_vector w_regioninterface;
};



class PMI_EXTERNAL des_data {
  friend class PMIDeviceSupport;

public:
  typedef enum { vertex, edge, element, rivertex } des_location;

  // element-edge coefficients: [element index][local element-edge index]
  const double*const* ReadCoefficient ();

  // element-vertex measures: [element index][local element-vertex index]
  const double*const* ReadMeasure ();

  // interface-vertex surface measures: [interface index][local interface-vertex index]
  const double*const* ReadSurfaceMeasure ();

  // return scalar data: [location index]
  const double* ReadScalar (des_location location, std::string name);

  // return vector data: [mesh dimension index][location index]
  const double*const* ReadVector (des_location location, std::string name);

  // set scalar data to values obtained from newvalue: [location index]
  void WriteScalar (des_location location, std::string name, const double* newvalue);

  // return gradient of variable 'name'
  // NOTE: Actual implementation works for vertex-based datasets only
  const double*const* ReadGradient (des_location location, std::string name);

  // return surface integral of the gradient of variable 'name' taken over the
  // boundary box divided by the box volume
  // NOTE: Actual implementation works for vertex-based datasets only
  const double* ReadFlux (des_location location, std::string name);

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
                           double*const* values) const;

  double ReadeSHEDistribution (des_bulk* r, des_vertex* v, double energy) const;
  double ReadhSHEDistribution (des_bulk* r, des_vertex* v, double energy) const;
  double ReadeSHETotalDOS (des_bulk* r, double energy) const;
  double ReadhSHETotalDOS (des_bulk* r, double energy) const;
  double ReadeSHETotalGSV (des_bulk* r, double energy) const;
  double ReadhSHETotalGSV (des_bulk* r, double energy) const;


  // Synopsys internal use only
  class des_data_identifier {
  public:
    des_location w_location;
    std::string w_name;
    des_data_identifier ();  // only used in std::map
    des_data_identifier (des_location location, std::string name);

    bool operator () (const des_data_identifier id1, const des_data_identifier id2) const;
  };

private:
  des_data (const PMI_Device_Environment& env);
  ~des_data ();

  typedef std::map <des_data_identifier, double*, des_data_identifier> des_scalar_map;
  typedef std::map <des_data_identifier, double**, des_data_identifier> des_vector_map;

  void delete_maps ();

  PMI_Device_Environment* w_env;

  double** w_coefficient;
  double** w_measure;
  double** w_surface_measure;

  des_scalar_map w_scalar_map;
  des_vector_map w_vector_map;

  des_vector_map w_grad_map;
  des_scalar_map w_flux_map;
};



class PMI_EXTERNAL PMI_Device_Common_Base {
  friend class PMI_Device_Input_Base;

protected:
  PMI_Device_Environment* pmi_env;
  PMIModel* pmi_model;

public:
  PMI_Device_Common_Base (const PMI_Device_Environment& env);

  virtual ~PMI_Device_Common_Base ();

  // run time support

  // return name of this PMI model
  const char* Name () const;

  // read parameter from Sentaurus Device parameter file
  const PMIBaseParam* ReadParameter (const char* name) const;

  // initialize parameter from Sentaurus Device parameter file or from default value
  double InitParameter (const char* name, double defaultvalue) const;

  // device mesh
  const des_mesh* Mesh () const;
};



class PMI_EXTERNAL PMI_Device_Interface : public PMI_Device_Common_Base {

public:
  PMI_Device_Interface (const PMI_Device_Environment& env);

  virtual ~PMI_Device_Interface ();

  // run time support

  // time in seconds during transient simulations
  double ReadTime () const;

  // step size in seconds during transient simulations
  double ReadTransientStepSize () const;

  // step type during transient simulations
  PMI_StepType ReadTransientStepType () const;

  // device data
  des_data* Data () const;
};



class PMI_EXTERNAL PMI_Recombination : public PMI_Vertex_Interface {

public:
  PMI_Recombination (const PMI_Environment& env);

  virtual ~PMI_Recombination ();

  // methods to be implemented by user
  virtual void Compute_r
    (const double t,        // lattice temperature
     const double n,        // electron density
     const double p,        // hole density
     const double nie,      // effective intrinsic density
     const double f,        // absolute value of electric field
     double& r) = 0;        // recombination rate

  virtual void Compute_drdt
    (const double t,        // lattice temperature
     const double n,        // electron density
     const double p,        // hole density
     const double nie,      // effective intrinsic density
     const double f,        // absolute value of electric field
     double& drdt) = 0;     // derivative of recombination rate
                            // with respect to lattice temperature

  virtual void Compute_drdn
    (const double t,        // lattice temperature
     const double n,        // electron density
     const double p,        // hole density
     const double nie,      // effective intrinsic density
     const double f,        // absolute value of electric field
     double& drdn) = 0;     // derivative of recombination rate
                            // with respect to electron density

  virtual void Compute_drdp
    (const double t,        // lattice temperature
     const double n,        // electron density
     const double p,        // hole density
     const double nie,      // effective intrinsic density
     const double f,        // absolute value of electric field
     double& drdp) = 0;     // derivative of recombination rate
                            // with respect to hole density

  virtual void Compute_drdnie
    (const double t,        // lattice temperature
     const double n,        // electron density
     const double p,        // hole density
     const double nie,      // effective intrinsic density
     const double f,        // absolute value of electric field
     double& drdnie) = 0;   // derivative of recombination rate
                            // with respect to effective intrinsic density

  virtual void Compute_drdf
    (const double t,        // lattice temperature
     const double n,        // electron density
     const double p,        // hole density
     const double nie,      // effective intrinsic density
     const double f,        // absolute value of electric field
     double& drdf) = 0;     // derivative of recombination rate
                            // with respect to absolute value of electric field
};



PMI_EXTERNAL enum PMI_AvalancheDrivingForce {
  PMI_AvalancheElectricField,
  PMI_AvalancheParallelElectricField,
  PMI_AvalancheGradQuasiFermi,
  PMI_AvalancheCarrierTemperatureCanali
};



class PMI_EXTERNAL PMI_Avalanche : public PMI_Vertex_Interface {

private:
  const PMI_AvalancheDrivingForce drivingForce;

public:
  PMI_Avalanche (const PMI_Environment& env,
                 const PMI_AvalancheDrivingForce force);

  virtual ~PMI_Avalanche ();

  PMI_AvalancheDrivingForce AvalancheDrivingForce () const { return drivingForce; }

  // methods to be implemented by user
  virtual void Compute_alpha
    (const double F,                      // driving force
     const double t,                      // lattice temperature
     const double bg,                     // bandgap
     const double ct,                     // carrier temperature
     const double currentWoMob[3],        // current density (without mobility)
     double& alpha) = 0;                  // ionization coefficient

  virtual void Compute_dalphadF
    (const double F,                      // driving force
     const double t,                      // lattice temperature
     const double bg,                     // bandgap
     const double ct,                     // carrier temperature
     const double currentWoMob[3],        // current density (without mobility)
     double& dalphadF) = 0;               // derivative of ionization coefficient
                                          // with respect to driving force

  virtual void Compute_dalphadt
    (const double F,                      // driving force
     const double t,                      // lattice temperature
     const double bg,                     // bandgap
     const double ct,                     // carrier temperature
     const double currentWoMob[3],        // current density (without mobility)
     double& dalphadt) = 0;               // derivative of ionization coefficient
                                          // with respect to lattice temperature

  virtual void Compute_dalphadbg
    (const double F,                      // driving force
     const double t,                      // lattice temperature
     const double bg,                     // bandgap
     const double ct,                     // carrier temperature
     const double currentWoMob[3],        // current density (without mobility)
     double& dalphadbg) = 0;              // derivative of ionization coefficient
                                          // with respect to bandgap

  virtual void Compute_dalphadct
    (const double F,                      // driving force
     const double t,                      // lattice temperature
     const double bg,                     // bandgap
     const double ct,                     // carrier temperature
     const double currentWoMob[3],        // current density (without mobility)
     double& dalphadct) = 0;              // derivative of ionization coefficient
                                          // with respect to carrier temperature

  virtual void Compute_dalphadcurrentWoMob
    (const double F,                      // driving force
     const double t,                      // lattice temperature
     const double bg,                     // bandgap
     const double ct,                     // carrier temperature
     const double currentWoMob[3],        // current density (without mobility)
     double dalphadcurrentWoMob[3]) = 0;  // derivative of ionization coefficient
                                          // with respect to currentWoMob
};



PMI_EXTERNAL enum PMI_AnisotropyType {
  PMI_Isotropic,
  PMI_Anisotropic
};



class PMI_EXTERNAL PMI_DopingDepMobility : public PMI_Vertex_Interface {

private:
  const PMI_AnisotropyType anisoType;

public:
  PMI_DopingDepMobility (const PMI_Environment& env,
                         const PMI_AnisotropyType anisotype);

  virtual ~PMI_DopingDepMobility ();

  PMI_AnisotropyType AnisotropyType () const { return anisoType; }

  // methods to be implemented by user
  virtual void Compute_m
    (const double n,        // electron density
     const double p,        // hole density
     const double t,        // lattice temperature
     double& m) = 0;        // doping dependent mobility

  virtual void Compute_dmdn
    (const double n,        // electron density
     const double p,        // hole density
     const double t,        // lattice temperature
     double& dmdn) = 0;     // derivative of doping dependent mobility
                            // with respect to electron density

  virtual void Compute_dmdp
    (const double n,        // electron density
     const double p,        // hole density
     const double t,        // lattice temperature
     double& dmdp) = 0;     // derivative of doping dependent mobility
                            // with respect to hole density

  virtual void Compute_dmdt
    (const double n,        // electron density
     const double p,        // hole density
     const double t,        // lattice temperature
     double& dmdt) = 0;     // derivative of doping dependent mobility
                            // with respect to lattice temperature

  // methods to optionally be overridden by the user
  virtual void Compute_dmdNa
    (const double n,        // electron density
     const double p,        // hole density
     const double t,        // lattice temperature
     double& dmdNa);        // derivative of doping dependent mobility
                            // with respect to acceptor concentration

  virtual void Compute_dmdNd
    (const double n,        // electron density
     const double p,        // hole density
     const double t,        // lattice temperature
     double& dmdNd);        // derivative of doping dependent mobility
                            // with respect to donor concentration
};



PMI_EXTERNAL enum PMI_EnormalType {
  PMI_EnormalToCurrent,
  PMI_EnormalToInterface
};



class PMI_EXTERNAL PMI_EnormalMobility : public PMI_Vertex_Interface {

private:
  const PMI_EnormalType enormalType;
  const PMI_AnisotropyType anisoType;

public:
  PMI_EnormalMobility (const PMI_Environment& env,
                       const PMI_EnormalType type,
                       const PMI_AnisotropyType anisotype);

  virtual ~PMI_EnormalMobility ();

  PMI_EnormalType EnormalType () const { return enormalType; }
  PMI_AnisotropyType AnisotropyType () const { return anisoType; }

  // methods to be implemented by user
  virtual void Compute_muinv
    (const double dist,          // distance to nearest interface
     const double pot,           // electrostatic potential
     const double enorm,         // normal electric field
     const double n,             // electron density
     const double p,             // hole density
     const double t,             // lattice temperature
     const double ct,            // carrier temperature
     double& muinv) = 0;         // inverse of mobility degradation

  virtual void Compute_dmuinvdpot
    (const double dist,          // distance to nearest interface
     const double pot,           // electrostatic potential
     const double enorm,         // normal electric field
     const double n,             // electron density
     const double p,             // hole density
     const double t,             // lattice temperature
     const double ct,            // carrier temperature
     double& dmuinvdpot) = 0;    // derivative of muinv
                                 // with respect to electrostatic potential

  virtual void Compute_dmuinvdenorm
    (const double dist,          // distance to nearest interface
     const double pot,           // electrostatic potential
     const double enorm,         // normal electric field
     const double n,             // electron density
     const double p,             // hole density
     const double t,             // lattice temperature
     const double ct,            // carrier temperature
     double& dmuinvdenorm) = 0;  // derivative of muinv
                                 // with respect to normal electric field

  virtual void Compute_dmuinvdn
    (const double dist,          // distance to nearest interface
     const double pot,           // electrostatic potential
     const double enorm,         // normal electric field
     const double n,             // electron density
     const double p,             // hole density
     const double t,             // lattice temperature
     const double ct,            // carrier temperature
     double& dmuinvdn) = 0;      // derivative of muinv
                                 // with respect to electron density

  virtual void Compute_dmuinvdp
    (const double dist,          // distance to nearest interface
     const double pot,           // electrostatic potential
     const double enorm,         // normal electric field
     const double n,             // electron density
     const double p,             // hole density
     const double t,             // lattice temperature
     const double ct,            // carrier temperature
     double& dmuinvdp) = 0;      // derivative of muinv
                                 // with respect to hole density

  virtual void Compute_dmuinvdt
    (const double dist,          // distance to nearest interface
     const double pot,           // electrostatic potential
     const double enorm,         // normal electric field
     const double n,             // electron density
     const double p,             // hole density
     const double t,             // lattice temperature
     const double ct,            // carrier temperature
     double& dmuinvdt) = 0;      // derivative of muinv
                                 // with respect to lattice temperature

  virtual void Compute_dmuinvdct
    (const double dist,          // distance to nearest interface
     const double pot,           // electrostatic potential
     const double enorm,         // normal electric field
     const double n,             // electron density
     const double p,             // hole density
     const double t,             // lattice temperature
     const double ct,            // carrier temperature
     double& dmuinvdct) = 0;     // derivative of muinv
                                 // with respect to carrier temperature

  // methods to optionally be overridden by the user
  virtual void Compute_dmuinvdNa
    (const double dist,          // distance to nearest interface
     const double pot,           // electrostatic potential
     const double enorm,         // normal electric field
     const double n,             // electron density
     const double p,             // hole density
     const double t,             // lattice temperature
     const double ct,            // carrier temperature
     double& dmuinvdNa);         // derivative of muinv
                                 // with respect to acceptor concentration

  virtual void Compute_dmuinvdNd
    (const double dist,          // distance to nearest interface
     const double pot,           // electrostatic potential
     const double enorm,         // normal electric field
     const double n,             // electron density
     const double p,             // hole density
     const double t,             // lattice temperature
     const double ct,            // carrier temperature
     double& dmuinvdNd);         // derivative of muinv
                                 // with respect to donor concentration
};



PMI_EXTERNAL enum PMI_HighFieldDrivingForce {
  PMI_HighFieldParallelElectricField,             // electric field parallel to current
  PMI_HighFieldParallelToInterfaceElectricField,  // electric field parallel to interface
  PMI_HighFieldGradQuasiFermi                     // gradient of quasi Fermi potential
};



class PMI_EXTERNAL PMI_HighFieldMobility : public PMI_Vertex_Interface {

private:
  const PMI_HighFieldDrivingForce drivingForce;
  const PMI_AnisotropyType anisoType;

public:
  PMI_HighFieldMobility (const PMI_Environment& env,
                         const PMI_HighFieldDrivingForce force,
                         const PMI_AnisotropyType anisotype);

  virtual ~PMI_HighFieldMobility ();

  PMI_HighFieldDrivingForce HighFieldDrivingForce () const { return drivingForce; }
  PMI_AnisotropyType AnisotropyType () const { return anisoType; }

  // methods to be implemented by user
  virtual void Compute_mu
    (const double pot,        // electrostatic potential
     const double n,          // electron density
     const double p,          // hole density
     const double t,          // lattice temperature
     const double ct,         // carrier temperature
     const double mulow,      // low field mobility
     const double F,          // driving force
     double& mu) = 0;         // mobility

  virtual void Compute_dmudpot
    (const double pot,        // electrostatic potential
     const double n,          // electron density
     const double p,          // hole density
     const double t,          // lattice temperature
     const double ct,         // carrier temperature
     const double mulow,      // low field mobility
     const double F,          // driving force
     double& dmudpot) = 0;    // derivative of mobility
                              // with respect to electrostatic potential

  virtual void Compute_dmudn
    (const double pot,        // electrostatic potential
     const double n,          // electron density
     const double p,          // hole density
     const double t,          // lattice temperature
     const double ct,         // carrier temperature
     const double mulow,      // low field mobility
     const double F,          // driving force
     double& dmudn) = 0;      // derivative of mobility
                              // with respect to electron density

  virtual void Compute_dmudp
    (const double pot,        // electrostatic potential
     const double n,          // electron density
     const double p,          // hole density
     const double t,          // lattice temperature
     const double ct,         // carrier temperature
     const double mulow,      // low field mobility
     const double F,          // driving force
     double& dmudp) = 0;      // derivative of mobility
                              // with respect to hole density

  virtual void Compute_dmudt
    (const double pot,        // electrostatic potential
     const double n,          // electron density
     const double p,          // hole density
     const double t,          // lattice temperature
     const double ct,         // carrier temperature
     const double mulow,      // low field mobility
     const double F,          // driving force
     double& dmudt) = 0;      // derivative of mobility
                              // with respect to lattice temperature

  virtual void Compute_dmudct
    (const double pot,        // electrostatic potential
     const double n,          // electron density
     const double p,          // hole density
     const double t,          // lattice temperature
     const double ct,         // carrier temperature
     const double mulow,      // low field mobility
     const double F,          // driving force
     double& dmudct) = 0;     // derivative of mobility
                              // with respect to carrier temperature

  virtual void Compute_dmudmulow
    (const double pot,        // electrostatic potential
     const double n,          // electron density
     const double p,          // hole density
     const double t,          // lattice temperature
     const double ct,         // carrier temperature
     const double mulow,      // low field mobility
     const double F,          // driving force
     double& dmudmulow) = 0;  // derivative of mobility
                              // with respect to low field mobility

  virtual void Compute_dmudF
    (const double pot,        // electrostatic potential
     const double n,          // electron density
     const double p,          // hole density
     const double t,          // lattice temperature
     const double ct,         // carrier temperature
     const double mulow,      // low field mobility
     const double F,          // driving force
     double& dmudF) = 0;      // derivative of mobility
                              // with respect to driving force

  // methods to optionally be overridden by the user
  virtual void Compute_dmudNa
    (const double pot,        // electrostatic potential
     const double n,          // electron density
     const double p,          // hole density
     const double t,          // lattice temperature
     const double ct,         // carrier temperature
     const double mulow,      // low field mobility
     const double F,          // driving force
     double& dmudNa);         // derivative of mobility
                              // with respect to acceptor concentration

  virtual void Compute_dmudNd
    (const double pot,        // electrostatic potential
     const double n,          // electron density
     const double p,          // hole density
     const double t,          // lattice temperature
     const double ct,         // carrier temperature
     const double mulow,      // low field mobility
     const double F,          // driving force
     double& dmudNd);         // derivative of mobility
                              // with respect to donor concentration
};


class PMI_EXTERNAL PMI_HighFieldMobility2 : public PMI_Vertex_Interface
{
private:
  const std::string mstring;
  const int mindex;
  const PMI_AnisotropyType anisoType;

public:
  // the input data coming from the simulator
  class idata {
  public:
    idata(const void*);
    double mulow() const;     // low field mobility
    double n() const;         // electron density
    double p() const;         // hole density
    double T() const;         // lattice temperature
    double cT() const;        // carrier temperature
    double Epar() const;      // parallel electric field
    double gradQF() const;    // gradient of quasi-Fermi potential
    double EprodQF() const;   // product gradient QF and electric field
    double Na0() const;       // acceptor concentration
    double Nd0() const;       // donor concentration
  private:
    const void* w;
  };

  // the results computed by the PMI
  class odata {
  public:
    odata(void*);
    double& val();            // mobility
    double& dval_dmulow();    // derivative wrt. low field mobility
    double& dval_dn();        // derivative wrt. electron density
    double& dval_dp();        // derivative wrt. hole density
    double& dval_dT();        // derivative wrt. lattice temperature
    double& dval_dcT();       // derivative wrt. carrier temperature
    double& dval_dEpar();     // derivative wrt. parallel electric field
    double& dval_dgradQF();   // derivative wrt. gradient of quasi-Fermi potential
    double& dval_dEprodQF();  // derivative wrt. product gradient QF and field
    double& dval_dNa0();      // derivative wrt. acceptor concentration
    double& dval_dNd0();      // derivative wrt. donor concentration
  private:
    void* w;
  };

  // constructor and destructor
  PMI_HighFieldMobility2(const PMI_Environment& env,
                         const int model_index,
                         const std::string& model_string,
                         const PMI_AnisotropyType anisotype);
  virtual ~PMI_HighFieldMobility2();

  PMI_AnisotropyType AnisotropyType() const { return anisoType; }
  int model_index() const { return mindex; }
  const std::string& model_string() const { return mstring; }

  // compute value and derivatives
  virtual void compute(const idata* id, odata* od ) = 0;
};



class PMI_EXTERNAL PMI_BandGap : public PMI_Vertex_Interface {

public:
  PMI_BandGap (const PMI_Environment& env);

  virtual ~PMI_BandGap ();

  // methods to be implemented by user
  virtual void Compute_bg
    (const double t,         // lattice temperature
     double& bg) = 0;        // band bap

  virtual void Compute_dbgdt
    (const double t,         // lattice temperature
     double& dbgdt) = 0;     // derivative of band bap
                             // with respect to lattice temperature
};



class PMI_EXTERNAL PMI_BandGapNarrowing : public PMI_Vertex_Interface {

public:
  PMI_BandGapNarrowing (const PMI_Environment& env);

  virtual ~PMI_BandGapNarrowing ();

  // methods to be implemented by user
  virtual void Compute_bgn
    (double& bgn) = 0;        // band bap narrowing

  // methods to optionally be overridden by the user
  virtual void Compute_dbgndNa
    (double& dbgndNa);        // derivative of band bap narrowing
                              // with respect to acceptor concentration

  virtual void Compute_dbgndNd
    (double& dbgndNd);        // derivative of band bap narrowing
                              // with respect to donor concentration
};



class PMI_EXTERNAL PMI_Affinity : public PMI_Vertex_Interface {

public:
  PMI_Affinity (const PMI_Environment& env);

  virtual ~PMI_Affinity ();

  // methods to be implemented by user
  virtual void Compute_affinity
    (const double t,              // lattice temperature
     double& affinity) = 0;       // electron affinity

  virtual void Compute_daffinitydt
    (const double t,              // lattice temperature
     double& daffinitydt) = 0;    // derivative of electron affinity
                                  // with respect to lattice temperature
};



class PMI_EXTERNAL PMI_EffectiveMass : public PMI_Vertex_Interface {

public:
  PMI_EffectiveMass (const PMI_Environment& env);

  virtual ~PMI_EffectiveMass ();

  // methods to be implemented by user
  virtual void Compute_m
    (const double t,              // lattice temperature
     const double bg,             // band gap
     double& m) = 0;              // effective mass

  virtual void Compute_dmdt
    (const double t,              // lattice temperature
     const double bg,             // band gap
     double& dmdt) = 0;           // derivative of effective mass
                                  // with respect to lattice temperature

  virtual void Compute_dmdbg
    (const double t,              // lattice temperature
     const double bg,             // band gap
     double& dmdbg) = 0;          // derivative of effective mass
                                  // with respect to band gap
};



class PMI_EXTERNAL PMI_EnergyRelaxationTime : public PMI_Vertex_Interface {

public:
  PMI_EnergyRelaxationTime (const PMI_Environment& env);

  virtual ~PMI_EnergyRelaxationTime ();

  // methods to be implemented by user
  virtual void Compute_tau
    (const double ct,            // carrier temperature
     double& tau) = 0;           // energy relaxation time

  virtual void Compute_dtaudct
    (const double ct,            // carrier temperature
     double& dtaudct) = 0;       // derivative of energy relaxation time
                                 // with respect to carrier temperature
};



PMI_EXTERNAL enum PMI_LifetimeModel {
  PMI_SRH,   // Shockley-Read-Hall
  PMI_CDL1,  // Coupled Defect Level 1
  PMI_CDL2   // Coupled Defect Level 2
};



class PMI_EXTERNAL PMI_Lifetime : public PMI_Vertex_Interface {

private:
  const PMI_LifetimeModel lifetimeModel;

public:
  PMI_Lifetime (const PMI_Environment& env,
                const PMI_LifetimeModel model);

  virtual ~PMI_Lifetime ();

  PMI_LifetimeModel LifetimeModel () const { return lifetimeModel; }

  // methods to be implemented by user
  virtual void Compute_tau
    (const double t,             // lattice temperature
     double& tau) = 0;           // lifetime

  virtual void Compute_dtaudt
    (const double t,             // lattice temperature
     double& dtaudt) = 0;        // derivative of lifetime
                                 // with respect to lattice temperature
};



class PMI_EXTERNAL PMI_ThermalConductivity : public PMI_Vertex_Interface {

private:
  const PMI_AnisotropyType anisoType;

public:
  PMI_ThermalConductivity (const PMI_Environment& env,
                           const PMI_AnisotropyType anisotype);

  virtual ~PMI_ThermalConductivity ();

  PMI_AnisotropyType AnisotropyType () const { return anisoType; }

  // methods to be implemented by user
  virtual void Compute_kappa
    (const double t,             // lattice temperature
     double& kappa) = 0;         // thermal conductivity

  virtual void Compute_dkappadt
    (const double t,             // lattice temperature
     double& dkappadt) = 0;      // derivative of thermal conductivity
                                 // with respect to lattice temperature
};


//  KV  08/18/09
class PMI_EXTERNAL PMI_MetalResistivity : public PMI_Vertex_Interface {

public:
  PMI_MetalResistivity (const PMI_Environment& env);

  virtual ~PMI_MetalResistivity ();

  // methods to be implemented by user
  virtual void Compute_Resist
    (const double t,             // lattice temperature
     double& resist) = 0;        // metal resistivity

  virtual void Compute_dResistdt
    (const double t,             // lattice temperature
     double& dresistdt) = 0;     // derivative of thermal conductivity
                                 // with respect to lattice temperature
};



class PMI_EXTERNAL PMI_HeatCapacity : public PMI_Vertex_Interface {

public:
  PMI_HeatCapacity (const PMI_Environment& env);

  virtual ~PMI_HeatCapacity ();

  // methods to be implemented by user
  virtual void Compute_c
    (const double t,             // lattice temperature
     double& c) = 0;             // heat capacity

  virtual void Compute_dcdt
    (const double t,             // lattice temperature
     double& dcdt) = 0;          // derivative of heat capacity
                                 // with respect to lattice temperature
};

class PMI_EXTERNAL PMI_MSC_HeatCapacity : public PMI_MSC_Vertex_Interface
{
public:
  // the input data coming from the simulator
  class idata {
  public:
    idata(const NS_PMI_MSC::idata*);
    double n () const;
    double p () const;
    double T () const;
    double eT () const;
    double hT () const;
    double s ( size_t ind ) const;
  private:
    const NS_PMI_MSC::idata* w;
  };

  // the output data written by the PMI
  class odata {
  public:
    odata(NS_PMI_MSC::odata*);
    double& val ();
    double& dval_dn ();
    double& dval_dp ();
    double& dval_dT ();
    double& dval_deT ();
    double& dval_dhT ();
    double& dval_ds ( size_t ind );
  private:
    NS_PMI_MSC::odata* w;
  };


  // constructor and destructor
  PMI_MSC_HeatCapacity (const PMI_Environment& env,
                        const std::string& msconfig_name,
                        const int model_index,
                        const std::string& model_string);
  virtual ~PMI_MSC_HeatCapacity ();

  // *** methods to be implemented by the user

  // compute value and derivatives
  virtual void compute (const idata* id,
                        odata* od ) = 0;
};


class PMI_EXTERNAL PMI_MSC_ThermalConductivity : public PMI_MSC_Vertex_Interface
{
private:
  const PMI_AnisotropyType anisoType;

public:
  // the input data coming from the simulator
  class idata {
  public:
    idata(const NS_PMI_MSC::idata*);
    double n () const;             // electron density
    double p () const;             // hole density
    double T () const;             // lattice temperature
    double eT () const;            // electron temperature
    double hT () const;            // electron temperature
    double s ( size_t ind ) const; // phase fraction
  private:
    const NS_PMI_MSC::idata* w;
  };

  // the results computed by the PMI
  class odata {
  public:
    odata(NS_PMI_MSC::odata*);
    double& val ();                  // thermal conductivity
    double& dval_dn ();              // derivative wrt. electron density
    double& dval_dp ();              // derivative wrt. hole density
    double& dval_dT ();              // derivative wrt. lattice temperature
    double& dval_deT ();             // derivative wrt. electron temperature
    double& dval_dhT ();             // derivative wrt. electron temperature
    double& dval_ds ( size_t ind );  // derivative wrt. phase fraction
  private:
    NS_PMI_MSC::odata* w;
  };


  // constructor and destructor
  PMI_MSC_ThermalConductivity(const PMI_Environment& env,
                              const std::string& msconfig_name,
                              const int model_index,
                              const std::string& model_string,
                              const PMI_AnisotropyType anisotype);
  virtual ~PMI_MSC_ThermalConductivity ();

  PMI_AnisotropyType AnisotropyType () const { return anisoType; }

  // compute value and derivatives
  virtual void compute(const idata* id,
                       odata* od ) = 0;
};


class PMI_EXTERNAL PMI_MSC_Mobility : public PMI_MSC_Vertex_Interface
{
private:
  const PMI_AnisotropyType anisoType;

public:
  // the input data coming from the simulator
  class idata {
  public:
    idata(const NS_PMI_MSC::idata*);
    double n () const;             // electron density
    double p () const;             // hole density
    double T () const;             // lattice temperature
    double eT () const;            // electron temperature
    double hT () const;            // electron temperature
    double s ( size_t ind ) const; // phase fraction
  private:
    const NS_PMI_MSC::idata* w;
  };

  // the results computed by the PMI
  class odata {
  public:
    odata(NS_PMI_MSC::odata*);
    double& val ();                  // mobility
    double& dval_dn ();              // derivative wrt. electron density
    double& dval_dp ();              // derivative wrt. hole density
    double& dval_dT ();              // derivative wrt. lattice temperature
    double& dval_deT ();             // derivative wrt. electron temperature
    double& dval_dhT ();             // derivative wrt. electron temperature
    double& dval_ds ( size_t ind );  // derivative wrt. phase fraction
  private:
    NS_PMI_MSC::odata* w;
  };


  // constructor and destructor
  PMI_MSC_Mobility(const PMI_Environment& env,
                   const std::string& msconfig_name,
                   const int model_index,
                   const std::string& model_string,
                   const PMI_AnisotropyType anisotype);
  virtual ~PMI_MSC_Mobility ();

  PMI_AnisotropyType AnisotropyType () const { return anisoType; }

  // compute value and derivatives
  virtual void compute(const idata* id,
                       odata* od ) = 0;
};


class PMI_EXTERNAL PMI_Absorption : public PMI_Vertex_Interface {

public:
  PMI_Absorption (const PMI_Environment& env);

  virtual ~PMI_Absorption ();

  // methods to be implemented by user
  virtual void Compute_alpha
    (const double energy,       // optical wave energy
     const double t,            // device or lattice temp
     double& alpha) = 0;        // optical absorption coefficient
};



class PMI_EXTERNAL PMI_RefractiveIndex : public PMI_Vertex_Interface {

public:
  PMI_RefractiveIndex (const PMI_Environment& env);

  virtual ~PMI_RefractiveIndex ();

  // methods to be implemented by user
  virtual void Compute_Refract
    (const double energy,       // optical wave energy
     const double t,            // device or lattice temp
     double& Refract) = 0;      // optical Refractive Index
};


class PMI_EXTERNAL PMI_RayTraceBoundary : public PMI_Vertex_Interface {

public:
  PMI_RayTraceBoundary (const PMI_Environment& env);
  virtual ~PMI_RayTraceBoundary();

  // methods to be implemented by user
  virtual void Compute_BoundaryParameters
       (// Non-changeable quantities
        const double wavelength,                // wavelength [cm]
        const double incident_angle,            // incident angle
        const double* incident_dirvec,          // direction vector of incident ray
        const double* polarvec,                 // polarization vector of incident ray
        const double* normalvec,                // normal to surface of impingment
        const double* intersectpoint,           // intersection point
        const char* region1_name,               // name of region 1
        const char* region2_name,               // name of region 2
        const double n1_real,                   // real part of refractive index 1
        const double n1_imag,                   // imag part of refractive index 1
        const double n2_real,                   // real part of refractive index 2
        const double n2_imag,                   // imag part of refractive index 2
        // User changeable quantities
        bool& is_reflectedangle_changed,        // is reflected angle changed by user?
        bool& is_reflecteddirvec_changed,       // is reflected direction changed?
        bool& is_transmittedangle_changed,      // is transmitted angled changed?
        bool& is_transmitteddirvec_changed,     // is transmitted direction changed?
        bool& is_reflected_new_startposition,   // is reflected start position changed?
        bool& is_transmitted_new_startposition, // is transmitted start position changed?
        double& reflected_angle,                // reflected angle
        double& transmitted_angle,              // transmitted angle
        double* reflected_dirvec,               // direction vec of reflected ray
        double* transmitted_dirvec,             // direction vec of transmitted ray
        double* reflected_startposition,        // starting position of reflected ray
        double* transmitted_startposition,      // starting position of transmitted ray
        double& R_TE,                           // power TE reflection coeff.
        double& T_TE,                           // power TE transmission coeff.
        double& R_TM,                           // power TM reflection coeff.
        double& T_TM                            // power TM transmission coeff.
       ) = 0;

  // Auxiliary functions for users
  void ReadComplexRefractiveIndex(std::string materialname,
                                  double wavelength,    // in microns
                                  double& n,
                                  double& k);
private:
  const PMI_Environment* thisenv;
};



class PMI_EXTERNAL PMI_Stress : public PMI_Vertex_Interface {

public:
  PMI_Stress (const PMI_Environment& env);

  virtual ~PMI_Stress ();

  // methods to be implemented by user
  virtual void Compute_StressXX
    (double& stress_xx) = 0;          // xx component of stress

  virtual void Compute_StressYY
    (double& stress_yy) = 0;          // yy component of stress

  virtual void Compute_StressZZ
    (double& stress_zz) = 0;          // zz component of stress

  virtual void Compute_StressYZ
    (double& stress_yz) = 0;          // yz component of stress

  virtual void Compute_StressXZ
    (double& stress_xz) = 0;          // xz component of stress

  virtual void Compute_StressXY
    (double& stress_xy) = 0;          // xy component of stress
};



class PMI_EXTERNAL PMI_SpaceFactor : public PMI_Vertex_Interface {

public:
  PMI_SpaceFactor (const PMI_Environment& env);

  virtual ~PMI_SpaceFactor ();

  // methods to be implemented by user
  virtual void Compute_spacefactor
    (double& spacefactor) = 0;        // trap space factor
};



class PMI_EXTERNAL PMI_Polarization : public PMI_Vertex_Interface {

public:
  PMI_Polarization (const PMI_Environment& env);

  virtual ~PMI_Polarization ();

  // methods to be implemented by user
  virtual void Compute_pol
    (double pol [3]) = 0;          // piezoelectric polarization
};



class PMI_EXTERNAL PMI_StimEmissionCoeff : public PMI_Vertex_Interface {

public:
  PMI_StimEmissionCoeff (const PMI_Environment& env);

  virtual ~PMI_StimEmissionCoeff ();

  // methods to be implemented by user
  virtual void Compute_rstim
  (double E,                    // transition energy
   double n,                    // electron density
   double p,                    // hole density
   double et,                   // electron temperature if hydro, lattice temperature otherwise
   double ht,                   // hole temperature if hydro, lattice temperature otherwise
   double& rstim) = 0;          // stimulated emission coefficient

  virtual void Compute_drstimdn
  (double E,                    // transition energy
   double n,                    // electron density
   double p,                    // hole density
   double et,                   // electron temperature if hydro, lattice temperature otherwise
   double ht,                   // hole temperature if hydro, lattice temperature otherwise
   double& drstimdn) = 0;       // derivative of stimulated emission coefficient
                                // with respect to electron density

  virtual void Compute_drstimdp
  (double E,                    // transition energy
   double n,                    // electron density
   double p,                    // hole density
   double et,                   // electron temperature if hydro, lattice temperature otherwise
   double ht,                   // hole temperature if hydro, lattice temperature otherwise
   double& drstimdp) = 0;       // derivative of stimulated emission coefficient
                                // with respect to hole density

  virtual void Compute_drstimdet
  (double E,                    // transition energy
   double n,                    // electron density
   double p,                    // hole density
   double et,                   // electron temperature if hydro, lattice temperature otherwise
   double ht,                   // hole temperature if hydro, lattice temperature otherwise
   double& drstimdet) = 0;      // derivative of stimulated emission coefficient
                                // with respect to electron or lattice temperature

  virtual void Compute_drstimdht
  (double E,                    // transition energy
   double n,                    // electron density
   double p,                    // hole density
   double et,                   // electron temperature if hydro, lattice temperature otherwise
   double ht,                   // hole temperature if hydro, lattice temperature otherwise
   double& drstimdht) = 0;      // derivative of stimulated emission coefficient
                                // with respect to hole or lattice temperature
};



class PMI_EXTERNAL PMI_PhotonPhaseCoeff : public PMI_Vertex_Interface {

public:
  PMI_PhotonPhaseCoeff (const PMI_Environment& env);

  virtual ~PMI_PhotonPhaseCoeff ();

  // methods to be implemented by user
  virtual void Compute_pphase
  (double E,                    // transition energy
   double n,                    // electron density
   double p,                    // hole density
   double et,                   // electron temperature if hydro, lattice temperature otherwise
   double ht,                   // hole temperature if hydro, lattice temperature otherwise
   double& pphase) = 0;         // photon phase coefficient

  virtual void Compute_dpphasedn
  (double E,                    // transition energy
   double n,                    // electron density
   double p,                    // hole density
   double et,                   // electron temperature if hydro, lattice temperature otherwise
   double ht,                   // hole temperature if hydro, lattice temperature otherwise
   double& dpphasedn) = 0;      // derivative of photon phase coefficient
                                // with respect to electron density

  virtual void Compute_dpphasedp
  (double E,                    // transition energy
   double n,                    // electron density
   double p,                    // hole density
   double et,                   // electron temperature if hydro, lattice temperature otherwise
   double ht,                   // hole temperature if hydro, lattice temperature otherwise
   double& dpphasedp) = 0;      // derivative of photon phase coefficient
                                // with respect to hole density

  virtual void Compute_dpphasedet
  (double E,                    // transition energy
   double n,                    // electron density
   double p,                    // hole density
   double et,                   // electron temperature if hydro, lattice temperature otherwise
   double ht,                   // hole temperature if hydro, lattice temperature otherwise
   double& dpphasedet) = 0;     // derivative of photon phase coefficient
                                // with respect to electron or lattice temperature

  virtual void Compute_dpphasedht
  (double E,                    // transition energy
   double n,                    // electron density
   double p,                    // hole density
   double et,                   // electron temperature if hydro, lattice temperature otherwise
   double ht,                   // hole temperature if hydro, lattice temperature otherwise
   double& dpphasedht) = 0;     // derivative of photon phase coefficient
                                // with respect to hole or lattice temperature
};



class PMI_EXTERNAL PMI_SponEmissionCoeff : public PMI_Vertex_Interface {

public:
  PMI_SponEmissionCoeff (const PMI_Environment& env);

  virtual ~PMI_SponEmissionCoeff ();

  // methods to be implemented by user
  virtual void Compute_rspon
  (double E,                    // transition energy
   double n,                    // electron density
   double p,                    // hole density
   double et,                   // electron temperature if hydro, lattice temperature otherwise
   double ht,                   // hole temperature if hydro, lattice temperature otherwise
   double& rspon) = 0;          // spontaneous emission coefficient

  virtual void Compute_drspondn
  (double E,                    // transition energy
   double n,                    // electron density
   double p,                    // hole density
   double et,                   // electron temperature if hydro, lattice temperature otherwise
   double ht,                   // hole temperature if hydro, lattice temperature otherwise
   double& drspondn) = 0;       // derivative of spontaneous emission coefficient
                                // with respect to electron density

  virtual void Compute_drspondp
  (double E,                    // transition energy
   double n,                    // electron density
   double p,                    // hole density
   double et,                   // electron temperature if hydro, lattice temperature otherwise
   double ht,                   // hole temperature if hydro, lattice temperature otherwise
   double& drspondp) = 0;       // derivative of spontaneous emission coefficient
                                // with respect to hole density

  virtual void Compute_drspondet
  (double E,                    // transition energy
   double n,                    // electron density
   double p,                    // hole density
   double et,                   // electron temperature if hydro, lattice temperature otherwise
   double ht,                   // hole temperature if hydro, lattice temperature otherwise
   double& drspondet) = 0;      // derivative of spontaneous emission coefficient
                                // with respect to electron or lattice temperature

  virtual void Compute_drspondht
  (double E,                    // transition energy
   double n,                    // electron density
   double p,                    // hole density
   double et,                   // electron temperature if hydro, lattice temperature otherwise
   double ht,                   // hole temperature if hydro, lattice temperature otherwise
   double& drspondht) = 0;      // derivative of spontaneous emission coefficient
                                // with respect to hole or lattice temperature
};



PMI_EXTERNAL enum PMI_SpeciesType {
  PMI_acceptor,
  PMI_donor
};

class PMI_EXTERNAL PMI_DistributionFunction : public PMI_Vertex_Interface {
//
//            NdPlus = Nd*f(T,Ef),
//  where
//            f(T,Ef) = 1/(1 +g(T)*exp((Ef - Ec)/kT) )
//   or
//
//            NaMinus = Na*f(T,Ef),
//  where
//            f(T,Ef) = 1/(1 +g(T)*exp((Ev - Ef)/kT) )
//

private:
  const PMI_SpeciesType speciesType;
  const char* speciesName;

public:
  PMI_DistributionFunction (const PMI_Environment& env,
                            const char* name,
                            const PMI_SpeciesType type = PMI_acceptor);

  virtual ~PMI_DistributionFunction ();

  PMI_SpeciesType SpeciesType () const { return speciesType; }
  const char* SpeciesName () const { return speciesName; }

  // read parameter from Sentaurus Device parameter file
  // (override for PMI_Vertex_Interface::ReadParameter)
  const PMIBaseParam* ReadParameter (const char* name) const;

  // initialize parameter from Sentaurus Device parameter file or from default value
  // (override for PMI_Vertex_Interface::InitParameter)
  double InitParameter (const char* name, double defaultvalue) const;

  // methods to be implemented by user
  virtual void Compute_g
    (const double t,        // lattice temperature
     double& g) = 0;        // g = g(t)

  virtual void Compute_dgdt
    (const double t,        // lattice temperature
     double& dgdt) = 0;     // dgdt = g'(t)
};


class PMI_EXTERNAL PMI_SpatialDistributionFunction: public PMI_Vertex_Interface {
//  * The spatial distribution function:
//  *
//  *   R(w,l,E)
//  *
//  * where
//  *   - l is the coordinate along the particle path [um];
//  *   - w is radial coordinate orthogonal to l [um];
//  *   - E is the kinematically particle energy;

private:
  const char* HeavyIonType;

public:
  PMI_SpatialDistributionFunction (const PMI_Environment& env, const char* name);

  virtual ~PMI_SpatialDistributionFunction ();

  const char* GetHeavyIonType () const { return HeavyIonType; }

  // methods to be implemented by user
  virtual void Compute_R(double& R, const double w, const double l = -1., const double E = -1.) = 0;

};


class PMI_EXTERNAL PMI_TrapCaptureEmission : public PMI_Vertex_Interface {
public:
    PMI_TrapCaptureEmission(const PMI_Environment& env);

    virtual ~PMI_TrapCaptureEmission();

    // methods to be implemented by user
    virtual void Compute_rates
    (const double n,            // electron density
     const double p,            // hole density
     const double t,            // lattice temperature
     const double tn,           // electron temperature
     const double tp,           // hole temperature
     const double f,            // absolute value of electric field
     double& capture,           // capture rate
     double& emission) = 0;     // emission rate

    virtual void Compute_dratesdn
    (const double n,            // electron density
     const double p,            // hole density
     const double t,            // lattice temperature
     const double tn,           // electron temperature
     const double tp,           // hole temperature
     const double f,            // absolute value of electric field
     double& dcapturedn,        // derivative of capture rate wrt electron density
     double& demissiondn) = 0;  // derivative of emission ratewrt electron density

    virtual void Compute_dratesdp
    (const double n,            // electron density
     const double p,            // hole density
     const double t,            // lattice temperature
     const double tn,           // electron temperature
     const double tp,           // hole temperature
     const double f,            // absolute value of electric field
     double& dcapturedp,        // derivative of capture rate wrt hole density
     double& demissiondp) = 0;  // derivative of emission ratewrt hole density

    virtual void Compute_dratesdt
    (const double n,            // electron density
     const double p,            // hole density
     const double t,            // lattice temperature
     const double tn,           // electron temperature
     const double tp,           // hole temperature
     const double f,            // absolute value of electric field
     double& dcapturedt,        // derivative of capture rate wrt lattice temperature
     double& demissiondt) = 0;  // derivative of emission ratewrt lattice temperature

    virtual void Compute_dratesdtn
    (const double n,            // electron density
     const double p,            // hole density
     const double t,            // lattice temperature
     const double tn,           // electron temperature
     const double tp,           // hole temperature
     const double f,            // absolute value of electric field
     double& dcapturedtn,       // derivative of capture rate wrt electron temperature
     double& demissiondtn) = 0; // derivative of emission ratewrt electron temperature

    virtual void Compute_dratesdtp
    (const double n,            // electron density
     const double p,            // hole density
     const double t,            // lattice temperature
     const double tn,           // electron temperature
     const double tp,           // hole temperature
     const double f,            // absolute value of electric field
     double& dcapturedtp,       // derivative of capture rate wrt hole temperature
     double& demissiondtp) = 0; // derivative of emission ratewrt hole temperature

    virtual void Compute_dratesdf
    (const double n,            // electron density
     const double p,            // hole density
     const double t,            // lattice temperature
     const double tn,           // electron temperature
     const double tp,           // hole temperature
     const double f,            // absolute value of electric field
     double& dcapturedf,        // derivative of capture rate wrt electric field
     double& demissiondf) = 0;  // derivative of emission ratewrt emission field
};


class PMI_EXTERNAL PMI_TrapEnergyShift : public PMI_Vertex_Interface {
public:
    PMI_TrapEnergyShift(const PMI_Environment& env);

    virtual ~PMI_TrapEnergyShift();

    // methods to be implemented by the user
    virtual void Compute_shift(
        const double f[3],      // electric field vector
        const double t,         // lattice temperature
        double& shift) = 0;     // trap energy shift
    virtual void Compute_dshiftdf(
        const double f[3],      // electric field vector
        const double t,         // lattice temperature
        double df[3]) = 0;      // derivative of energy shift wrt field components
    virtual void Compute_dshiftdt(
        const double f[3],      // electric field vector
        const double t,         // lattice temperature
        double& dt) = 0;        // derivative of energy shift wrt lattice temperature
};

class PMI_EXTERNAL PMI_ApparentBandEdgeShift : public PMI_Vertex_Interface {

public:
  PMI_ApparentBandEdgeShift (const PMI_Environment& env);

  virtual ~PMI_ApparentBandEdgeShift ();

  // methods to be implemented by user
  virtual void Compute_shift
    (const double n,            // electron density
     const double p,            // hole density
     const double t,            // lattice temperature
     const double f,            // absolute value of electric field
     double& shift) = 0;        // apparent band edge shift

  virtual void Compute_dshiftdn
    (const double n,            // electron density
     const double p,            // hole density
     const double t,            // lattice temperature
     const double f,            // absolute value of electric field
     double& dshiftdn) = 0;     // derivative of shift with respect to electron density

  virtual void Compute_dshiftdp
    (const double n,            // electron density
     const double p,            // hole density
     const double t,            // lattice temperature
     const double f,            // absolute value of electric field
     double& dshiftdp) = 0;     // derivative of shift with respect to hole density

  virtual void Compute_dshiftdt
    (const double n,            // electron density
     const double p,            // hole density
     const double t,            // lattice temperature
     const double f,            // absolute value of electric field
     double& dshiftdt) = 0;     // derivative with respect to lattice temperature

  virtual void Compute_dshiftdf
    (const double n,            // electron density
     const double p,            // hole density
     const double t,            // lattice temperature
     const double f,            // absolute value of electric field
     double& dshiftdf) = 0;     // derivative with respect to electric field
};

//======================================================================

// though the model is derived from PMI_MSC_Vertex_Interface
// it uses its original interface (by using e_var, input_data, output_data)
//
//  the classes idata/odata are defined only for internal reasons
//  and might build the interface in future versions

class PMI_EXTERNAL PMI_MSC_ApparentBandEdgeShift
    : public PMI_MSC_Vertex_Interface
{
public:
  // variables
  enum e_var {
        var_n,        // electron density
        var_p,        // hole density
        var_T,        // lattice temperature
        var_eT,       // electron temperature
        var_hT,       // hole temperature
        var_s,        // occupation rates of msconfig
        var_undefined // last entry
  };

  class input_data
  {
  public: // will be made private with the next version
          // please use the access functions val() and the
          // short variants n(), p(), etc.

    double w_n;         // electron density [ cm^{-1} ]
    double w_p;         // hole density [ cm^{-1} ]
    double w_T;         // lattice temperature [K]
    double w_eT;        // electron temperature [K]
    double w_hT;        // hole temperature [K]
    std::vector<double> w_s; // occupation rates of reference msconfig [1]

  public:
    input_data() : w_s(0) { w_n = w_p = w_T = w_eT = w_hT = 0.; }
    ~input_data() {}

    double& val ( e_var var, size_t ind )
    { switch ( var ) {
      case var_n: return w_n;
      case var_p: return w_p;
      case var_T: return w_T;
      case var_eT: return w_eT;
      case var_hT: return w_hT;
      case var_s: return w_s[ind];
      case var_undefined:
        std::cout << std::endl
                << "PMI_MSC_ApparentBandEdgeShift:: val: undefined var (1)";
        exit (1);
      }
      return w_s[1000]; // should cause an error
    }

    double val ( e_var var, size_t ind ) const
    { switch ( var ) {
      case var_n: return w_n;
      case var_p: return w_p;
      case var_T: return w_T;
      case var_eT: return w_eT;
      case var_hT: return w_hT;
      case var_s: return w_s[ind];
      case var_undefined:
        std::cout << std::endl
                << "PMI_MSC_ApparentBandEdgeShift:: val: undefined var (2)";
        exit (1);
      }
      return 1.e300; // should cause an error
    }

    double& n() { return w_n; }
    double  n() const { return w_n; }

    double& p() { return w_p; }
    double  p() const { return w_p; }

    double& T() { return w_T; }
    double  T() const { return w_T; }

    double& eT() { return w_eT; }
    double  eT() const { return w_eT; }

    double& hT() { return w_hT; }
    double  hT() const { return w_hT; }

    std::vector<double>& s() { return w_s; }
    const std::vector<double>& s() const { return w_s; }

    double& s(size_t ind) { return w_s[ind]; }
    double  s(size_t ind) const { return w_s[ind]; }
  };

private:
  // dependency control
  // default is 'true' for all variables
  std::vector<bool> w_dependency_used;

protected:
  // the actual status described essentially by the solution variables
  // at this vertex
  //
  // the user might provide its own control if the status has been
  // updated by overloading the set_actual_status function
  // and counting its calls
  input_data w_actual_status;

  // For backward compartibility only.  Do not use.
  int w_model_index;

public:

  // sets the actual status
  // the function might be overloaded to perform 'internal' computations
  //
  // the function is called once before the calls to the
  // compute_val and compute_dval_dX functions
  virtual void set_actual_status (
    const PMI_MSC_ApparentBandEdgeShift::input_data& id )
    { w_actual_status = id; }

public:
  PMI_MSC_ApparentBandEdgeShift ( const PMI_Environment& env,
    const std::string& msconfig_name, int model_index = 0,
    const std::string& model_string = std::string ( "" ) );
  virtual ~PMI_MSC_ApparentBandEdgeShift();

  // functions which should be provided by the user
  // if the corresponding variable is used (see set_dependency_used)
  virtual void compute_val ( double& val ) = 0;
  virtual void compute_dval_dn ( double& val );
  virtual void compute_dval_dp ( double& val );
  virtual void compute_dval_dT ( double& val );
  virtual void compute_dval_deT ( double& val );
  virtual void compute_dval_dhT ( double& val );
  virtual void compute_dval_ds ( std::vector<double>& val );

  // dependency control
  // for unused dependencies the corresponding compute_dval_dX
  // function is never called
  bool dependency_used ( PMI_MSC_ApparentBandEdgeShift::e_var var ) const
    { return w_dependency_used[var]; }

  // the user might set which dependencies are really used
  void set_dependency_used ( PMI_MSC_ApparentBandEdgeShift::e_var var,
    bool flag ) { w_dependency_used[var] = flag; }
};



class PMI_EXTERNAL PMI_CurrentPlot : public PMI_Device_Interface {

public:
  PMI_CurrentPlot (const PMI_Device_Environment& env);

  virtual ~PMI_CurrentPlot ();

  // methods to be implemented by user
  virtual void Compute_Dataset_Names
    (des_string_vector& dataset) = 0;   // array of dataset names

  virtual void Compute_Function_Names
    (des_string_vector& function) = 0;  // array of function names

  virtual void Compute_Plot_Values
    (des_double_vector& value) = 0;     // array of plot values
};



PMI_EXTERNAL enum PMI_CarrierType {
  PMI_Hole,
  PMI_Electron
};


class PMI_EXTERNAL PMI_HotCarrierInjection : public PMI_Device_Interface {

protected:
  const PMI_CarrierType cType;

public:
  PMI_HotCarrierInjection (const PMI_Device_Environment& env,
                           const PMI_CarrierType cType);

  virtual ~PMI_HotCarrierInjection ();

  // methods to be implemented by user
  virtual void Compute_gCurr
    (const des_regioninterface_vector& regioninterfaces, // region interfaces
                                                         // associated with model name
     des_array_vector& gCurr) = 0;                       // gate injection current in
                                                         // each vertex of specified
                                                         // region interfaces
};



class PMI_EXTERNAL PMI_PiezoresistanceFactor : public PMI_Vertex_Interface
{
public:
  PMI_PiezoresistanceFactor (const PMI_Environment& env);

  virtual ~PMI_PiezoresistanceFactor ();

  // methods to be implemented by user
  virtual bool IsPrefactor() = 0;

  virtual void Compute_Pij
    (const double Enormal,
     double& p11,
     double& p12,
     double& p44) = 0;

  virtual void Compute_DerPij
    (const double Enormal,
     double& dp11,
     double& dp12,
     double& dp44) = 0;
};



class PMI_EXTERNAL PMI_PostProcess : public PMI_Device_Interface {

public:
  PMI_PostProcess (const PMI_Device_Environment& env);

  virtual ~PMI_PostProcess ();

  // methods to be implemented by user
  virtual void Compute_PostProcess () = 0;

};



extern "C" {

typedef PMI_PiezoresistanceFactor* new_PMI_PiezoresistanceFactor_func
  (const PMI_Environment& env);

typedef PMI_DistributionFunction* new_PMI_DistributionFunction_func
  (const PMI_Environment& env, const char* name, const PMI_SpeciesType type);

typedef PMI_SpatialDistributionFunction* new_PMI_SpatialDistributionFunction_func
  (const PMI_Environment& env, const char* HeavyIonName);

typedef PMI_Recombination* new_PMI_Recombination_func
  (const PMI_Environment& env);

typedef int PMI_Recombination_ElectricField_func ();

typedef PMI_Avalanche* new_PMI_Avalanche_func
  (const PMI_Environment& env, const PMI_AvalancheDrivingForce force);

typedef PMI_DopingDepMobility* new_PMI_DopingDepMobility_func
  (const PMI_Environment& env, const PMI_AnisotropyType anisotype);

typedef PMI_EnormalMobility* new_PMI_EnormalMobility_func
  (const PMI_Environment& env, const PMI_EnormalType type,
   const PMI_AnisotropyType anisotype);

typedef PMI_HighFieldMobility* new_PMI_HighFieldMobility_func
  (const PMI_Environment& env, const PMI_HighFieldDrivingForce force,
   const PMI_AnisotropyType anisotype);

typedef PMI_HighFieldMobility2* new_PMI_HighFieldMobility2_func
  (const PMI_Environment& env, int model_index,
   const std::string& model_string,
   const PMI_AnisotropyType anisotype);

typedef PMI_BandGap* new_PMI_BandGap_func
  (const PMI_Environment& env);

typedef PMI_BandGapNarrowing* new_PMI_BandGapNarrowing_func
  (const PMI_Environment& env);

typedef PMI_Affinity* new_PMI_Affinity_func
  (const PMI_Environment& env);

typedef PMI_EffectiveMass* new_PMI_EffectiveMass_func
  (const PMI_Environment& env);

typedef PMI_EnergyRelaxationTime* new_PMI_EnergyRelaxationTime_func
  (const PMI_Environment& env);

typedef PMI_Lifetime* new_PMI_Lifetime_func
  (const PMI_Environment& env, const PMI_LifetimeModel model);

typedef PMI_ThermalConductivity* new_PMI_ThermalConductivity_func
  (const PMI_Environment& env, const PMI_AnisotropyType anisotype);

typedef PMI_MetalResistivity* new_PMI_MetalResistivity_func
  (const PMI_Environment& env);

typedef PMI_HeatCapacity* new_PMI_HeatCapacity_func
  (const PMI_Environment& env);

typedef PMI_MSC_HeatCapacity* new_PMI_MSC_HeatCapacity_func
  (const PMI_Environment& env, const std::string& msconfig_name,
        int model_index, const std::string& model_string );

typedef PMI_MSC_ThermalConductivity* new_PMI_MSC_ThermalConductivity_func
  (const PMI_Environment& env, const std::string& msconfig_name,
   int model_index, const std::string& model_string,
   const PMI_AnisotropyType anisotype);

typedef PMI_MSC_Mobility* new_PMI_MSC_Mobility_func
  (const PMI_Environment& env, const std::string& msconfig_name,
   int model_index, const std::string& model_string,
   const PMI_AnisotropyType anisotype);

typedef PMI_Absorption* new_PMI_Absorption_func
  (const PMI_Environment& env);

typedef PMI_RefractiveIndex* new_PMI_RefractiveIndex_func
  (const PMI_Environment& env);

typedef PMI_RayTraceBoundary* new_PMI_RayTraceBoundary_func
  (const PMI_Environment& env);

typedef PMI_Stress* new_PMI_Stress_func
  (const PMI_Environment& env);

typedef PMI_SpaceFactor* new_PMI_SpaceFactor_func
  (const PMI_Environment& env);

typedef PMI_Polarization* new_PMI_Polarization_func
  (const PMI_Environment& env);

typedef PMI_StimEmissionCoeff* new_PMI_StimEmissionCoeff_func
  (const PMI_Environment& env);

typedef PMI_PhotonPhaseCoeff* new_PMI_PhotonPhaseCoeff_func
  (const PMI_Environment& env);

typedef PMI_SponEmissionCoeff* new_PMI_SponEmissionCoeff_func
  (const PMI_Environment& env);

typedef PMI_TrapCaptureEmission* new_PMI_TrapCaptureEmission_func
  (const PMI_Environment& env, int id);

typedef PMI_TrapEnergyShift* new_PMI_TrapEnergyShift_func
  (const PMI_Environment& env, int id);

typedef PMI_ApparentBandEdgeShift* new_PMI_ApparentBandEdgeShift_func
  (const PMI_Environment& env);

typedef PMI_MSC_ApparentBandEdgeShift* new_PMI_MSC_ApparentBandEdgeShift_func
  (const PMI_Environment& env, const std::string& msconfig_name, int id);

typedef PMI_CurrentPlot* new_PMI_CurrentPlot_func
  (const PMI_Device_Environment& env);

typedef PMI_HotCarrierInjection* new_PMI_HotCarrierInjection_func
  (const PMI_Device_Environment& env, const PMI_CarrierType carType);

typedef PMI_PostProcess* new_PMI_PostProcess_func
  (const PMI_Device_Environment& env);

#ifndef PMI_DESSIS_SRC

// functions to be implemented by user

#undef PMI_EXTERNAL
#if defined(_MSC_VER)
#define PMI_EXTERNAL __declspec (dllexport)
#else
#define PMI_EXTERNAL
#endif


// virtual constructor
PMI_EXTERNAL new_PMI_Recombination_func new_PMI_Recombination;

// does the recombination model depend on the electric field?
// optional function (by default we assume a dependence on the electric field)
PMI_EXTERNAL PMI_Recombination_ElectricField_func PMI_Recombination_ElectricField;

// virtual constructor
PMI_EXTERNAL new_PMI_Avalanche_func new_PMI_e_Avalanche;
PMI_EXTERNAL new_PMI_Avalanche_func new_PMI_h_Avalanche;

// virtual constructor
PMI_EXTERNAL new_PMI_HotCarrierInjection_func new_PMI_e_HotCarrierInjection;
PMI_EXTERNAL new_PMI_HotCarrierInjection_func new_PMI_h_HotCarrierInjection;

// virtual constructor
PMI_EXTERNAL new_PMI_DopingDepMobility_func new_PMI_DopingDep_e_Mobility;
PMI_EXTERNAL new_PMI_DopingDepMobility_func new_PMI_DopingDep_h_Mobility;

// virtual constructor
PMI_EXTERNAL new_PMI_EnormalMobility_func new_PMI_Enormal_e_Mobility;
PMI_EXTERNAL new_PMI_EnormalMobility_func new_PMI_Enormal_h_Mobility;

// virtual constructor
PMI_EXTERNAL new_PMI_HighFieldMobility_func new_PMI_HighField_e_Mobility;
PMI_EXTERNAL new_PMI_HighFieldMobility_func new_PMI_HighField_h_Mobility;

// virtual constructor
PMI_EXTERNAL new_PMI_HighFieldMobility2_func new_PMI_HighFieldMobility2;

// virtual constructor
PMI_EXTERNAL new_PMI_BandGap_func new_PMI_BandGap;

// virtual constructor
PMI_EXTERNAL new_PMI_BandGapNarrowing_func new_PMI_BandGapNarrowing;

// virtual constructor
PMI_EXTERNAL new_PMI_Affinity_func new_PMI_Affinity;

// virtual constructor
PMI_EXTERNAL new_PMI_EffectiveMass_func new_PMI_e_EffectiveMass;
PMI_EXTERNAL new_PMI_EffectiveMass_func new_PMI_h_EffectiveMass;

// virtual constructor
PMI_EXTERNAL new_PMI_EnergyRelaxationTime_func new_PMI_e_EnergyRelaxationTime;
PMI_EXTERNAL new_PMI_EnergyRelaxationTime_func new_PMI_h_EnergyRelaxationTime;

// virtual constructor
PMI_EXTERNAL new_PMI_Lifetime_func new_PMI_e_Lifetime;
PMI_EXTERNAL new_PMI_Lifetime_func new_PMI_h_Lifetime;

// virtual constructor
PMI_EXTERNAL new_PMI_ThermalConductivity_func new_PMI_ThermalConductivity;

// virtual constructor
PMI_EXTERNAL new_PMI_MetalResistivity_func new_PMI_MetalResistivity;

// virtual constructor
PMI_EXTERNAL new_PMI_HeatCapacity_func new_PMI_HeatCapacity;

// virtual constructor
PMI_EXTERNAL new_PMI_MSC_HeatCapacity_func new_PMI_MSC_HeatCapacity;

// virtual constructor
PMI_EXTERNAL new_PMI_MSC_ThermalConductivity_func new_PMI_MSC_ThermalConductivity;

// virtual constructor
PMI_EXTERNAL new_PMI_MSC_Mobility_func new_PMI_MSC_Mobility;

// virtual constructor
PMI_EXTERNAL new_PMI_Absorption_func new_PMI_Absorption;

// virtual constructor
PMI_EXTERNAL new_PMI_RefractiveIndex_func new_PMI_RefractiveIndex;

// virtual constructor
PMI_EXTERNAL new_PMI_RayTraceBoundary_func new_PMI_RayTraceBoundary;

// virtual constructor
PMI_EXTERNAL new_PMI_Stress_func new_PMI_Stress;

// virtual constructor
PMI_EXTERNAL new_PMI_SpaceFactor_func new_PMI_SpaceFactor;

// virtual constructor
PMI_EXTERNAL new_PMI_Polarization_func new_PMI_Polarization;

// virtual constructor
PMI_EXTERNAL new_PMI_StimEmissionCoeff_func new_PMI_StimEmissionCoeff;

// virtual constructor
PMI_EXTERNAL new_PMI_PhotonPhaseCoeff_func new_PMI_PhotonPhaseCoeff;

// virtual constructor
PMI_EXTERNAL new_PMI_SponEmissionCoeff_func new_PMI_SponEmissionCoeff;

// virtual constructor
PMI_EXTERNAL new_PMI_PiezoresistanceFactor_func new_PMI_ePiezoresistanceFactor;
PMI_EXTERNAL new_PMI_PiezoresistanceFactor_func new_PMI_hPiezoresistanceFactor;

// virtual constructor
PMI_EXTERNAL new_PMI_DistributionFunction_func new_PMI_DistributionFunction;

// virtual constructor
PMI_EXTERNAL new_PMI_SpatialDistributionFunction_func new_PMI_SpatialDistributionFunction;

// virtual constructor
PMI_EXTERNAL new_PMI_TrapCaptureEmission_func new_PMI_TrapCaptureEmission;

// virtual constructor
PMI_EXTERNAL new_PMI_TrapEnergyShift_func new_PMI_TrapEnergyShift;

// virtual constructor
PMI_EXTERNAL new_PMI_ApparentBandEdgeShift_func new_PMI_ApparentBandEdgeShift;

// virtual constructor
PMI_EXTERNAL new_PMI_MSC_ApparentBandEdgeShift_func new_PMI_MSC_ApparentBandEdgeShift;

// virtual constructor
PMI_EXTERNAL new_PMI_CurrentPlot_func new_PMI_CurrentPlot;

// virtual constructor
PMI_EXTERNAL new_PMI_PostProcess_func new_PMI_PostProcess;

#endif

} // extern "C"

#undef PMI_EXTERNAL

#endif
