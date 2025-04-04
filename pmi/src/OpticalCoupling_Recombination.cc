/*------------------------------------------------------------------------
 * File: NonlocalRec.cc
 * Description: model for optical coupling in sentaurus device.
 *-----------------------------------------------------------------------*/ 

#include "PMIModels.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"

#include <assert.h>
#include <iostream>
using namespace std;

#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>

#define DEBUG 0

/*-------------------------------------------------------------------------

 sdevice calls the PMI functions in this order:

   During setup, for each instance:
       ::OpticalCoupling_Recombination()    (the constructor)
       ::DefineDependencies()
       ::DefineJacobians()

   Then during the solution:
       ::Compute_parallel() is called following each Newton iteration, at least once per instance.
 *-----------------------------------------------------------------------*/


class OpticalCoupling_Recombination : public PMI_NonLocal_Recombination
{
   private:
     // Model parameters.  By declaring them static, they are shared between all instances of this class.
     // This saves time and memory because LCMatrix will be huge once it's fully initialized.
     static bool initialized;
     bool is_firstinstance;
     static std::string LCMatrix_file;
     static des_jacobian* LCMatrix;
     static std::vector<double> Brad_list;
     
     //  these are not static -- different for each class instance.
     double Brad;
     int layer_id;

     void Read_LCMatrix();
     
   public:
// Constructor
     OpticalCoupling_Recombination (const PMI_Device_Environment& env); 
    ~OpticalCoupling_Recombination (); // Destructor
     
     void DefineDependencies(std::vector<des_data::des_id>& dependencies);
     void DefineJacobians(des_id_to_jacobian_map& J_elec,
                          des_id_to_jacobian_map& J_hole);
     void Compute_parallel (const PMI_NonLocal_Recombination::Input& input,
				PMI_NonLocal_Recombination::Output& output); // main calculation

// implement this method to force reevaluation of dependencies
//     bool NeedNewEdges();

};

// Static class members must be declared outside the class definition.
bool OpticalCoupling_Recombination::initialized;
std::string OpticalCoupling_Recombination::LCMatrix_file;
des_jacobian* OpticalCoupling_Recombination::LCMatrix;

std::vector<double> OpticalCoupling_Recombination::Brad_list;


/*-----------------------------------------------------------------
 * Constructor: OpticalCoupling_Recombination
 *--------------------------------------------------------------------*/
// Constructor for the class OpticalCoupling_Recombination.  Retrieve the name of the LCMatric file
// and assemble the matrix in memory.

OpticalCoupling_Recombination::OpticalCoupling_Recombination (const PMI_Device_Environment& env):PMI_NonLocal_Recombination (env)
{
	const des_mesh* mesh = Mesh ();
	const int n_vertices = mesh->size_vertex ();

	// This will be initialized in compute_parallel once we know what region this instance applies to.
	layer_id =-1;        

	// This will only run once for the first instance of the class.
        if (initialized == false)
	{
		initialized = true;
		is_firstinstance = true;

	 	// fill Brad_list with empty elements
	 	for (int i =0; i < mesh->size_region(); i++)
			Brad_list.push_back(0.0);

		cout << "** Initializing Optical Coupling data" << endl;

		LCMatrix_file = InitStringParameter("MatrixFile", "");
	        std::cout << "MatrixFile " << LCMatrix_file << endl;

		// Create a sparse matrix to contain the coupling coefficients
	        LCMatrix = new des_jacobian(n_vertices, n_vertices, 1,1);

	        // Load optical generation matrix and compute vertex measures
		Read_LCMatrix();



	} else {
		is_firstinstance=false;
               cout << "** created new OpticalCoupling class instance, class is already initialized"  << endl;
	}
}

// Destructor.  This will be called after the solution finishes, to free up memory allocated by the class.
OpticalCoupling_Recombination::~OpticalCoupling_Recombination () { 
	cout << "OpticalCoupling_Recombination destructor";

}

void OpticalCoupling_Recombination::DefineDependencies (std::vector<des_data::des_id>& dependencies)
{ 
     cout << "DefineDependencies" << endl;
     dependencies.push_back (des_data::des_id (des_data::scalar,des_data::vertex,"eDensity"));
     dependencies.push_back (des_data::des_id (des_data::scalar,des_data::vertex,"hDensity"));
}

/*-----------------------------------------------------------------
 * DefineJacobians
 *--------------------------------------------------------------------*/
// Create Jacobian matrices and allocate memory for all non-zero elements.  These will
// correspond to the non-zero elements in LCMatrix.  We only create 2 matrices, as J_elec_n = J_hole_n
// and J_elec_p = J_hole_p. They are static class members, which makes the Jacobians accessible to all
// instances of class OpticalCoupling_Recombination.

void OpticalCoupling_Recombination::DefineJacobians (des_id_to_jacobian_map& J_elec, des_id_to_jacobian_map& J_hole)
{

}

// find the index of the mesh region corresponding to region_name.
int find_layer(const des_mesh *mesh, std::string region_name) 
{
      int i=0;
      while (i< mesh->size_region()) {
               if (region_name == mesh->region(i)->name()) {
			return i;
               }
               i++;
      }
      cout << "region " << region_name << " not found" << endl; 
      assert( 1==0 && "find_layer: specified layer not found");
}


/*-----------------------------------------------------------------
 * Function: Compute_parallel
 *--------------------------------------------------------------------*/
// This function is called at least once for each region in the mesh.  If running with multiple
// threads it may be called multiple times for subsets of each region.  The function
// must only update the vertices (and Jacobian rows) identified in input.


void OpticalCoupling_Recombination::Compute_parallel (const PMI_NonLocal_Recombination::Input& input,
				PMI_NonLocal_Recombination::Output& output)
{
	const des_mesh* mesh = Mesh ();
	des_data* data = Data (); 

	// This would be better done in the constructor, but this is the first time we know which region
 	// this instance is working on.
        if (layer_id == -1) {
                layer_id = find_layer(mesh, input.region->name());
                assert (layer_id >= 0 && "assertion failed - layer_id");
		Brad = Brad_list.at(layer_id);
		//cout << layer_id << ", Set Brad " << Brad << endl;
        }

	// do matrix multiplication
	// update Jacobians

	const double* eDensity = data->ReadScalar(des_data::vertex, "eDensity");
	const double* hDensity = data->ReadScalar(des_data::vertex, "hDensity");
	const double* nie = data->ReadScalar(des_data::vertex, "EffectiveIntrinsicDensity");

	for (size_t it = 0; it < input.vertices.size(); it++)
	{
		size_t id = input.vertices.at(it);

		// zero these first...
		output.elec.at(id) = 0.0;
		output.hole.at(id) = 0.0;

		des_jacobian_iterator it_begin = LCMatrix->lower_bound(id, 0);
		des_jacobian_iterator it_end   = LCMatrix->lower_bound(id+1, 0);

		// iterate over all non-zero elements in row [id] of LCMatrix
		for( des_jacobian_iterator it2 = it_begin; it2!=it_end; it2++) {
			size_t col = it2.col();
			double newterm = *(LCMatrix->element(id,col))*( hDensity[col]*eDensity[col] - nie[col]*nie[col] );
		
			// generation at B due to recombination at A. Note Brad and vertex measure have already been applied in
			// the LCMatrix element.
			output.elec.at(id) -= newterm;
			output.hole.at(id) -= newterm;

		}



	}
}

/*-----------------------------------------------------------------
 * Function: new_PMI_NonLocal_Recombination
 * Parameters: env - not documented by Sentaurus.
 * Description: Function which creates an object of the
 *              class OpticalCoupling_Recombination.
 *--------------------------------------------------------------------*/
extern "C" PMI_NonLocal_Recombination* new_PMI_NonLocal_Recombination(const PMI_Device_Environment& env)
{
    return new OpticalCoupling_Recombination (env);
}



/* --------------------------------------------------------------
 * These are used to make a vector of vertex indices sorted by y-coodinate
 *--------------------------------------------------------------*/

typedef struct {
 size_t index;
 float y;
 double measure; } vert_y;

bool mycompare (const vert_y &lhs, const vert_y &rhs){
  return lhs.y < rhs.y ;
}

/*---------------------------------------------------------------------
 * Function: Read_LCMatrix
 *---------------------------------------------------------------------*/
void OpticalCoupling_Recombination::Read_LCMatrix()
{
  const des_mesh* mesh = Mesh ();
  des_data* data = Data ();

  std::string line;
  std::vector<std::string> items;
  std::string item;

  int li, lk;
  int elem_count = 0;

  // Build a vector of vertices sorted by y-coord for each region
  vector< vector<vert_y> > yvals;
  for (int i = 0; i< mesh->size_region(); i++) {
         vector< vert_y > ys;
         for ( int k =0; k < mesh->region(i)->size_vertex(); k++) {
////         	cout << "** vertex is" << mesh->region(i)->vertex(k) << "for k=" << i <<endl;
                des_vertex *vert = mesh->region(i)->vertex(k);
		const double *c = vert->coord();
                vert_y* vert_struct = new vert_y;
                vert_struct->index = vert->index() ;
                vert_struct->y = c[1];
		// TODO: Add measures!

		ys.push_back( *vert_struct );
          }
          sort(ys.begin(),ys.end(),mycompare);
          yvals.push_back(ys);
  }

  ifstream inFile(LCMatrix_file.c_str(), ifstream::in);
  if (inFile.is_open() )
  {

     // Read a line to a vector of items
     while( getline(inFile, line) )
     { 
         if (line[line.size() -1] == '\r')
            line.erase(line.size()-1);
	    stringstream linestream(line);
            //cout << linestream.str() << endl;
            while( getline(linestream,item,',')) {
               items.push_back(item);
         }

       	 //  Collect Brad values
	 if (items[0] == "Brad") {
            string region_name = items[1];
            int i=find_layer(mesh, region_name);
            sscanf(items[3].c_str(), "%lf", &Brad_list[i]);
            cout << region_name << ' ' << Brad_list[i] << endl;
         }

         //  New layer pair
         else if (items[0] == "**") {
              li = find_layer(mesh, items[2]);
              lk = find_layer(mesh, items[4]);
              cout << "Elements so far: " << elem_count << endl;
              cout << "Coupling between regions " << mesh->region(li)->name() << " " << li << " " << mesh->region(lk)->name() << " " << lk << endl;
         }
         else if ( items.size() == 3) {
	      // we have a line defining a matrix element LC(yA, yB) = P.
              float yA, yB;
              double v1, v2, P;
              sscanf( items[0].c_str(), "%lf", &v1);
              sscanf( items[1].c_str(), "%lf", &v2);
              sscanf( items[2].c_str(), "%lf", &P);
              yA = v1; yB = v2;
	      // yA and yB must be read from the file as doubles, but compared with sdevice as single floats.
	      double Brad = Brad_list[li];

	      // Throw out coupling less than threshold.
	      if (P > 8e-3) {
                // declare a reference struct for comparisons
  	        vert_y ref;

	        // find first and last vertices with y=yB in layer lk.
                ref.y = yB;
                vector<vert_y>::iterator  lowB = lower_bound(yvals[lk].begin(), yvals[lk].end(), ref, mycompare);
                vector<vert_y>::iterator highB = upper_bound(lowB, yvals[lk].end(), ref, mycompare);

	        // find first and last vertices with y=yA in layer li.
	        ref.y = yA;
                vector<vert_y>::iterator  lowA = lower_bound(yvals[li].begin(), yvals[li].end(), ref, mycompare);
                vector<vert_y>::iterator highA = upper_bound(lowA, yvals[li].end(), ref, mycompare);

		const double*const* measure = data->ReadMeasure ();
	        for (vector<vert_y>::iterator itA = lowA; itA < highA; itA++ ) {
			// Calculate the measure for vertex A.
			double measureA = 0.0;
			des_vertex* vertexA = mesh->vertex (itA->index);

			// iterate over elements associated with vertex A
			for (size_t eli = 0; eli < vertexA->size_element (); eli++) {
				des_element* el = vertexA->element (eli);
				if (el->bulk()->name() == mesh->region(li)->name() ) {
					for (size_t vi = 0; vi < el->size_vertex (); vi++) {
					   des_vertex* v = el->vertex (vi);
					   if (v == vertexA) {
					       measureA += measure [el->index()][vi];
					   }
					}
				}
			}

                      int inner_count =0;
  		      for (vector<vert_y>::iterator itB = lowB; itB < highB; itB++ ) {
                            LCMatrix->define_element(itB->index, itA->index);
                            *(LCMatrix->element(itB->index, itA->index)) = P*Brad*measureA;
			    assert (P > 0.0 && "P <=0!");
			    assert (Brad > 0.0 && "Brad <=0!");
			    assert (measureA > 0.0 && "measureA <=0!");
		//	    cout << "P " << P << " Brad " << Brad << " measureA " << measureA << endl;
                            elem_count++;
			    inner_count++;
                      }
                   //  cout <<  "yB=" << yB << " inner count " << inner_count << endl;
	        }

             }

         }
        items.clear();
  
     }
  inFile.close();
  }
  yvals.clear();
  cout << "Read_LCMatrix() finished. Element count "<< elem_count << endl;
  cout << "Size of LCMatrix " << sizeof(*LCMatrix);
}
