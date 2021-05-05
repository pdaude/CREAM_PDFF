#include <optnet/optnet.hxx>
#include <optnet/_base/array_ref.hxx>
#include <exception>
#include <iostream>
#include <stdlib.h>

#include <fstream>                         // file I/O support
#include <cstdlib>                         // support for exit()
#include <cmath>
#include <vector>

#include <string>
#include <sstream>
#include <time.h>

#include "mex.h"
#include "string.h"
# define PI 3.14159265
using namespace std;
bool detect_one_surface(int*			elemap,
						int*			node_costs,
						int				x_size,
						int				y_size,
						int				z_size,
						int				x_smooth,
						int				y_smooth,
						bool			x_circular,
						bool			y_circular);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Read in parameters from header file.*/
    int ifield, nfields;
    mwIndex    jstruct;
    mwSize     NStructElems;
    const char *fnames;  /* pointers to field names */

    size_t buflen;
    mwSize sizebuf;
    int status;

    mxArray *tmp;
    string Dataname, fname;
    double DELTA, stepsize, nodecost_min, nodecost_max;
    int freqRange, Numlayer, fatpeaksNo, z, Date, smooth_constraint;
    //* get input arguments */
    nfields = mxGetNumberOfFields(prhs[0]);
    NStructElems = mxGetNumberOfElements(prhs[0]);

    for (ifield=0; ifield< nfields; ifield++){
    fnames = mxGetFieldNameByNumber(prhs[0],ifield);
    tmp = mxGetFieldByNumber(prhs[0],0,ifield);
    stringstream ss;
        ss << fnames;
        ss >> fname;

    if (fname == "Dataname"){
         buflen = mxGetN(tmp)*sizeof(mxChar)+1;
         char* buf = new char[buflen];
         status = mxGetString(tmp, buf, (mwSize)buflen);
         stringstream ss2;
         ss2 << buf;
         ss2 >> Dataname;
         cout << "Dataset: " << Dataname << endl;
         delete [] buf;
        }

    if (fname == "Delta"){
         DELTA= mxGetScalar(tmp);
//          double* t;
//          t = (double*)mxGetData(tmp);
         cout << "The TE period is " << DELTA << endl;
    }

    if (fname == "freqRange"){
         freqRange = mxGetScalar(tmp);
         cout << "The fieldmap search range is " << 2*freqRange << endl;
    }

    if (fname == "fatpeaksNo"){
         fatpeaksNo = mxGetScalar(tmp);
         cout << "Number of fat peaks used in the model is " << fatpeaksNo << endl;
    }

    if (fname == "stepsize"){
         stepsize = mxGetScalar(tmp);
         cout << "The stepsize for frequency change between layers is " << stepsize << endl;
    }

    if (fname == "NumLayer"){
         Numlayer = mxGetScalar(tmp);
         cout << "The number of layers (z-direction) in the graph is " << Numlayer << endl;
    }

    if (fname == "Z"){
         z = mxGetScalar(tmp);
         cout << "This is No." << z << " slice in 3rd/z dimension of original data." << endl;
    }

    if (fname == "Date"){
         Date = mxGetScalar(tmp);
         cout << "Date:" << Date << endl;
    }

    if (fname == "constraint"){
         smooth_constraint = mxGetScalar(tmp);
         cout << "The smooth constraint for graph is:" << smooth_constraint << endl;
    }

    if (fname == "residueMin"){
         nodecost_min = mxGetScalar(tmp);
         cout << "The minimum nodecost in graph is:" << nodecost_min << endl;
    }

    if (fname == "residueMax"){
         nodecost_max = mxGetScalar(tmp);
         cout << "The maximum nodecost in graph is:" << nodecost_max << endl;
    }
    }

   /* 2. Read in nodecost for graph. (Residue in Matlab) */
   mwSize elements, numOfDim, count;
   int HEIGHT, WIDTH;
   const int *dimArray;
   numOfDim = mxGetNumberOfDimensions(prhs[1]);
   dimArray = mxGetDimensions(prhs[1]);
   HEIGHT = dimArray[0];  //x -- 1st dimension in Matlab
   WIDTH = dimArray[1];//y -- 2nd dimension in Matlab
   elements = mxGetNumberOfElements(prhs[1]);
   cout << "The number of elements read in: " << elements << endl;
   double *costfunc ;
   costfunc = (double*)mxGetPr(prhs[1]); // costfunc now is a const pointer, cannot be deleted;
   // rescale cost to 0~1000;
   int *NodeCost = new int[HEIGHT*WIDTH*Numlayer];
   for (int z = 0; z < Numlayer; z++)
   {
    for (int y = 0; y < WIDTH; y++)
    {
     for (int x = 0; x < HEIGHT; x++)
     {
         NodeCost[z*HEIGHT*WIDTH + y*HEIGHT + x] = int((costfunc[z*HEIGHT*WIDTH + y*HEIGHT + x] - nodecost_min)*((1000)-0)/(nodecost_max - nodecost_min));
     } // for x--WIDTH
    } // for y--HEIGHT
   } // for z--layers
   
   cout << "First node cost = "<< costfunc[0] << endl;
   cout << "Last node cost = "<< costfunc[Numlayer*HEIGHT*WIDTH - 1] << endl;

   /* 3. Define fieldmap as the output, run graph algorithm. */
   double* surface;
   plhs[0] = mxCreateDoubleMatrix(HEIGHT, WIDTH, mxREAL);
   surface = mxGetPr(plhs[0]);
   int* tempmap = new int[HEIGHT*WIDTH];
   cout << "FINDING OPTIMAL SURFACE..." << endl;
   detect_one_surface(tempmap, NodeCost, HEIGHT, WIDTH, Numlayer, smooth_constraint, smooth_constraint, 0, 0);
   cout << "GRAPH SEARCH FINISHED..." << endl;
   for (int y = 0; y < WIDTH; y++)
    {
     for (int x = 0; x < HEIGHT; x++)
      {
//           fieldmap[y*HEIGHT+x] = -freqRange + tempmap[y*HEIGHT+x]*stepsize;
         surface[y*HEIGHT+x] = tempmap[y*HEIGHT+x]+1.0;
      }
    }

   delete [] NodeCost;
   delete [] tempmap;
}


bool detect_one_surface(int*			elemap,
						int*			node_costs,
						int				x_size,
						int				y_size,
						int				z_size,
						int				x_smooth,
						int				y_smooth,
						bool			x_circular,
						bool			y_circular)

{
	typedef optnet::optnet_ia_3d<int, long int, optnet::net_f_xy> SolverType;
    //typedef optnet::optnet_pr_3d<float, double, optnet::net_f_xy> SolverType;

	try
	{
		// Create a reference object to the input node costs array.
		// The size of this array is x_size * y_size * z_size.
		SolverType::cost_array_ref_type node_costs_ref(node_costs,
													   x_size,
													   y_size,
													   z_size);

		// Create a reference object to the output result.
		// The size of the result is x_size * y_size
		SolverType::net_ref_type elemap_ref(elemap,
											x_size,
											y_size,
											1);

		// Declare a solver object.
		SolverType solver;

		// Create and initialize the solver with the node costs.
		solver.create(node_costs_ref);

		// Set smoothness and circularity parameters
		solver.set_params(x_smooth, y_smooth, x_circular, y_circular);

		// Solve the problem.
		solver.solve(elemap_ref,
					 NULL); // The output maximum flow parameter has no use
							// to us, simply set it to NULL.

	}
	catch (std::exception& e)
	{
		// Handle any errors here
		std::cout << e.what() << std::endl;

		return false;
	}

	return true;
}

