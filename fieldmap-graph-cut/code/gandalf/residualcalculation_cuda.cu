#include "mex.h"
#include "matrix.h"
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <cuda_runtime.h>
#include <string.h>

__global__
void calcresKernel(double *resi, double *P_mat_r, double *P_mat_i, double *i_mage_r, double *i_mage_i, int xsize, int ysize, int nr_fms, int nr_kr, int nr_echo)
{
    //extern __shared__ double s_data[];

    int ix = threadIdx.x + blockDim.x * blockIdx.x;
    int iy = threadIdx.y + blockDim.y * blockIdx.y;
    int ifms = threadIdx.z + blockDim.z * blockIdx.z;

    int ikr,iecho,k,l;
    double sum, imagsum, realsum, a1, a2, b1, b2;
    double minimum;    
    
    /*
    for (int ism1 = 0; ism1 < nr_echo; ism1++)
    {
        s_data[ism1*2] = i_mage_r[ix+iy*xsize+ism1*xsize*ysize];
        s_data[ism1*2 + 1] = i_mage_i[ix+iy*xsize+ism1*xsize*ysize];
    }
    __syncthreads;
    */
    
    //Initialize with maximum double 
    minimum = DBL_MAX;
    
    if (ix < xsize && iy < ysize && ifms < nr_fms)  
    {
        for (ikr=0; ikr<nr_kr; ikr++){
            sum = 0;

            // sum is calculated over k and l
            for (k=0;k<nr_echo;k++){
                realsum = 0;
                imagsum = 0;

                for (l=0;l<nr_echo;l++){
    //                          sum = sum +P_mat[k][l][ifms][ikr]*image[ix][iy][l]; but complex multiplication
                    a1 = P_mat_r[k+l*nr_echo+ifms*nr_echo*nr_echo+ikr*nr_echo*nr_echo*nr_fms];
                    b1 = P_mat_i[k+l*nr_echo+ifms*nr_echo*nr_echo+ikr*nr_echo*nr_echo*nr_fms];

                    //a2 = s_data[2*l];
                    //b2 = s_data[2*l + 1];
                    a2 = i_mage_r[ix+iy*xsize+l*xsize*ysize];
                    b2 = i_mage_i[ix+iy*xsize+l*xsize*ysize];
                    realsum = realsum + a1*a2 - b1*b2;
                    imagsum = imagsum + a1*b2 + a2*b1;              
                }
                // calculation of square of absolute value of complex product:
                sum = sum + realsum*realsum + imagsum*imagsum;
            }
            if (sum<minimum){
                minimum = sum;
    //                         printf("sum is smaller!\n");  
            }
        }
        //printf("assigning value to residual\n");
        resi[ifms + ix*nr_fms + iy*nr_fms*xsize] = minimum;
    }
}

static void HandleError( cudaError_t err,
                const char *file,
                int line )
{
if (err != cudaSuccess) {
    printf( "%s in %s at line %d\n", 
            cudaGetErrorString( err ),
            file, line );
    // exit( EXIT_FAILURE );
}
}

#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

inline dim3 computeGrid1D(const dim3 &block, const int w)
{
    int num = (w + block.x -1) / (block.x);
    //int num = w + block.x - 1;
    return dim3(num, 1, 1);   // TODO (3.2) compute 1D grid size from block size
}

inline dim3 computeGrid2D(const dim3 &block, const int w, const int h)
{
    int num1 = (w + block.x -1) / (block.x);
    int num2 = (h + block.y -1) / (block.y);
    //int num = 1;
    return dim3(num1, num2, 1);   // TODO (3.2) compute 2D grid size from block size
}

inline dim3 computeGrid3D(const dim3 &block, const int w, const int h, const int s)
{
    int num1 = (w + block.x -1) / (block.x);
    int num2 = (h + block.y -1) / (block.y);
    int num3 = (s + block.z -1) / (block.z);
    //int num = 1;
    return dim3(num1, num2, num3);   // TODO (3.2) compute 3D grid size from block size
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
//example out- and input matrix sizes: [921,400,400] = calcResidualInC(Phi[20,20,921,26],Image[400,400,20])
{
    //====================  Variable definition  ======================
    
    int number_dims_out;
    const mwSize *dim_phi, *dim_array_image;
    mwSize dims_out[] = {1,2,3};
    
    //Inputs
    int num_fms;
    int num_kr;
    int num_echo;
    int xres, yres;
    double *P_matrix_r;
    double *P_matrix_i;
    double *in_image_r;
    double *in_image_i;
    
    //Output
    double *residual; //mxArray could also be a data type
    
    //CUDA arrays
    double *d_P_matrix_r = NULL;
    double *d_P_matrix_i = NULL;
    double *d_in_image_r = NULL;
    double *d_in_image_i = NULL;
    double *d_residual = NULL;
    
    int fmscnt;
    int ix,iy,ifms,nx;
    
    //==============================================================
    
    //for output initialization, 3 dimensions are created
    number_dims_out = 3;
    //size of Psi matrix is extracted: 20x20x921x26
    dim_phi = mxGetDimensions(prhs[0]);
    num_fms = dim_phi[2];
    num_kr = dim_phi[3];
    //size of input image is extracted: 400x400x20 xres*yres*nrecho
    dim_array_image = mxGetDimensions(prhs[2]);
    num_echo = dim_array_image[2];
    xres = dim_array_image[0];
    yres = dim_array_image[1];  
    
    // ============= Data Initialization and Allocation ==============
    //Output
    dims_out[0] = num_fms;
    dims_out[1] = xres;
    dims_out[2] = yres;
    plhs[0] = mxCreateNumericArray(number_dims_out, dims_out, mxDOUBLE_CLASS, mxREAL);
    residual = mxGetPr(plhs[0]);
    
    //Get input data
    P_matrix_r = mxGetPr(prhs[0]);
    P_matrix_i = mxGetPr(prhs[1]);
    in_image_r = mxGetPr(prhs[2]);
    in_image_i = mxGetPr(prhs[3]);
    
    // Allocate cuda memory for CUDA arrays
    HANDLE_ERROR(cudaMalloc( &d_P_matrix_r, dim_phi[0]*dim_phi[1]*dim_phi[2]*dim_phi[3]*sizeof(double))); //CUDA_CHECK;
    HANDLE_ERROR(cudaMalloc( &d_P_matrix_i, dim_phi[0]*dim_phi[1]*dim_phi[2]*dim_phi[3]*sizeof(double))); //CUDA_CHECK;
    HANDLE_ERROR(cudaMalloc( &d_in_image_r, num_echo*xres*yres*sizeof(double))); //CUDA_CHECK;
    HANDLE_ERROR(cudaMalloc( &d_in_image_i, num_echo*xres*yres*sizeof(double))); //CUDA_CHECK;
    HANDLE_ERROR(cudaMalloc( &d_residual, num_fms*xres*yres*sizeof(double))); //CUDA_CHECK;
    
    //Copy memory
    HANDLE_ERROR(cudaMemcpy(d_P_matrix_r, P_matrix_r, dim_phi[0]*dim_phi[1]*dim_phi[2]*dim_phi[3]*sizeof(double), cudaMemcpyHostToDevice)); //CUDA_CHECK; 
    HANDLE_ERROR(cudaMemcpy(d_P_matrix_i, P_matrix_i, dim_phi[0]*dim_phi[1]*dim_phi[2]*dim_phi[3]*sizeof(double), cudaMemcpyHostToDevice)); //CUDA_CHECK;
    HANDLE_ERROR(cudaMemcpy(d_in_image_r, in_image_r, num_echo*xres*yres*sizeof(double), cudaMemcpyHostToDevice)); //CUDA_CHECK;
    HANDLE_ERROR(cudaMemcpy(d_in_image_i, in_image_i, num_echo*xres*yres*sizeof(double), cudaMemcpyHostToDevice)); //CUDA_CHECK;
//     cudaMemcpy(d_residual, residual, num_fms*xres*yres*sizeof(double), cudaMemcpyHostToDevice)); //CUDA_CHECK;
    
    //================ Calculation of residual function ============
    //calcres(residual,P_matrix_r,P_matrix_i,in_image_r,in_image_i,xres,yres,num_fms,num_kr,num_echo);
    
    // calculate block and grid size
    dim3 block = dim3(8, 8, 16);     // Specify suitable block size
    dim3 grid = computeGrid3D(block, xres, yres, num_fms);
    size_t smBytes = 2 * num_echo * sizeof(double) ;
    //printf("%d, %d, %d", grid.x, grid.y, grid.z);

    // run CUDA kernel
    calcresKernel <<<grid,block>>> (d_residual, d_P_matrix_r, d_P_matrix_i, d_in_image_r, d_in_image_i,xres,yres,num_fms,num_kr,num_echo); //CUDA_CHECK;
    HANDLE_ERROR(cudaMemcpy(residual, d_residual, num_fms*xres*yres*sizeof(double), cudaMemcpyDeviceToHost)); //CUDA_CHECK;
    
    #if MX_HAS_INTERLEAVED_COMPLEX
        //printf("it has interleave!");
    #endif 
    
    HANDLE_ERROR(cudaFree(d_P_matrix_r));
    HANDLE_ERROR(cudaFree(d_P_matrix_i));
    HANDLE_ERROR(cudaFree(d_in_image_r));
    HANDLE_ERROR(cudaFree(d_in_image_i));
    HANDLE_ERROR(cudaFree(d_residual));
    
//     delete[] imgIn;
//     delete[] imgOut;    
        
//     do not call mxDestroyArray or mxFree on an mxArray returned in a left-side argument of a mex-file
//     mxDestroyDouble(in_image);
//     mxDestroyDouble(P_matrix);
        
}