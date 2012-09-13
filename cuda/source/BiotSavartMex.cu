#include "mex.h"
#include "matrix.h"
#include "cutil_math.h"
#include <cuda.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <omp.h>
//#include "cuPrintf.cu"

// Input Arguments
#define F1_in       prhs[0]
#define F2_in       prhs[1]
#define P_in        prhs[2]
#define GAMMA_in    prhs[3]
#define RC_in       prhs[4]
#define D_in        prhs[5]
#define CMOD_in     prhs[6]
#define CO_in       prhs[7]
#define TYPE_in     prhs[8]
#define GPU_in      prhs[9]
#define GPUHW_in      prhs[10]

// Output Arguments
#define UIND_out    plhs[0]
#define L_out       plhs[1]


// Precision
#define precis      double

// CUDA stuff
#define pi          3.14159265358979324
#define tpb         128

typedef struct
        {
                double x;
                double y;
                double z;
        } vector3;

// CUDA get errors
void checkCUDAError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err)
    {
        mexPrintf("Cuda error: %s: %s.\n", msg,
                                  cudaGetErrorString( err) );
    }                        
}

__constant__ int dev_np;
__constant__ int dev_nf;
__constant__ double cutoff;
__constant__ double cored;
__constant__ int threadsperblock;

// CUDA kernel for Biot-Savart Particle Interaction (Direct Sum Version)
__device__ void biot_p2p(double2 bfx, double2 bfy, double2 bfz, double2 brg, double2 bpxy, double1 bpz, double* uindx, double* uindy, double* uindz)
{

    double r1, r2, len, ldx, ldy, ldz, pxx1, pxx2, pyy1, pyy2, pzz1, pzz2, r1dr2,ubar, r1tr2, den;
    
    ldx=bfx.x-bfx.y;
    ldy=bfy.x-bfy.y;
    ldz=bfz.x-bfz.y;

    pxx1=bpxy.x-bfx.x;
    pxx2=bpxy.x-bfx.y;
    pyy1=bpxy.y-bfy.x;
    pyy2=bpxy.y-bfy.y;
    pzz1=bpz.x-bfz.x;
    pzz2=bpz.x-bfz.y;

    r1 = sqrt(pxx1*pxx1+pyy1*pyy1+pzz1*pzz1);
    r2 = sqrt(pxx2*pxx2+pyy2*pyy2+pzz2*pzz2);

    r1dr2=pxx1*pxx2+pyy1*pyy2+pzz1*pzz2;
    r1tr2=r1*r2;

    len=ldx*ldx + ldy*ldy + ldz*ldz; //L^2
    den=r1tr2*(r1tr2 + r1dr2) + cored*len; //cored already squared
    ubar=(brg.x*(r1+r2))*(1.000/(4.000*pi));
    ubar = ubar/den;
 
    uindx[0] += ubar*(pyy1*pzz2-pzz1*pyy2);
    uindy[0] += ubar*(pzz1*pxx2-pxx1*pzz2);
    uindz[0] += ubar*(pxx1*pyy2-pyy1*pxx2);
}


// CUDA kernel for Biot-Savart Particle Interaction (Shared Version)
__device__ vector3 biot_p2p_sh(double2 bfx, double2 bfy, double2 bfz, double2 brg, double2 bpxy, double1 bpz, vector3 uind)
{

    double r1, r2, len, ldx, ldy, ldz, pxx1, pxx2, pyy1, pyy2, pzz1, pzz2, r1dr2,ubar, r1tr2, den;
    
    ldx=bfx.x-bfx.y;
    ldy=bfy.x-bfy.y;
    ldz=bfz.x-bfz.y;

    pxx1=bpxy.x-bfx.x;
    pxx2=bpxy.x-bfx.y;
    pyy1=bpxy.y-bfy.x;
    pyy2=bpxy.y-bfy.y;
    pzz1=bpz.x-bfz.x;
    pzz2=bpz.x-bfz.y;

    r1 = sqrt(pxx1*pxx1+pyy1*pyy1+pzz1*pzz1);
    r2 = sqrt(pxx2*pxx2+pyy2*pyy2+pzz2*pzz2);

    r1dr2=pxx1*pxx2+pyy1*pyy2+pzz1*pzz2;
    r1tr2=r1*r2;

    len=ldx*ldx + ldy*ldy + ldz*ldz; //L^2
    den=r1tr2*(r1tr2 + r1dr2) + cored*len; //cored already squared
    ubar=(brg.x*(r1+r2))*(1.000/(4.000*pi));
    ubar = ubar/den;
 
    uind.x += ubar*(pyy1*pzz2-pzz1*pyy2);
    uind.y += ubar*(pzz1*pxx2-pxx1*pzz2);
    uind.z += ubar*(pxx1*pyy2-pyy1*pxx2);

    return uind;
}

__global__ void BiotSavart_naive(double2 *p, double1 *pz, double2 *fx, double2 *fy, double2 *fz, double2 *rg, double2 *uind, double1 *uindz)
{
//Get thread's global index
int k;
int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if(idx<dev_np){
    //Loop over all source filaments
    for(k=0; k < dev_nf; k++)
    {
      biot_p2p(fx[k], fy[k], fz[k], rg[k], p[idx], pz[idx], &uind[idx].x, &uind[idx].y, &uindz[idx].x);
    }
  }

}

__global__ void BiotSavart_sh(double2 *p, double1 *pz, double2 *fx, double2 *fy, double2 *fz, double2 *rg, double2 *uind, double1 *uindz)
{
  __shared__ double2 shFx[tpb];
  __shared__ double2 shFy[tpb];
  __shared__ double2 shFz[tpb];
  __shared__ double2 shRg[tpb];

//Get thread's global index
int k,i,j,tile;
int gid = blockIdx.x * blockDim.x + threadIdx.x;
int tilesize = blockDim.x;
int tileloc = threadIdx.x;

//Local Vars
vector3 myUind;
int fid, tileid;


    //Loop over all tiles of source filaments (tile size = block size)
    for(k=0, tile=0; k < dev_nf; k+= tilesize, tile++)
    {
      tileid = tile*tilesize+tileloc;
      fid = tile*tilesize+k;
         
      shFx[threadIdx.x] = fx[tileid];
      shFy[threadIdx.x] = fy[tileid];
      shFz[threadIdx.x] = fz[tileid];
      shRg[threadIdx.x] = rg[tileid];
     

      // Syncronize threads before using shared mem
      __syncthreads();
 
      if(gid<dev_np){
        for(i=0; i<tilesize; i++)
        {
        fid = tile*tilesize+i;
        if(fid<dev_nf)
        myUind = biot_p2p_sh(shFx[i], shFy[i], shFz[i], shRg[i], p[gid], pz[gid], myUind);
        }
      }
      __syncthreads(); 
    }
  uind[gid].x = myUind.x;
  uind[gid].y = myUind.y;
  uindz[gid].x = myUind.z;
}


// CUDA kernel that calculates segment length
__global__ void calcLengthOnly(double2 *a, double2 *b, double2 *c, double1 *d)
{
    //Get thread's global index
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if(idx<dev_np){
       d[idx].x = sqrt((a[idx].x-a[idx].y)*(a[idx].x-a[idx].y) + (b[idx].x-b[idx].y)*(b[idx].x-b[idx].y) + (c[idx].x-c[idx].y)*(c[idx].x-c[idx].y));      
    }

}


// CUDA Device properties and problem size calc
static int cuda_setup(int probsize, int numthreads, int numcards, int devicenum)
{
    cudaDeviceProp prop;
    double maxglobal;
    double maxshared;
    int blocksize;

    cudaGetDeviceProperties(&prop, devicenum);
    //mexPrintf("Device Name: %s \n",prop.name);
    //mexPrintf("Global Memory: %d \n",prop.totalGlobalMem);
    //mexPrintf("Shared Memory Per Block: %d \n",prop.sharedMemPerBlock);
    //mexPrintf("Number of MPs: %d \n",prop.multiProcessorCount);   
 
    maxglobal = prop.totalGlobalMem/sizeof(double);
    //mexPrintf("Max number of elements in Global: %3.9f \n",maxglobal);
    maxshared = prop.sharedMemPerBlock/sizeof(double);
    //mexPrintf("Max number of elements in Shared Per Block: %3.9f \n",maxshared);

    blocksize = (probsize+(numthreads-1))/numthreads;
    //mexPrintf("Given thread count: %d \n",numthreads);
    //mexPrintf("Computed Block Size: %d \n",blocksize);  
    return blocksize;
}



void biot_p2p_nogpu(double2 bfx, double2 bfy, double2 bfz, double2 brg, double2 bpxy, double1 bpz, double* uindx, double* uindy, double* uindz, double cd, double co)
{
    //mexPrintf("Cored: %3.9f \n",cd);
    double r1, r2, len, ldx, ldy, ldz, pxx1, pxx2, pyy1, pyy2, pzz1, pzz2, r1dr2,ubar, r1tr2, den;
    
    ldx=bfx.x-bfx.y;
    ldy=bfy.x-bfy.y;
    ldz=bfz.x-bfz.y;

    pxx1=bpxy.x-bfx.x;
    pxx2=bpxy.x-bfx.y;
    pyy1=bpxy.y-bfy.x;
    pyy2=bpxy.y-bfy.y;
    pzz1=bpz.x-bfz.x;
    pzz2=bpz.x-bfz.y;

    r1 = sqrt(pxx1*pxx1+pyy1*pyy1+pzz1*pzz1);
    r2 = sqrt(pxx2*pxx2+pyy2*pyy2+pzz2*pzz2);

    r1dr2=pxx1*pxx2+pyy1*pyy2+pzz1*pzz2;
    r1tr2=r1*r2;

    len=ldx*ldx + ldy*ldy + ldz*ldz; //L^2
    den=r1tr2*(r1tr2 + r1dr2) + cd*len; //cored already squared
    ubar=(brg.x*(r1+r2))*(1.000/(4.000*pi));
    ubar = ubar/den;
 
    uindx[0] += ubar*(pyy1*pzz2-pzz1*pyy2);
    uindy[0] += ubar*(pzz1*pxx2-pxx1*pzz2);
    uindz[0] += ubar*(pxx1*pyy2-pyy1*pxx2);
}

void BiotSavart_nogpu(double2 *p, double1 *pz, double2 *fx, double2 *fy, double2 *fz, double2 *rg, double2 *uind, double1 *uindz, int n, int f_n, double cd, double co)
{
omp_set_num_threads(4);
int j,k,omptnum;
  // Loop over all particles  
  for(j=0; j < n; j++)
  {
    //Loop over all source filaments
    for(k=0; k < f_n; k++)
    {
      biot_p2p_nogpu(fx[k], fy[k], fz[k], rg[k], p[j], pz[j], &uind[j].x, &uind[j].y, &uindz[j].x, cd, co);
    }
  }
  
}



// Resolve arrays to working 3D matrix
static void crunch_array(double *arrayin, double3 *arrayout, mwSize *di)
{
   
    mwSize ns = di[0];
    mwSize nd = di[1];
    mwSize nt = di[2];
    mwSize nb = di[3];

    int i = 0;
    int m,n,p;
    
    

        for(p=0; p<nb;p++){ 
          for(n=0; n<nt; n++){ 
             for(m=0; m<ns; m++){
                 
            arrayout[i].x = arrayin[m+0*ns+n*ns*nd+p*ns*nt*nd];
            arrayout[i].y = arrayin[m+1*ns+n*ns*nd+p*ns*nt*nd];
            arrayout[i].z = arrayin[m+2*ns+n*ns*nd+p*ns*nt*nd];
          
              i++;
            }
        }
    }
}


// Two input arrays gamma and rc stuffed into double2
static void crunch_rg(double *arrayin1, double *arrayin2, double2 *arrayout, mwSize *di)
{
    mwSize ns = di[0];
    mwSize nd = 1;
    mwSize nt = di[2];
    mwSize nb = di[3];
    int i = 0;
    int m,n,p;
        for(p=0; p<nb;p++){ 
          for(n=0; n<nt; n++){ 
            for(m=0; m<ns; m++){     
            arrayout[i].x = arrayin1[m+0*ns+n*ns*nd+p*ns*nt*nd];
            arrayout[i].y = arrayin2[m+0*ns+n*ns*nd+p*ns*nt*nd];
            i++;
            }
        }
    }
}

// Two input arrays F1 and F2 stuffed into 3x double2
static void crunch_f(double *arrayin1, double *arrayin2, double2 *arrayout1,double2 *arrayout2,double2 *arrayout3, mwSize *di)
{
    mwSize ns = di[0];
    mwSize nd = di[1];
    mwSize nt = di[2];
    mwSize nb = di[3];
    int i = 0;
    int m,n,p;
        for(p=0; p<nb;p++){ 
          for(n=0; n<nt; n++){ 
            for(m=0; m<ns; m++){     
            arrayout1[i].x = arrayin1[m+0*ns+n*ns*nd+p*ns*nt*nd];
            arrayout1[i].y = arrayin2[m+0*ns+n*ns*nd+p*ns*nt*nd];
            arrayout2[i].x = arrayin1[m+1*ns+n*ns*nd+p*ns*nt*nd];
            arrayout2[i].y = arrayin2[m+1*ns+n*ns*nd+p*ns*nt*nd];
            arrayout3[i].x = arrayin1[m+2*ns+n*ns*nd+p*ns*nt*nd];
            arrayout3[i].y = arrayin2[m+2*ns+n*ns*nd+p*ns*nt*nd];
            i++;
            }
        }
    }
}


// One input array P stuffed into 1x double2 1x double1
static void crunch_p(double *arrayin1, double2 *arrayout1,double1 *arrayout2,mwSize *di)
{
    mwSize ns = di[0];
    mwSize nd = 3;
    mwSize nt = di[2];
    mwSize nb = di[3];
    int i = 0;
    int m,n,p;
        for(p=0; p<nb;p++){ 
          for(n=0; n<nt; n++){ 
            for(m=0; m<ns; m++){     
            arrayout1[i].x = arrayin1[m+0*ns+n*ns*nd+p*ns*nt*nd];
            arrayout1[i].y = arrayin1[m+1*ns+n*ns*nd+p*ns*nt*nd];
            arrayout2[i].x = arrayin1[m+2*ns+n*ns*nd+p*ns*nt*nd];
            i++;
            }
        }
    }
}


// Reassemble uind
static void reassem_uind(double *arrayout, double2 *arrayin1,double1 *arrayin2,mwSize *di)
{
    mwSize ns = di[0];
    mwSize nd = 3;
    mwSize nt = di[2];
    mwSize nb = di[3];
    int i = 0;
    int m,n,p;
        for(p=0; p<nb;p++){ 
          for(n=0; n<nt; n++){ 
            for(m=0; m<ns; m++){     
            arrayout[m+0*ns+n*ns*nd+p*ns*nt*nd] = arrayin1[i].x;
            arrayout[m+1*ns+n*ns*nd+p*ns*nt*nd] = arrayin1[i].y;
            arrayout[m+2*ns+n*ns*nd+p*ns*nt*nd] = arrayin2[i].x;
            i++;
            }
        }
    }
}

// Reassemble double3 into 4d array
static void reassem_array3(double3 arrayin[], double arrayout[], mwSize *di)
{
    mwSize ns = di[0];
    mwSize nd = di[1];
    mwSize nt = di[2];
    mwSize nb = di[3];

    int i = 0;
    int m,n,p;
    
    

        for(p=0; p<nb;p++){ 
          for(n=0; n<nt; n++){ 
             for(m=0; m<ns; m++){
                 
            arrayout[m+0*ns+n*ns*nd+p*ns*nt*nd] = arrayin[i].x;

              i++;
            }
        }
    }
}

// Reassemble double3 into 4d array
static void reassem_array1(double1 arrayin[], double arrayout[], mwSize *di)
{
    mwSize ns = di[0];
    mwSize nd = 1;
    mwSize nt = di[2];
    mwSize nb = di[3];

    int i = 0;
    int m,n,p;

        for(p=0; p<nb;p++){ 
          for(n=0; n<nt; n++){ 
             for(m=0; m<ns; m++){
                 
            arrayout[m+0*ns+n*ns*nd+p*ns*nt*nd] = arrayin[i].x;

              i++;
            }
        }
    }
}

// Biot Savart Launching Routine
static void biot_calc(double *uind ,double *L, double *F1, double *F2,
                   double *P, double *gamma, double *rc, double d, char *cm, 
                   double co, char *casetype, mwSize ndi, mwSize *di, 
                   mwSize *fdi, char *usegpu, double *gpuinfo)
    {
      
    // Dimensions
    mwSize ns = di[0];
    mwSize nd = di[1];
    mwSize nt = di[2];
    mwSize nb = di[3];

    mwSize f_ns = fdi[0];
    mwSize f_nd = fdi[1];
    mwSize f_nt = fdi[2];
    mwSize f_nb = fdi[3];

    // Various problem constants and indices
    mwSize cartprobsize = ns*nt*nb;
    mwSize f_cartprobsize = f_ns*f_nt*f_nb;
    int np = (int)cartprobsize;
    int nf = (int)f_cartprobsize;
    int i,j,k;

    // Create vars to send to GPU (double2 is max native)
    double2 *pxy;
    double1 *pz;
    double2 *fx, *fy, *fz ;
    double2 *rg;
    double2 *uxy;
    double1 *uz;
    double1 *L_new;
    
    // Test whether this is a GPU enabled simulation  
    int gpusim = 0;
    if(strncmp("true",usegpu,4)==0){
       gpusim = 1;
    }

    // GPU Info
    int threadcount = (int)gpuinfo[0];
    int blockcount = 0;
    int numcards = (int)gpuinfo[1]; // Number of GPUs
    int deviceind = 0; // Only one GPU now, change for multiple
    size_t memfree, memtotal;
    clock_t t0,t1;
     
    // Create vars that will be allocated on GPU
    double2 *dev_pxy;
    double1 *dev_pz;
    double2 *dev_fx,*dev_fy, *dev_fz;
    double2 *dev_rg;
    double1 *dev_len;
    double2 *dev_uxy;
    double1 *dev_uz;
    

    //GPU Error catcher
    cudaError_t cudaError;

    // Allocate memory for local host variables
    pxy = (double2 *)malloc(cartprobsize*sizeof(double2));
    pz = (double1 *)malloc(cartprobsize*sizeof(double1));
    uxy = (double2 *)malloc(cartprobsize*sizeof(double2));
    uz = (double1 *)malloc(cartprobsize*sizeof(double1));
    fx = (double2 *)malloc(f_cartprobsize*sizeof(double2));
    fy = (double2 *)malloc(f_cartprobsize*sizeof(double2));
    fz = (double2 *)malloc(f_cartprobsize*sizeof(double2));
    rg = (double2 *)malloc(f_cartprobsize*sizeof(double2));
    L_new = (double1 *)malloc(cartprobsize*sizeof(double1));
   
    // Pack arrays into double2 structures
    crunch_p(P,pxy,pz,di);
    crunch_f(F1,F2,fx,fy,fz,fdi);
    crunch_rg(gamma,rc,rg,fdi);

        
     
    if(gpusim == 1){
        // Check for any CUDA errors
        checkCUDAError("memcpy");
    
        // Calculate number of blocks
        blockcount = cuda_setup(f_cartprobsize, threadcount, numcards, deviceind);
      
        // Allocate memory on GPU
        cudaError = cudaMalloc( (void**)&dev_fx, f_cartprobsize*sizeof(double2) );
        if (cudaError != cudaSuccess) { mexErrMsgTxt("Out of Nvidia device memory."); }
        cudaError = cudaMalloc( (void**)&dev_fy, f_cartprobsize*sizeof(double2) );
        if (cudaError != cudaSuccess) { mexErrMsgTxt("Out of Nvidia device memory."); }
        cudaError = cudaMalloc( (void**)&dev_fz, f_cartprobsize*sizeof(double2) );
        if (cudaError != cudaSuccess) { mexErrMsgTxt("Out of Nvidia device memory."); }
        cudaError = cudaMalloc( (void**)&dev_len, cartprobsize*sizeof(double1) );
        if (cudaError != cudaSuccess) { mexErrMsgTxt("Out of Nvidia device memory."); }
        cudaError = cudaMalloc( (void**)&dev_pxy, cartprobsize*sizeof(double2) );
        if (cudaError != cudaSuccess) { mexErrMsgTxt("Out of Nvidia device memory."); }
        cudaError = cudaMalloc( (void**)&dev_pz, cartprobsize*sizeof(double1) );
        if (cudaError != cudaSuccess) { mexErrMsgTxt("Out of Nvidia device memory."); }
        cudaError = cudaMalloc( (void**)&dev_rg, f_cartprobsize*sizeof(double2) );
        if (cudaError != cudaSuccess) { mexErrMsgTxt("Out of Nvidia device memory."); }
        cudaError = cudaMalloc( (void**)&dev_uxy, cartprobsize*sizeof(double2) );
        if (cudaError != cudaSuccess) { mexErrMsgTxt("Out of Nvidia device memory."); }
        cudaError = cudaMalloc( (void**)&dev_uz, cartprobsize*sizeof(double1) );
        if (cudaError != cudaSuccess) { mexErrMsgTxt("Out of Nvidia device memory."); }
        
        // Reset GPU memory blocks to 0
        cudaMemset(dev_uxy, 0.00, cartprobsize*sizeof(double2));
        cudaMemset(dev_uz, 0.00, cartprobsize*sizeof(double1));
        cudaMemset(dev_fx, 0.00, f_cartprobsize*sizeof(double2));
        cudaMemset(dev_fy, 0.00, f_cartprobsize*sizeof(double2));
        cudaMemset(dev_fz, 0.00, f_cartprobsize*sizeof(double2));

        // Copy host memory vars to device memory vars
        cudaMemcpy(dev_fx, fx, f_cartprobsize*sizeof(double2), cudaMemcpyHostToDevice);
        cudaMemcpy(dev_fy, fy, f_cartprobsize*sizeof(double2), cudaMemcpyHostToDevice);  
        cudaMemcpy(dev_fz, fz, f_cartprobsize*sizeof(double2), cudaMemcpyHostToDevice); 
        cudaMemcpy(dev_pxy, pxy, cartprobsize*sizeof(double2), cudaMemcpyHostToDevice); 
        cudaMemcpy(dev_pz, pz, cartprobsize*sizeof(double1), cudaMemcpyHostToDevice); 
        cudaMemcpy(dev_rg, rg, f_cartprobsize*sizeof(double2), cudaMemcpyHostToDevice); 
        cudaMemcpyToSymbol(cutoff, &co, sizeof(double));
        cudaMemcpyToSymbol(cored, &d, sizeof(double));
        cudaMemcpyToSymbol(dev_np, &np, sizeof(int));
        cudaMemcpyToSymbol(dev_nf, &nf, sizeof(int));
        cudaMemcpyToSymbol(threadsperblock, &blockcount, sizeof(int));


        // Check for any CUDA errors
        checkCUDAError("memcpy");
       
        // Length calculation (Testing only...not used)
        if(strncmp("leng",casetype,4)==0){
        //mexPrintf("Calling the length cuda kernel!\n");
        dim3 dimGrid(blockcount);
        dim3 dimBlock(threadcount);
            calcLengthOnly<<<dimGrid,dimBlock>>>(dev_fx, dev_fy, dev_fz, dev_len);       
        }
        
        // Biot-Savart naive kernel call
        double c0 = 0; 
        if(strncmp("fuln",casetype,4)==0){
        //mexPrintf("Calling the full cuda kernel!\n");
        dim3 dimGrid(blockcount);
        dim3 dimBlock(threadcount);
        c0 = omp_get_wtime( );
            BiotSavart_naive<<<dimGrid,dimBlock>>>(dev_pxy, dev_pz, dev_fx, dev_fy, dev_fz, dev_rg, dev_uxy, dev_uz);       
         
        }
        
        // Biot-Savart shared memory kernel
        if(strncmp("fuls",casetype,4)==0){
        //mexPrintf("Calling the full cuda kernel!\n");
        dim3 dimGrid(blockcount);
        dim3 dimBlock(threadcount);
        c0 = omp_get_wtime( );
            BiotSavart_sh<<<dimGrid,dimBlock>>>(dev_pxy, dev_pz, dev_fx, dev_fy, dev_fz, dev_rg, dev_uxy, dev_uz);       
         
        }
               
        // block until the device has completed
        cudaThreadSynchronize();
        double c1 = omp_get_wtime( );
        //mexPrintf ("Elapsed wall clock time: %3.9f seconds\n", c1-c0); 

        // Check for any CUDA errors
        checkCUDAError("kernel invocation");
 
        if(strncmp("leng",casetype,4)==0){
        cudaMemcpy(L_new, dev_len, cartprobsize*sizeof(double1), cudaMemcpyDeviceToHost); 
        }

        if(strncmp("ful",casetype,3)==0){
        cudaMemcpy(uxy, dev_uxy, cartprobsize*sizeof(double2), cudaMemcpyDeviceToHost);
        cudaMemcpy(uz, dev_uz, cartprobsize*sizeof(double1), cudaMemcpyDeviceToHost);  
         }


        // Check for any CUDA errors
        checkCUDAError("memcpy");


    } //Close if GPU Sim

    if(gpusim == 0){
         
         
         memset(uxy, 0.00, cartprobsize*sizeof(double2));
         memset(uz, 0.00, cartprobsize*sizeof(double1));

         BiotSavart_nogpu(pxy, pz, fx, fy, fz, rg, uxy, uz, cartprobsize, f_cartprobsize, d, co);
         
         
    }
    
    if(strncmp("leng",casetype,4)==0){
    reassem_array1(L_new,L,di);
    }
    
    if(strncmp("ful",casetype,3)==0){
    reassem_uind(uind,uxy,uz,di);  
    }
    
    
 
    if(gpusim == 1){
    cudaFree(dev_uxy);
    cudaFree(dev_uz);
    cudaFree(dev_fx);
    cudaFree(dev_fy);
    cudaFree(dev_fz);
    cudaFree(dev_len);
    cudaFree(dev_pxy);
    cudaFree(dev_pz);
    cudaFree(dev_rg);
    }
    
    free(pxy);
    free(pz);
    free(uxy);
    free(uz);
    free(fx);
    free(fy);
    free(fz);
    free(L_new);
    free(rg);
    
  
    if(gpusim == 1){
    cuMemGetInfo(&memfree, &memtotal);
    //mexPrintf("Free Memory: %d \n",memfree);
    //mexPrintf("Total Memory: %d \n",memtotal);
    }

    return;
    }

// Function which interfaces with MATLAB
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    mwSize ndims, ndimsin;
    mwSize dims[4],lendims[4], fdims[4];
    const mwSize *indims, *findims;
    double *uind_ptr, *l_ptr;
    double *f1_ptr, *f2_ptr, *p_ptr, *gamma_ptr;
    double *rc_ptr, *co_ptr, *d_ptr;
    double *gpuhw_ptr;
    char *type_ptr, *cm_ptr, *gpuflag_ptr;
    const char *cmnum_ptr;
    mwSize cm_buflen, type_buflen, gpuflag_buflen, cmnum_buflen;

    cm_buflen = mxGetNumberOfElements(CMOD_in) + 1; 
    type_buflen = mxGetNumberOfElements(TYPE_in) + 1; 
    gpuflag_buflen = mxGetNumberOfElements(GPU_in) + 1; 
    cmnum_buflen = cm_buflen - 4; 
 
    cm_ptr = (char *)mxCalloc(cm_buflen, sizeof(char));
    type_ptr = (char *)mxCalloc(type_buflen, sizeof(char));
    gpuflag_ptr = (char *)mxCalloc(gpuflag_buflen, sizeof(char));
    cmnum_ptr = (const char *)mxCalloc(cmnum_buflen, sizeof(char));
    

    if (mxGetString(CMOD_in, cm_ptr, cm_buflen))
        mexErrMsgTxt("Couldn't make Core Model string.");
    if (mxGetString(TYPE_in, type_ptr, type_buflen))
        mexErrMsgTxt("Couldn't make TYPE string.");
    if (mxGetString(GPU_in, gpuflag_ptr, gpuflag_buflen))
        mexErrMsgTxt("Couldn't make GPU string.");

    if (nrhs != 11)
        mexErrMsgTxt("11 inputs expected.");
    if (nlhs != 2)
        mexErrMsgTxt("2 outputs expected."); 
    
    //Copy Core model Number to Int
    
    

    ndimsin = mxGetNumberOfDimensions(P_in);
    //mexPrintf("Num Dims In: %i \n",ndimsin);


    //mexPrintf("1st Type is: %s\n", type_ptr);
    
    if(ndimsin != 4){
      //mexPrintf("Changing Number of Dimensions to 4!\n");
      ndims = 4;
    }else{
    ndims = ndimsin;
    }
    
    indims = mxGetDimensions(P_in);
    findims = mxGetDimensions(F1_in);    

    if(ndimsin == 1){
         mexErrMsgTxt("Need more dimensions in P.");
    }else if(ndimsin == 2){
         dims[0] = indims[0];
         dims[1] = indims[1];
         dims[2] = 1;
         dims[3] = 1;
    }else if(ndimsin == 3){
         dims[0] = indims[0];
         dims[1] = indims[1];
         dims[2] = indims[2];
         dims[3] = 1;
    }else if(ndimsin == 4){
         dims[0] = indims[0];
         dims[1] = indims[1];
         dims[2] = indims[2];
         dims[3] = indims[3]; 
    }else{
         mexErrMsgTxt("Wrong dimensions in P.");
    } 

    fdims[0] = findims[0];
    fdims[1] = findims[1];
    fdims[2] = findims[2];
    fdims[3] = findims[3]; 

    //mexPrintf("Dims: %i %i %i %i \n",dims[0],dims[1],dims[2],dims[3]);
    //mexPrintf("F Dims: %i %i %i %i \n",fdims[0],fdims[1],fdims[2],fdims[3]);
    //mexPrintf("Num Dims: %i \n",ndims);

    // Create a matrix for the return arguments
    UIND_out = mxCreateNumericArray(ndims, dims, mxDOUBLE_CLASS, mxREAL);
    
         lendims[0] = dims[0];
         lendims[1] = 1;
         lendims[2] = dims[2];
         lendims[3] = dims[3]; 
  
    L_out = mxCreateNumericArray(ndims, lendims, mxDOUBLE_CLASS, mxREAL);
   
    // Assign pointers to the various parameters
    uind_ptr = mxGetPr(UIND_out);
    l_ptr = mxGetPr(L_out);     
    f1_ptr = mxGetPr(F1_in); 
    f2_ptr = mxGetPr(F2_in);
    p_ptr = mxGetPr(P_in);
    gamma_ptr = mxGetPr(GAMMA_in);
    rc_ptr = mxGetPr(RC_in);
    d_ptr = mxGetPr(D_in);
    co_ptr = mxGetPr(CO_in);
    gpuhw_ptr = mxGetPr(GPUHW_in);
    

 
        
    // Do the actual computations in the biot subroutine
    biot_calc(uind_ptr,l_ptr,f1_ptr,f2_ptr,p_ptr,
              gamma_ptr,rc_ptr,d_ptr[0],cm_ptr,co_ptr[0],
              type_ptr, ndims, dims, fdims, gpuflag_ptr, gpuhw_ptr); 
    

    
}

    