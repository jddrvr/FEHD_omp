#include <iostream>
#include "FEHD.h"
#include "dataContainers.h"
#include "mkAR.h"
#include "mkARGPU.h"
#include <algorithm>
#include <math.h>
#include <fstream>
//#include "GCopenmp.h"
#include "kernels.h"
#include "GC.h"
#include "utility.h"
#include <chrono>
#include <cblas.h>
#include "timeSeriesOPs.h"

void compGradient(ARmodel A, std::vector<float> &gradient ,std::vector<float> GCvalsBASE,std::vector<float> angleArray,paramContainer params, int numComps,
		  std::vector<int> lagList)
{
  const int numVars = numComps-1;
  const float  h_val = 0.001f; // This is for the gradient spacing.
  
  std::vector<float> angle(angleArray); // Copy this
  std::vector<float> GCvalsUTIL(params.numParticles,0);
  
  for(int varIndex=0;varIndex<numVars;varIndex++)
    {
      //printf("varIndex = %i \n", varIndex);
      
      for(int particle=0;particle<params.numParticles;particle++)
	{
	  angle[particle*numVars+varIndex] += h_val;
	}

      granger(A,angle,GCvalsUTIL,params,numComps,lagList);

      for(int particle=0;particle<params.numParticles;particle++)
	{
	  gradient[particle*numVars+varIndex] = (GCvalsUTIL[particle]-
							GCvalsBASE[particle])/h_val;
	  angle[particle*numVars+varIndex] -= h_val;
	  
	}      
    }  
}


void runFEHD(dataList dataArray, std::vector<float> &Lmat, paramContainer params)
{
  // Set the parameters for sgemm. 
  float alpha=1.0f;
  float beta=0.0f;

  
  std::vector<float> bestAngle;
  matrix Rdecor; // Another example - it's the pca function.
  std::vector<float> Q; // The rotation matrix, resized at each step.
  std::vector<float> T(params.numPCs*params.numPCs,0); // The "work" transformation
  std::vector<float> oneArrayData; // Holds the data without epoch boundaries.
  std::vector<float> transformedData; // Holds the new data without epoch boundaries.  
  std::vector<float> newTrans(params.numPCs*params.numChannels); // Another worker.

  for(int numComps = params.numPCs;numComps>1;numComps--)
    {
      
      Rdecor.elements.clear();
      bestAngle.resize(numComps-1);
      
      runFEHDstep(bestAngle, Rdecor, dataArray, params, numComps);
      //printf("got here \n");
      
      Q.resize(numComps*numComps);

      singleQ(Q,bestAngle);

      // Compute the total transformation for this step
      std::fill(T.begin(),T.end(), 0.0f);

      for(int i=numComps;i<params.numPCs;i++)
	T[i*params.numPCs+i] = 1.0f;

      cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
		  numComps,numComps,numComps,
		  alpha,Q.data(),numComps,
		  Rdecor.elements.data(),numComps,
		  beta,T.data(),params.numPCs);

      oneArrayData.resize(numComps*params.numEpochs*params.epochPts);
      convertDataListToRawArray(dataArray, oneArrayData.data());

      // Transform the data

      
      transformedData.resize(numComps*params.numEpochs*params.epochPts);

      cblas_sgemm(CblasColMajor, CblasNoTrans,CblasNoTrans,
		  numComps,params.numEpochs*params.epochPts,numComps,
		  alpha,T.data(),params.numPCs,
		  oneArrayData.data(), numComps,
		  beta,transformedData.data(), numComps);

      // Convert it back into dataArray

      convertRawArrayToDataList(transformedData.data(), dataArray, numComps, params.epochPts, params.numEpochs);

      // Remove the last component
      removeComponent(dataArray, numComps-1);

      // Update the transformation

      cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
		  params.numPCs,params.numChannels,params.numPCs,
		  alpha, T.data(),params.numPCs,
		  Lmat.data(),params.numPCs,
		  beta, newTrans.data(), params.numPCs);

      Lmat = newTrans;
    }
  
  return;
}

// Runs a step of the algorithm. Randomly initializes the rotation angles,
// rotates the AR model and computes the Granger Causality of the last component
// on those above. The component found to have the minimum influence
// is subsequently removed from the time series. 
void runFEHDstep(std::vector<float> &bestAngle, matrix &L, dataList dataArray ,paramContainer params,int numComps)
{

  // Create random initial conditions
  std::vector<float> angleArray;

  srand((unsigned)time(0));

  for(int indx=0;indx<params.numParticles*(numComps-1);indx++)
    angleArray.push_back((float)(rand()%314-157)/100.0f);

  //angleArray.resize((numComps-1)*params.numParticles);
  //std::fill(angleArray.begin(),angleArray.end(),0.0);
  ARmodel A;
  dataList residuals;
  std::vector<int> lagList(params.lagList); // Just do this for sequential lags for now.

  // Assume you are dealing with a sorted lagList
  //mkAR(dataArray, lagList, A, residuals);


  mkARGPU(dataArray, lagList, A, residuals);

  
  // Need to orthonormalize the residuals.
  dataList ortho_residuals;

  orthonormalizeR(residuals, ortho_residuals, L);

  
  rotate_model(A, L);


  printf("size %i \n",sizeof(int));
  // Some numerical parameters. At this point I am not sure how often opne would need tofiddle with these.

  // Number of iterations
  int numIts = 50;
  std::vector<float> h = {0.01f, 1.0f, 10.0f};
  // Declare a couple of arrays to hold the GCs

  std::vector<float> GCvals(params.numParticles,0);
  std::vector<float> GCvals1(params.numParticles,0);
  std::vector<float> GCvals2(params.numParticles,0);
  std::vector<float> GCvals3(params.numParticles,0);
  
  // And some more angle arrays

  std::vector<float> angleArray1(params.numParticles*(numComps-1),0);
  std::vector<float> angleArray2(params.numParticles*(numComps-1),0);
  std::vector<float> angleArray3(params.numParticles*(numComps-1),0);

  // And a gradient

  std::vector<float> gradient(params.numParticles*(numComps-1),0);

  // Begin the loop by evaluating the GC at angleArray
  //auto t0 = std::chrono::steady_clock::now();
  granger(A,angleArray,GCvals, params,numComps,lagList);
  //exit(0); // Debugging!!!
  //auto t1 = std::chrono::steady_clock::now();

  //std::chrono::duration<float> time_span = std::chrono::duration_cast<std::chrono::duration<float>>(t1-t0);
  
  //printf("time elapsed = %f \n",time_span.count());
  //printf("The smallest element is %f \n",*std::min_element(GCvals.begin(),GCvals.end()));

  std::vector<float> candidates(4,0);
  int minIndx;
  for(int iter=0;iter<numIts;iter++)
    {
      //auto t2 = std::chrono::steady_clock::now();
      compGradient(A, gradient, GCvals, angleArray, params, numComps,lagList);
      //auto t3 = std::chrono::steady_clock::now();
      //std::chrono::duration<float> time_span2 = std::chrono::duration_cast<std::chrono::duration<float>>(t3-t2);
  
      //printf("time elapsed = %f \n",time_span2.count());
      
      angleArray1 = angleArray;
      angleArray2 = angleArray;
      angleArray3 = angleArray;
      cblas_saxpy(params.numParticles*(numComps-1), -h[0], gradient.data(), 1, angleArray1.data(),1);
      cblas_saxpy(params.numParticles*(numComps-1), -h[1], gradient.data(), 1, angleArray2.data(),1);
      cblas_saxpy(params.numParticles*(numComps-1), -h[2], gradient.data(), 1, angleArray3.data(),1);

      // Evaluate GC for each of the shifted angle arrays
      granger(A,angleArray1,GCvals1,params,numComps,lagList);
      granger(A,angleArray2,GCvals2,params,numComps,lagList);
      granger(A,angleArray3,GCvals3,params,numComps,lagList);
      
      // This loop determines the minimum of four values for each particle.
      for(int particle=0;particle<params.numParticles;particle++)
	{

	  // Check for negative values here - if it is negative just
	  // choose a new angle and disregard. 
	  // Some error checking here would be good. For example, make sure the GC is positive
	  // for each particle. 
	  candidates[0] = GCvals[particle];
	  candidates[1] = GCvals1[particle];
	  candidates[2] = GCvals2[particle];
	  candidates[3] = GCvals3[particle];

	  minIndx = std::distance(candidates.begin(),min_element(candidates.begin(),candidates.end()));

	  // Here, if minIndx = 0, I am going to include a recycler, which resetS the angle so the search can be restarted.
	  // This could be condensed.
	  if(minIndx == 1)
	    {
	      GCvals[particle] = GCvals1[particle];
	      std::copy(angleArray1.data()+particle*(numComps-1),angleArray1.data()+particle*(numComps-1)+numComps-1,angleArray.data()+particle*(numComps-1));
	    }
	  if(minIndx == 2)
	    {
	      GCvals[particle] = GCvals2[particle];
	      std::copy(angleArray2.data()+particle*(numComps-1),angleArray2.data()+particle*(numComps-1)+numComps-1,angleArray.data()+particle*(numComps-1));
	    }
	  if(minIndx == 3)
	    {
	      GCvals[particle] = GCvals3[particle];
	      std::copy(angleArray3.data()+particle*(numComps-1),angleArray3.data()+particle*(numComps-1)+numComps-1,angleArray.data()+particle*(numComps-1));
	    }
	}
      
      //if(verbose)
      	printf("iteration = %i, particle = %li, value = %e \n",
             iter,std::min_element(GCvals.begin(),GCvals.end())-GCvals.begin(),
             GCvals[std::min_element(GCvals.begin(),GCvals.end())-GCvals.begin()]);
  
    }
  
  // Return the best angle.

  long unsigned int indexVal = std::min_element(GCvals.begin(),GCvals.end())-GCvals.begin();

  //printf("%li \n",indexVal);

  std::copy(angleArray.data()+indexVal*(numComps-1),angleArray.data()+indexVal*(numComps-1)+numComps-1,bestAngle.begin());

  return;
 
}



void singleQ(std::vector<float> &Q, std::vector<float> angle)
{

  int numVars = angle.size();
  
  float Qcol1;
  float Qcol2;

  float sinVal;
  float cosVal;

  // Create an identity matrix.
  for(int indx=0;indx<(numVars+1)*(numVars+1);indx++)
    {
      Q[indx] = 0;
    }
  for(int row=0;row<(numVars+1);row++)
    {
      Q[row*(numVars+1)+row] = 1.0f;
    }

  // Multiply the individual rotations together. 
  for(int varIndx=0; varIndx<numVars; varIndx++) // Cycle through the angles. 
    {
      sinVal = sinf(angle[varIndx]);// Assign the cos and sin to variables.
      cosVal = cosf(angle[varIndx]);

      for(int k=0;k<numVars+1;k++) // Do the matrix multiplication 
	{
	  Qcol1 = Q[k*(numVars+1)+varIndx]; //Q(particle,row i, column k
	  Qcol2 = Q[k*(numVars+1)+numVars];//(particle,row M-1,column k
	  Q[k*(numVars+1)+varIndx] = cosVal*Qcol1-sinVal*Qcol2;
	  Q[k*(numVars+1)+numVars] = sinVal*Qcol1+cosVal*Qcol2;
	}

    }
  
  return;
  
}
