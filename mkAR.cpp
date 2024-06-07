#include "timeSeriesOPs.h"
#include "dataContainers.h"
#include "mkAR.h"
#include <vector>
#include <cblas.h>
#include <lapacke.h>
#include <algorithm>


void mkAR(dataList dataArray, std::vector<int> lagList, ARmodel &A, dataList &R)
{
  int numEpochs = dataArray.epochArray.size();
  int epochPts = dataArray.epochArray[0].timePointArray.size();
  int numComps = dataArray.epochArray[0].timePointArray[0].dataVector.size();
  int numPoints = numEpochs*epochPts;
  
  matrix LHS;
  matrix RHS;

  // We want maxLag, not numLags
  
  int maxLag = *std::max_element(lagList.begin(),lagList.end());
  int tpoint;
  // Sort the lag list in reverse order, so everything matches up without changing a lot.

  std::sort(lagList.begin(),lagList.end(), std::greater<int>());
  
  for(int epoch=0;epoch<numEpochs;epoch++)
      {
      for(int tp=maxLag;tp<epochPts;tp++)
	{
	  RHS.elements.insert(RHS.elements.end(),
			      dataArray.epochArray[epoch].timePointArray[tp].dataVector.begin(),
			      dataArray.epochArray[epoch].timePointArray[tp].dataVector.end());
	}

      for(int tp=0;tp<epochPts-maxLag;tp++)
	{
	  for(int lagIndx=0;lagIndx<lagList.size();lagIndx++)
	    {
	      tpoint = tp + maxLag-lagList[lagIndx];
	      LHS.elements.insert(LHS.elements.end(),
				  dataArray.epochArray[epoch].timePointArray[tpoint].dataVector.begin(),
				  dataArray.epochArray[epoch].timePointArray[tpoint].dataVector.end());
	    }
	}

      }

  // Make some copies (probably some of them are unnecessary.
  std::vector<float> LHSvec(LHS.elements);
  std::vector<float> RHSvec(RHS.elements);
  std::vector<float> LHSvecBAK(LHS.elements);
  std::vector<float> RHSvecBAK(RHS.elements);

  int numLags = lagList.size();
  
  int mval = numLags*numComps;
  int nval = numLags*numComps;
  int kval = (epochPts-maxLag)*numEpochs;

  const float alpha = 1.0f;
  const float beta = 0.0f;

  float *LHScov = new float[numComps*numLags*numComps*numLags];
  float *RHScov = new float[numComps*numLags*numComps];

  // This can be altered to use a rank update blas routine
  cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans,
	      mval, nval, kval, alpha, LHSvec.data(), mval, LHSvec.data(), mval, beta, LHScov, mval);
  cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans,
	      mval, numComps, kval, alpha, LHSvec.data(), mval, RHSvec.data(), numComps, beta, RHScov, mval);

  // Solve LHScov*A=RHScov for A
  // LHScov is a symmetric matrix, so we can use a non-ge lapack routine.
  int info;
  int IPIV[numComps*numLags];
  
  info=LAPACKE_ssysv(LAPACK_COL_MAJOR, 'U', numComps*numLags, numComps, LHScov, numComps*numLags,IPIV,
		     RHScov,numComps*numLags);

  matrix lagTmp;

  for(int lag=numLags-1;lag>=0;lag--)
    {
      for(int col=0;col<numComps;col++)
	for(int row=0;row<numComps;row++)
	  {
	    // Transpose and copy.
	    lagTmp.elements.push_back(RHScov[row*numLags*numComps+col+lag*numComps]);
	  }

      A.lagMatrices.push_back(lagTmp);
      lagTmp.elements.clear();
    }

  // Compute the residuals.
  //float *modelvals = new float[numComps*(epochPts-numLags)*numEpochs];
  //cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
  //	      numComps,kval,mval,alpha,RHScov,mval,LHSvecBAK.data(),mval,beta,modelvals,numComps);
  const float alphaRes=-1.0;
  const float betaRes=1.0;

  cblas_sgemm(CblasColMajor,CblasTrans, CblasNoTrans,
	      numComps,kval,mval,alphaRes,RHScov,mval,LHSvecBAK.data(),mval,betaRes,RHSvecBAK.data(),numComps);
  //const float alphaneg = -1.0f;
  //int mval2 = RHS.elements.size();
  //cblas_saxpy(mval2, alphaneg, modelvals, 1, RHSvecBAK.data(), 1); // I think this is done right, just long.
  convertRawArrayToDataList(RHSvecBAK.data(), R, numComps, epochPts-maxLag, numEpochs);
  
  //  delete [] modelvals;
  delete [] RHScov;
  delete [] LHScov;
  
  
  
  return;
}

void orthonormalizeR(dataList residuals, dataList &ortho_residuals, matrix &L)
{
  PCA(residuals, ortho_residuals, L);
}

void rotate_model(ARmodel &A, matrix L)
{

  int M = sqrt(L.elements.size()); // Matrix dimension
  int numLags = A.lagMatrices.size(); // number of lags
  int info1,info2; // Error checking
  matrix LBAK; // The original will be changing, it is easiest to just make a copy. 
  LBAK.elements = L.elements;

  // Invert the matrix
  // Uses LAPACK to LU=A (trf) and back solve (tri)
  std::vector<int> ipiv(M,0);
  
  info1 = LAPACKE_sgetrf(LAPACK_COL_MAJOR,M,M,L.elements.data(),M,ipiv.data());
  info2 = LAPACKE_sgetri(LAPACK_COL_MAJOR,M,L.elements.data(),M,ipiv.data());
  if(info1 != 0 || info2 != 0)
    {
      // Put some diagnostics here - 

      printf("Error inverting the residual transformation matrix \n");
      exit(0);
    }
  /*  printf("The inverse \n");
  for(int row=0;row<M;row++)
    {
      for(int col=0; col<M;col++)
	{
	  printf("%f ",L.elements[col*M+row]);

	}
      printf("\n");
    }
  */
  // Multiply L A L^-1 for each lag matrix in the AR model
  const float alpha=1.0f;
  const float beta=0.0f;

  std::vector<float> tmp(M*M,0);
  for(int lag=0; lag<numLags; lag++)
    {
      
      cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
		  M,M,M,alpha,LBAK.elements.data(), M,
		  A.lagMatrices[lag].elements.data(),M,
		  beta, tmp.data(),M);
      

      cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
		  M,M,M, alpha,tmp.data(),M,
		  L.elements.data(),M,
		  beta, A.lagMatrices[lag].elements.data(),M);
    }
  
  
  return;
}
		  
  
