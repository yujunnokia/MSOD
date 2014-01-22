#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>

const double EPSILON = 1e-10;

// test
void hello(int *n)
{
    int i;
    for(i=0; i < *n; i++) {
        Rprintf("Hello, world!\n");
    }
}

void printMatrix(double *n, int *nElements)
{
    int i;
    for(i=0; i < *nElements; i++) {
        Rprintf("%f,",n[i]);
    }
    Rprintf("\n");
}

// 
void kernel_smooth(double *x, int *n, double *xpts, int *nxpts,
                   double *h, double *result)
{
    int i, j;
    double d, ksum;
    for(i=0; i < *nxpts; i++) {
        ksum = 0;
        for(j=0; j < *n; j++) {
            d = xpts[i] - x[j];
            ksum += dnorm(d / *h, 0, 1, 0);
        }
        result[i] = ksum / ((*n) * (*h));
    }
}

// logistic function
double logistic(double t)
{
    return( 1 / ( 1 + exp(-t)) );
}

// vector inner product
double innerProduct(double *a, double *b, int n)
{
    double sum = 0;
    for (int i = 0; i < n; i++) {
        sum = sum+a[i]*b[i];
    }
    return(sum);
}
double innerProductInt(int *a, double *b, int n)
{
    double sum = 0;
    for (int i = 0; i < n; i++) {
        sum = sum+a[i]*b[i];
    }
    return(sum);
}

// save multiple
double multiple(double a, double b)
{
    if ((a == 0 && b == -INFINITY) || (a == -INFINITY && b == 0)) {
        return(0);
    }
    return(a*b);
}

// compute idx from Z
int findIdx(int *Z, int nZ)
{
    int sum = 0;
    for (int i = 0; i < nZ; i++)
    {
        sum = sum + Z[i] * pow(2, nZ-i-1);
    }
    return(sum);
}

// compute probability of detection for NoisyOR model
double ComputeProbDet(int s, int *nSpecies, int *nDetCovs, double *betas, double *beta0s, double *Xdit, int *Z, int *confusion, double isLeak) 
{
    double probDet = 1;
    
    double *beta0 = beta0s + s*(*nDetCovs);
    probDet = probDet * pow(1-logistic(innerProduct(Xdit,beta0,*nDetCovs)), isLeak);
    
    for (int ss = 0; ss < *nSpecies; ss++) {
        double *beta = betas+s*(*nSpecies)*(*nDetCovs)+ss*(*nDetCovs);
        probDet = probDet * pow(1-logistic(innerProduct(Xdit,beta,*nDetCovs)), Z[ss]*confusion[s*(*nSpecies)+ss]);
    } // ss
    probDet = 1-probDet;
        
    if (probDet == 0) { probDet = EPSILON; }
    if (probDet == 1) { probDet = 1-EPSILON; }
    
    return probDet;
}

///////////////////////////////////////////////
// MSOD model with Noisy-OR parameterization //
///////////////////////////////////////////////

// compute log likelihood of site i for Noisy-OR model
void ComputeiLLNoisyOR(int *nOccCovs, int *nDetCovs, int *nSpecies, double *alphas, double *betas, double *beta0s, int* confusion, int *Yi, double *Xoi, double *Xdi, int *visit, int *Zs, int *nZs, int *leak, double* iLL)
{
    *iLL = 0;
    
    // check if leak is true
    double isLeak = 0;
    if (leak[0] == 1 || leak[0] == 2) {
        isLeak = 1;
    }
    
    // compute occ prob for all species
    double probOccs[1024];
    for (int s = 0; s < (*nSpecies); s++) {
        double *alpha = alphas + s * (*nOccCovs);
        probOccs[s] = logistic(innerProduct(Xoi,alpha,*nOccCovs));
    } // s
        
    // sum over all possible Zs
    double iL = 0;
    for (int k = 0; k < *nZs; k++) {
        int *Z = Zs + k * (*nSpecies);
        
        double isL = 1;
        for (int s = 0; s < *nSpecies; s++) {
            isL = isL * pow(probOccs[s], Z[s]) * pow(1-probOccs[s], 1-Z[s]);
            double *beta0 = beta0s+s*(*nDetCovs);
            
            for (int t = 0; t < *visit; t++) {
                double *Xdit = Xdi + t*(*nDetCovs);
                
                // compute det prob
                double probDet = ComputeProbDet(s,nSpecies,nDetCovs,betas,beta0s,Xdit,Z,confusion,isLeak);

                /*
                double probDet = 1;
                probDet = probDet * pow(1-logistic(innerProduct(Xdit,beta0,*nDetCovs)), isLeak);
                for (int ss = 0; ss < *nSpecies; ss++) {
                    double *beta = betas+s*(*nSpecies)*(*nDetCovs)+ss*(*nDetCovs);
                    probDet = probDet * pow(1-logistic(innerProduct(Xdit,beta,*nDetCovs)), Z[ss]*confusion[s*(*nSpecies)+ss]);
                } // ss
                probDet = 1-probDet;
                if (probDet == 0) { probDet = EPSILON; }
                if (probDet == 1) { probDet = 1-EPSILON; }
                */
                isL = isL * pow(probDet, Yi[t*(*nSpecies)+s]) * pow(1-probDet, 1-Yi[t*(*nSpecies)+s]);
            } // t
        } // s
        iL = iL + isL;
    } // k
    *iLL = log(iL);
} // ComputeiLLNoisyOR

// compute expected joint log-likelihood
void ComputeEJLLNoisyOR(int *nSites, int *nVisits, int *nOccCovs, int *nDetCovs, int *nSpecies, double *alphas, double *betas, double *beta0s, double *probExpectedOccs, int *confusion, int *Y, double *Xo, double *Xd, int *visits, int *Zs, int *nZs, int *leak, double *EJLL)
{
    *EJLL = 0;
    
    // check if leak is true
    double isLeak = 0;
    if (leak[0] == 1 || leak[0] == 2) {
        isLeak = 1;
    }
    
    for (int i = 0; i < *nSites; i++) {
        double *Xoi = Xo+i*(*nOccCovs);
        
        // compute occ prob for all species
        double probOccs[1024];
        for (int s = 0; s < (*nSpecies); s++) {
            double *alpha = alphas + s * (*nOccCovs);
            probOccs[s] = logistic(innerProduct(Xoi,alpha,*nOccCovs));
        } // s
        
        // compute EJLL for site i
        double iEJLL = 0;
        for (int k = 0; k < *nZs; k++) {
            int *Z = Zs + k * (*nSpecies);
            int idx = findIdx(Z,*nSpecies);
            double probZ = probExpectedOccs[idx*(*nSites) + i];
            
            if (probZ == 0) { continue; }
            
            double isEJLL = 0;
            for (int s = 0; s < *nSpecies; s++) {
                double *beta0 = beta0s+s*(*nDetCovs);
                
                isEJLL = isEJLL + multiple(Z[s], log(probOccs[s])) + multiple((1-Z[s]), log(1-probOccs[s]));
                for (int t = 0; t < visits[i]; t++) {
                    double *Xdit = Xd + i*(*nVisits)*(*nDetCovs) + t*(*nDetCovs);
                    int Yits = Y[i*(*nVisits)*(*nSpecies) + t*(*nSpecies) + s];
                                      
                    // compute det prob
                    double probDet = ComputeProbDet(s,nSpecies,nDetCovs,betas,beta0s,Xdit,Z,confusion,isLeak);
                    isEJLL = isEJLL + multiple(Yits, log(probDet)) + multiple((1-Yits), log(1-probDet));
                } // t
            } // s
            
            iEJLL = iEJLL + multiple(probZ, isEJLL);
        } // k
        
        *EJLL = *EJLL + iEJLL;
    } // i
} // ComputeEJLLNoisyOR


// compute derivs of EJLL for NoisyOR model
void ComputeDerivsOfEJLLNoisyOR(int *nSites, int *nVisits, int *nOccCovs, int *nDetCovs, int *nSpecies, double *alphas, double *betas, double *beta0s, double *probExpectedOccs, int *confusion, int *Y, double *Xo, double *Xd, int *visits, int *Zs, int *nZs, int *leak, double *derivsOfEJLL)
{
    double dQdalpha[1024];
    double dQdbeta0[1024];
    double dQdbeta[4096];
    
    for (int i = 0; i < 1024; i++) {
        dQdalpha[i] = 0;
        dQdbeta0[i] = 0;
    }
    for (int i = 0; i < 4096; i++) {
        dQdbeta[i] = 0;
    }
    
    // check if leak is true
    double isLeak = 0;
    if (leak[0] == 1 || leak[0] == 2) {
        isLeak = 1;
    }

    for (int s = 0; s < (*nSpecies); s++) {
        double *alpha = alphas + s * (*nOccCovs);
        double *beta0 = beta0s + s * (*nDetCovs);
        
        // compute derivs of alpha
        for (int i = 0; i < (*nSites); i++) {
            double *Xoi = Xo+i*(*nOccCovs);
            // compute prob occ
            double probOcc = logistic(innerProduct(Xoi,alpha,*nOccCovs));
            
            double probExpectedOccis = 0;
            for (int k = 0; k < (*nZs); k++) {
                int *Z = Zs + k * (*nSpecies);
                
                if (Z[s] == 0) { continue; }
                int idx = findIdx(Z,*nSpecies);
                probExpectedOccis = probExpectedOccis + probExpectedOccs[idx*(*nSites)+i];
            } // k
            for (int j = 0; j < (*nOccCovs); j++) {
                dQdalpha[s*(*nOccCovs)+j] = dQdalpha[s*(*nOccCovs)+j] + (probExpectedOccis - probOcc) * Xoi[j];
            } // j
        } // i

        
        // if leak probabality is on
        if (leak[0] == 1 || leak[0] == 2) {
            
            // compute derivs of beta0
            for (int i = 0; i < (*nSites); i++) {
                for (int k = 0; k < (*nZs); k++) {
                    int *Z = Zs + k * (*nSpecies);
                    int idx = findIdx(Z,*nSpecies);
                    double probZ = probExpectedOccs[idx*(*nSites) + i];
                    
                    double term[1024];
                    for (int j = 0; j < 1024; j++) { term[j] = 0; }
                    for (int t = 0; t < visits[i]; t++) {
                        double *Xdit = Xd + i*(*nVisits)*(*nDetCovs) + t*(*nDetCovs);
                        int Yits = Y[i*(*nVisits)*(*nSpecies) + t*(*nSpecies) + s];
                        
                        // compute det prob
                        double probDet = ComputeProbDet(s,nSpecies,nDetCovs,betas,beta0s,Xdit,Z,confusion,isLeak);
                        double probDetss = logistic(innerProduct(Xdit,beta0,*nDetCovs));
                        if (Yits == 1) {
                            for (int j = 0; j < (*nDetCovs); j++) {
                                term[j] = term[j] + (1-probDet) / probDet * probDetss * Xdit[j];
                            }
                        } else {
                            for (int j = 0; j < (*nDetCovs); j++) {
                                term[j] = term[j] - probDetss * Xdit[j];
                            }
                        }
                    } // t
                    
                    for (int j = 0; j < (*nDetCovs); j++) {
                        dQdbeta0[s*(*nDetCovs) + j] = dQdbeta0[s*(*nDetCovs) + j] + probZ * term[j];
                    } // j
                } // k
            } // i
        }  // if
        
        
        // compute derivs of beta
        for (int ss = 0; ss < (*nSpecies); ss++) {
            if (confusion[s*(*nSpecies)+ss] == 0) { continue; }
            for (int i = 0; i < (*nSites); i++) {
                for (int k = 0; k < (*nZs); k++) {
                    int *Z = Zs + k * (*nSpecies);
                    int idx = findIdx(Z,*nSpecies);
                    double probZ = probExpectedOccs[idx*(*nSites) + i];
                                        
                    if (Z[ss] == 0) { continue; }
                    
                    double term[1024];
                    for (int j = 0; j < 1024; j++) { term[j] = 0; }
                    for (int t = 0; t < visits[i]; t++) {
                        double *Xdit = Xd + i*(*nVisits)*(*nDetCovs) + t*(*nDetCovs);
                        int Yits = Y[i*(*nVisits)*(*nSpecies) + t*(*nSpecies) + s];
                        
                        // compute det prob
                        double probDet = ComputeProbDet(s,nSpecies,nDetCovs,betas,beta0s,Xdit,Z,confusion,isLeak);
                        
                        double *beta = betas+s*(*nSpecies)*(*nDetCovs)+ss*(*nDetCovs);
                        double probDetss = logistic(innerProduct(Xdit,beta,*nDetCovs));
                        
                        for (int j = 0; j < (*nDetCovs); j++) {
                            term[j] = term[j] + (1/probDet * Yits - 1) * probDetss * Xdit[j] * confusion[s*(*nSpecies)+ss];
                        }
                    } // t

                    for (int j = 0; j < (*nDetCovs); j++) {
                        dQdbeta[s*(*nSpecies)*(*nDetCovs) + ss*(*nDetCovs) + j] = dQdbeta[s*(*nSpecies)*(*nDetCovs) + ss*(*nDetCovs) + j] + probZ * term[j];
                    } // j
                } // k
            } // i
        } // ss
        
    } // s
    
    int count = 0;
    for (int s = 0; s < (*nSpecies); s++) {
        for (int j = 0; j < (*nOccCovs); j++) {
            derivsOfEJLL[count+j] = dQdalpha[s*(*nOccCovs)+j];
        }
        count = count + (*nOccCovs);
        
        if (leak[0] == 1) {
            for (int j = 0; j < (*nDetCovs); j++) {
                derivsOfEJLL[count+j] = dQdbeta0[s*(*nDetCovs)+j];
            }
            count = count + (*nDetCovs);
        }
        if (leak[0] == 2) {
            derivsOfEJLL[count] = dQdbeta0[s*(*nDetCovs)];
            count = count + 1;
        }
        
        for (int ss = 0; ss < (*nSpecies); ss++) {
            if (confusion[s*(*nSpecies)+ss] == 0) { continue; }
            for (int j = 0; j < (*nDetCovs); j++) {
                derivsOfEJLL[count+j] = dQdbeta[s*(*nSpecies)*(*nDetCovs)+ss*(*nDetCovs)+j];    
            }
            count = count + (*nDetCovs);
        } // ss
    } // s

} // ComputeDerivsOfEJLL


// Predict Detection
void PredictDetNoisyOR(int *nSites, int *nVisits, int *nOccCovs, int *nDetCovs, int *nPa, double *alphas, double *betas, double *beta0s, int *pi, double *Xo, double *Xd, int *visits, int *Zs, int *nZs, int *leak, double *predictedDet)
{
    // check if leak is true
    double isLeak = 0;
    if (leak[0] == 1 || leak[0] == 2) {
        isLeak = 1;
    }
    
    for (int i = 0; i < (*nSites); i++) {
        double *Xoi = Xo+i*(*nOccCovs);

        for (int t = 0; t < visits[i]; t++) {
            double *Xdit = Xd + i*(*nVisits)*(*nDetCovs) + t*(*nDetCovs);
            
            predictedDet[i*(*nVisits)+t] = 0;
            for (int k = 0; k < (*nZs); k++) {
                int *Z = Zs + k * (*nPa);
                
                // compute occ prob
                double probZOcc = 1;
                for (int s = 0; s < (*nPa); s++) {
                    double *alpha = alphas + s * (*nOccCovs);
                    double probOcc = logistic(innerProduct(Xoi,alpha,*nOccCovs));
                    probZOcc = probZOcc * pow(probOcc, Z[s]) * pow(1-probOcc, 1-Z[s]);
                } // s
                
                // compute det prob
                double probZDet = 1;
                probZDet = probZDet * pow(1-logistic(innerProduct(Xdit,beta0s,*nDetCovs)), isLeak);
                for (int s = 0; s < (*nPa); s++) {
                    double *beta = betas+s*(*nDetCovs);
                    probZDet = probZDet * pow(1-logistic(innerProduct(Xdit,beta,*nDetCovs)), Z[s]);
                } // s
                probZDet = 1 - probZDet;
                
                predictedDet[i*(*nVisits)+t] = predictedDet[i*(*nVisits)+t] + probZOcc*probZDet;
            } // k
        } // t
    } // i

} // PredictDet


///////////////////////////////////////////////
// MSOD model with Additive parameterization //
///////////////////////////////////////////////


// compute log likelihood of site i for Additive model
// needs to update
void ComputeiLLAdditive(int *nOccCovs, int *nDetCovs, int *nSpecies, double *alphas, double *betas, int* confusion, int *Yi, double *Xoi, double *Xdi, int *visit, int *Zs, int *nZs, double* iLL)
{
    *iLL = 0;
    
    // compute occ prob for all species
    double probOccs[1024];
    for (int s = 0; s < (*nSpecies); s++) {
        double *alpha = alphas + s * (*nOccCovs);
        probOccs[s] = logistic(innerProduct(Xoi,alpha,*nOccCovs));
    } // s
    
    // sum over all possible Zs
    double iL = 0;
    for (int k = 0; k < *nZs; k++) {
        int *Z = Zs + k * (*nSpecies);
        
        double isL = 1;
        for (int s = 0; s < *nSpecies; s++) {
            isL = isL * pow(probOccs[s], Z[s]) * pow(1-probOccs[s], 1-Z[s]);
            double *beta = betas+s*((*nDetCovs)+(*nSpecies));
            double *gamma = beta+(*nDetCovs);
            
            for (int t = 0; t < *visit; t++) {
                // compute det prob
                double probDet = logistic(innerProduct(Xdi+t*(*nDetCovs),beta,*nDetCovs)+innerProduct((double*)Z,gamma,*nSpecies));
                if (probDet == 0) { probDet = EPSILON; }
                if (probDet == 1) { probDet = 1-EPSILON; }
                
                isL = isL * pow(probDet, Yi[t*(*nSpecies)+s]) * pow(1-probDet, 1-Yi[t*(*nSpecies)+s]);
            } // t
        } // s
        
        iL = iL + isL;
    } // k
    
    *iLL = log(iL);
} // ComputeiLLAdditive

// compute expected joint log-likelihood
// needs to update
void ComputeEJLLAdditive(int *nSites, int *nVisits, int *nOccCovs, int *nDetCovs, int *nSpecies, double *alphas, double *betas, double *probExpectedOccs, int *confusion, int *Y, double *Xo, double *Xd, int *visits, int *Zs, int *nZs, int *regType, double *lambdaO, double *lambdaD, double *EJLL)
{
    *EJLL = 0;
    
    for (int i = 0; i < *nSites; i++) {
        double *Xoi = Xo+i*(*nOccCovs);
        
        // compute occ prob for all species
        double probOccs[1024];
        for (int s = 0; s < (*nSpecies); s++) {
            double *alpha = alphas + s * (*nOccCovs);
            probOccs[s] = logistic(innerProduct(Xoi,alpha,*nOccCovs));
        } // s
        
        // compute EJLL for site i
        double iEJLL = 0;
        for (int k = 0; k < *nZs; k++) {
            int *Z = Zs + k * (*nSpecies);
            int idx = findIdx(Z,*nSpecies);
            double probZ = probExpectedOccs[idx*(*nSites) + i];
            
            double isEJLL = 0;
            for (int s = 0; s < *nSpecies; s++) {
                if (Z[s] != 0 || log(probOccs[s]) != -INFINITY) {
                    isEJLL = isEJLL + Z[s] * log(probOccs[s]);
                }
                if ((1-Z[s]) != 0 || log(1-probOccs[s]) != -INFINITY) {
                    isEJLL = isEJLL + (1-Z[s]) * log(1-probOccs[s]);
                }
                //                isEJLL = isEJLL + Z[s] * log(probOccs[s]) + (1-Z[s]) * log(1-probOccs[s]);
                
                double *beta = betas+s*(*nDetCovs+*nSpecies);
                double *gamma = betas+s*(*nDetCovs+*nSpecies)+(*nDetCovs);
                for (int t = 0; t < visits[i]; t++) {
                    double *Xdit = Xd + i*(*nVisits)*(*nDetCovs) + t*(*nDetCovs);                    
                    double probDet = logistic(innerProduct(Xdit,beta,*nDetCovs)+innerProduct((double*)Z,gamma,*nSpecies));
                    if (probDet == 0) { probDet = 1e-5; }
                    if (probDet == 1) { probDet = 1-1e-5; }
                    
                    int Yits = Y[i*(*nVisits)*(*nSpecies) + t*(*nSpecies) + s];
                    
                    if (Yits != 0 || log(probDet) != -INFINITY) {
                        isEJLL = isEJLL + Yits * log(probDet);    
                    }
                    if ((1-Yits) != 0 || log(1-probDet) != -INFINITY) {
                        isEJLL = isEJLL + (1-Yits) * log(1-probDet);
                    }
                    //                    isEJLL = isEJLL + Yits * log(probDet) + (1-Yits) * log(1-probDet);
                } // t
            } // s
            
            if (probZ != 0 || isEJLL != -INFINITY) {
                iEJLL = iEJLL + probZ * isEJLL;
            }
        } // k
        
        *EJLL = *EJLL + iEJLL;
    } // i
    
    // regularization
    if (regType[0] == 2) {
        double regO = 0;
        for (int i = 0; i < (*nSpecies)*(*nOccCovs); i++) {
            if (i % (*nOccCovs) == 0) { continue; }
            regO = regO + pow(alphas[i],2);
        }
        *EJLL = *EJLL - 0.5 * lambdaO[0] * regO;
        
        double regD = 0;
        for (int s = 0; s < (*nSpecies); s++) {
            for (int j = 1; j < (*nDetCovs)+(*nSpecies); j++) {
                regD = regD + pow(betas[s*((*nDetCovs) + (*nSpecies)) + j],2);
            } // j
        } // s
        *EJLL = *EJLL - 0.5 * lambdaD[0] * regD;
    }
    
    *EJLL = - (*EJLL);
} // ComputeEJLLAdditive

// compute derivs of EJLL
// needs to update
void ComputeDerivsOfEJLLAdditive(int *nSites, int *nVisits, int *nOccCovs, int *nDetCovs, int *nSpecies, double *alphas, double *betas, double *probExpectedOccs, int *confusion, int *Y, double *Xo, double *Xd, int *visits, int *Zs, int *nZs, int *regType, double *lambdaO, double *lambdaD, double *derivsOfEJLL)
{
    double dQdalpha[1024];
    double dQdbeta[4096];
    
    for (int i = 0; i < 1024; i++) {
        dQdalpha[i] = 0;
    }
    for (int i = 0; i < 4096; i++) {
        dQdbeta[i] = 0;
    }
    
    for (int s = 0; s < (*nSpecies); s++) {
        double *alpha = alphas + s * (*nOccCovs);
        
        // compute derivs of alpha
        for (int i = 0; i < (*nSites); i++) {
            double *Xoi = Xo+i*(*nOccCovs);
            // compute prob occ
            double probOcc = logistic(innerProduct(Xoi,alpha,*nOccCovs));
            
            double probExpectedOccis = 0;
            for (int k = 0; k < (*nZs); k++) {
                int *Z = Zs + k * (*nSpecies);
                
                if (Z[s] == 0) { continue; }
                int idx = findIdx(Z,*nSpecies);
                probExpectedOccis = probExpectedOccis + probExpectedOccs[idx*(*nSites)+i];
            } // k
            for (int j = 0; j < (*nOccCovs); j++) {
                dQdalpha[s*(*nOccCovs)+j] = dQdalpha[s*(*nOccCovs)+j] + (probExpectedOccis - probOcc) * Xoi[j];
            } // j
        } // i
        
        // compute derivs of beta
        for (int i = 0; i < (*nSites); i++) {
            for (int k = 0; k < (*nZs); k++) {
                int *Z = Zs + k * (*nSpecies);
                int idx = findIdx(Z,*nSpecies);
                double probZ = probExpectedOccs[idx*(*nSites) + i];
                
                double term[1024];
                for (int j = 0; j < 1024; j++) { term[j] = 0; }
                
                double *beta = betas+s*(*nDetCovs+*nSpecies);
                double *gamma = betas+s*(*nDetCovs+*nSpecies)+(*nDetCovs);
                for (int t = 0; t < visits[i]; t++) {
                    double *Xdit = Xd + i*(*nVisits)*(*nDetCovs) + t*(*nDetCovs);     
                    double probDet = logistic(innerProduct(Xdit,beta,*nDetCovs)+innerProduct((double*)Z,gamma,*nSpecies));
                    if (probDet == 0) { probDet = 1e-5; }
                    if (probDet == 1) { probDet = 1-1e-5; }
                    //                        Rprintf("probDet is %f.\n",probDet);
                    int Yits = Y[i*(*nVisits)*(*nSpecies) + t*(*nSpecies) + s];
                    
                    for (int j = 0; j < (*nDetCovs); j++) {
                        term[j] = term[j] + (Yits - probDet) * Xdit[j];
                    } // j
                    for (int j = 0; j < (*nSpecies); j++) {
                        term[(*nDetCovs)+j] = term[(*nDetCovs)+j] + (Yits - probDet) * Z[j];
                    }
                    
                } // t
                
                for (int j = 0; j < (*nDetCovs)+(*nSpecies); j++) {
                    dQdbeta[s*((*nDetCovs)+(*nSpecies)) + j] = dQdbeta[s*((*nDetCovs)+(*nSpecies)) + j] + probZ * term[j];
                } // j
            } // k
        } // i
        
    } // s
    
    // regularization
    if (regType[0] == 2) {
        for (int s = 0; s < (*nSpecies); s++) {
            for (int j = 1; j < (*nOccCovs); j++) {
                dQdalpha[s*(*nOccCovs)+j] = dQdalpha[s*(*nOccCovs)+j] - lambdaO[0] * alphas[s*(*nOccCovs)+j];
            } // j
            
            for (int j = 1; j < (*nDetCovs)+(*nSpecies); j++) {
                dQdbeta[s*((*nDetCovs)+(*nSpecies))+j] = dQdbeta[s*((*nDetCovs)+(*nSpecies))+j] - lambdaD[0] * betas[s*((*nDetCovs)+(*nSpecies))+j];
            } // j
        } // s
    }
    
    int count = 0;
    for (int s = 0; s < (*nSpecies); s++) {
        for (int j = 0; j < (*nOccCovs); j++) {
            derivsOfEJLL[count+j] = -dQdalpha[s*(*nOccCovs)+j];
        } // j
        count = count + (*nOccCovs);
        
        for (int j = 0; j < (*nDetCovs); j++) {
            derivsOfEJLL[count+j] = -dQdbeta[s*((*nSpecies)+(*nDetCovs))+j];    
        } // j
        count = count + (*nDetCovs);
        
        for (int ss = 0; ss < (*nSpecies); ss++) {
            if (confusion[s*(*nSpecies)+ss] == 0) { continue; }
            derivsOfEJLL[count] = -dQdbeta[s*((*nSpecies)+(*nDetCovs))+(*nDetCovs)+ss];  
            count = count + 1;
        } // ss
    } // s
    
} // ComputeDerivsOfEJLLAdditive


// Predict Detection
void PredictDetAdditive(int *nSites, int *nVisits, int *nOccCovs, int *nDetCovs, int *nPa, double *alphas, double *betas, int *pi, double *Xo, double *Xd, int *visits, int *Zs, int *nZs, double *predictedDet)
{
    double *beta  = betas;
    double *gamma = betas + (*nDetCovs);
    for (int i = 0; i < (*nSites); i++) {
        double *Xoi = Xo+i*(*nOccCovs);
        
        for (int t = 0; t < visits[i]; t++) {
            double *Xdit = Xd + i*(*nVisits)*(*nDetCovs) + t*(*nDetCovs);

            predictedDet[i*(*nVisits)+t] = 0;
            for (int k = 0; k < (*nZs); k++) {
                int *Z = Zs + k * (*nPa);
                
                // compute occ prob
                double probZOcc = 1;
                for (int s = 0; s < (*nPa); s++) {
                    double *alpha = alphas + s * (*nOccCovs);
                    double probOcc = logistic(innerProduct(Xoi,alpha,*nOccCovs));
                    //                    Rprintf("probOcc is %f.\n",probOcc);
                    probZOcc = probZOcc * pow(probOcc, Z[s]) * pow(1-probOcc, 1-Z[s]);
                } // s
                //                Rprintf("probZocc is %f.\n",probZOcc);
                
                // compute det prob
                double probZDet = logistic(innerProduct(Xdit,beta,*nDetCovs)+innerProduct((double*)Z,gamma,*nPa));
                predictedDet[i*(*nVisits)+t] = predictedDet[i*(*nVisits)+t] + probZOcc*probZDet;
            } // k
        } // t
    } // i
    
} // PredictDet







