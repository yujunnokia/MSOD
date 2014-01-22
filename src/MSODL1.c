#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>

const double EPSILON = 1e-10;

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
    //Rprintf("Sum is %f.\n",sum);
    return(sum);
}

double innerProductInt(int *a, double *b, int n)
{
    double sum = 0;
    for (int i = 0; i < n; i++) {
        sum = sum+a[i]*b[i];
    }
    //Rprintf("Sum is %f.\n",sum);
    return(sum);
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

double multiply(double a, double b)
{
    if (a == 0 && b == -INFINITY) {
        return(0);
    }
    return(a*b);
}

// compute log likelihood of site i
void ComputeiLLNoisyORL1(int *nOccCovs, int *nDetCovs, int *nSpecies, double *alphas, double *betas, double *beta0s, double* confusion, int *Yi, double *Xoi, double *Xdi, int *visit, int *Zs, int *nZs, double* iLL)
{
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
            for (int t = 0; t < *visit; t++) {
                // compute det prob
                double probDet = 1;
                probDet = probDet * (1-beta0s[s]);
                for (int ss = 0; ss < *nSpecies; ss++) {
                    double *beta = betas+s*(*nSpecies)*(*nDetCovs)+ss*(*nDetCovs);
                    probDet = probDet * pow(1-logistic(innerProduct(Xdi+t*(*nDetCovs),beta,*nDetCovs)), Z[ss]*confusion[s*(*nSpecies)+ss]);
                } // ss
                probDet = 1-probDet;
                
                if (probDet == 0) { probDet = EPSILON; }
                if (probDet == 1) { probDet = 1-EPSILON; }
                
                isL = isL * pow(probDet, Yi[t*(*nSpecies)+s]) * pow(1-probDet, 1-Yi[t*(*nSpecies)+s]);
            } // t
        } // s
        iL = iL + isL;
    } // k
    
    *iLL = log(iL);
} // ComputeiLL


// compute expected joint log-likelihood
void ComputeEJLLNoisyORL1(int *nSites, int *nVisits, int *nOccCovs, int *nDetCovs, int *nSpecies, double *alphas, double *betas, double *beta0s, double *probExpectedOccs, double *confusion, int *Y, double *Xo, double *Xd, int *visits, int *Zs, int *nZs, double *EJLL)
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

            if (probZ == 0) { continue; }
            
            double isEJLL = 0;
            for (int s = 0; s < *nSpecies; s++) {
                isEJLL = isEJLL + multiply(Z[s], log(probOccs[s])) + multiply((1-Z[s]), log(1-probOccs[s]));
                    
                for (int t = 0; t < visits[i]; t++) {
                    double *Xdit = Xd + i*(*nVisits)*(*nDetCovs) + t*(*nDetCovs);
                    int Yits = Y[i*(*nVisits)*(*nSpecies) + t*(*nSpecies) + s];
                    
                    // compute det prob
                    double probDet = 1;
                    probDet = probDet * (1-beta0s[s]);
                    for (int ss = 0; ss < *nSpecies; ss++) {
                        double *beta = betas+s*(*nSpecies)*(*nDetCovs)+ss*(*nDetCovs);
                        probDet = probDet * pow(1-logistic(innerProduct(Xdit,beta,*nDetCovs)), Z[ss] * confusion[s*(*nSpecies)+ss]);
                    } // ss
                    probDet = 1 - probDet;

                    if (probDet == 0) { probDet = EPSILON; }
                    if (probDet == 1) { probDet = 1-EPSILON; }
 
                    isEJLL = isEJLL + multiply(Yits, log(probDet)) + multiply((1-Yits), log(1-probDet));
                } // t
            } // s
            iEJLL = iEJLL + multiply(probZ, isEJLL);
        } // k
        *EJLL = *EJLL + iEJLL;
    } // i
} // ComputeEJLL


// compute derivs of EJLL
void ComputeDerivsOfEJLLNoisyORL1(int *nSites, int *nVisits, int *nOccCovs, int *nDetCovs, int *nSpecies, double *alphas, double *betas, double *beta0s, double *probExpectedOccs, double *confusion, int *Y, double *Xo, double *Xd, int *visits, int *Zs, int *nZs, int *fixedLeak, double *derivsOfEJLL)
{
    double dQdalpha[1024];
    double dQdbeta0[1024];
    double dQdbeta[4096];
    double dQdgamma[1024];
    
    for (int i = 0; i < 1024; i++) {
        dQdalpha[i] = 0;
        dQdbeta0[i] = 0;
        dQdgamma[i] = 0;
    }
    for (int i = 0; i < 4096; i++) {
        dQdbeta[i] = 0;
    }

    for (int s = 0; s < (*nSpecies); s++) {
        double *alpha = alphas + s * (*nOccCovs);
        
        // compute derivs of alpha
        for (int i = 0; i < (*nSites); i++) {
            double *Xoi = Xo+i*(*nOccCovs);
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
        
        // compute derivs of beta0
        if (*fixedLeak == 0) {
            for (int i = 0; i < (*nSites); i++) {
                for (int k = 0; k < (*nZs); k++) {
                    int *Z = Zs + k * (*nSpecies);
                    int idx = findIdx(Z,*nSpecies);
                    double probZ = probExpectedOccs[idx*(*nSites) + i];
                    
                    double term = 0;
                    for (int t = 0; t < visits[i]; t++) {
                        double *Xdit = Xd + i*(*nVisits)*(*nDetCovs) + t*(*nDetCovs);
                        
                        // compute prob det
                        double probDet = 1;
                        probDet = probDet * (1-beta0s[s]);
                        for (int sss = 0; sss < *nSpecies; sss++) {
                            double *beta = betas+s*(*nSpecies)*(*nDetCovs)+sss*(*nDetCovs);
                            probDet = probDet * pow(1-logistic(innerProduct(Xdit,beta,*nDetCovs)), Z[sss] * confusion[s*(*nSpecies)+sss]);
                        } // sss
                        probDet = 1 - probDet;
                        
                        if (probDet == 0) { probDet = EPSILON; }
                        if (probDet == 1) { probDet = 1-EPSILON; }
                        
                        int Yits = Y[i*(*nVisits)*(*nSpecies) + t*(*nSpecies) + s];
                        
                        term = term + (Yits / probDet - 1);
                    } // t
                    
                    if (beta0s[s] != 1) {
                        dQdbeta0[s] = dQdbeta0[s] + probZ * term / (1-beta0s[s]);
                    } else {
                        dQdbeta0[s] = dQdbeta0[s] + probZ * term / EPSILON;
                    }
                } // k
            } // i
        }
        
        // compute derivs of beta and gamma
        for (int ss = 0; ss < (*nSpecies); ss++) {
            if (confusion[s*(*nSpecies)+ss] == 0) { continue; }
            for (int i = 0; i < (*nSites); i++) {
                for (int k = 0; k < (*nZs); k++) {
                    int *Z = Zs + k * (*nSpecies);
                    int idx = findIdx(Z,*nSpecies);
                    double probZ = probExpectedOccs[idx*(*nSites) + i];
                                        
                    if (Z[ss] == 0) { continue; }
                    
                    double term[1024];
                    double termGamma = 0;
                    for (int j = 0; j < 1024; j++) { term[j] = 0; }
                    for (int t = 0; t < visits[i]; t++) {
                        double *Xdit = Xd + i*(*nVisits)*(*nDetCovs) + t*(*nDetCovs);
                        
                        // compute prob det
                        double probDet = 1;
                        probDet = probDet * (1-beta0s[s]);
                        for (int sss = 0; sss < *nSpecies; sss++) {
                            double *beta = betas+s*(*nSpecies)*(*nDetCovs)+sss*(*nDetCovs);
                            probDet = probDet * pow(1-logistic(innerProduct(Xdit,beta,*nDetCovs)), Z[sss] * confusion[s*(*nSpecies)+sss]);
                        } // sss
                        probDet = 1 - probDet;
                        
                        if (probDet == 0) { probDet = EPSILON; }
                        if (probDet == 1) { probDet = 1-EPSILON; }
                        
                        double *beta = betas+s*(*nSpecies)*(*nDetCovs)+ss*(*nDetCovs);
                        int Yits = Y[i*(*nVisits)*(*nSpecies) + t*(*nSpecies) + s];
                        double probDetss = logistic(innerProduct(Xdit,beta,*nDetCovs));
                        if (probDetss == 0) { probDetss = EPSILON; }
                        if (probDetss == 1) { probDetss = 1-EPSILON; }
                        
                        // term for beta
                        for (int j = 0; j < (*nDetCovs); j++) {
                            term[j] = term[j] + (Yits/probDet - 1) * probDetss * Xdit[j] * confusion[s*(*nSpecies)+ss];
                        }  
                        
                        // term for gamma
                        termGamma = termGamma + (1 - Yits/probDet) * log(1-probDetss);
                    } // t

                    for (int j = 0; j < (*nDetCovs); j++) {
                        dQdbeta[s*(*nSpecies)*(*nDetCovs) + ss*(*nDetCovs) + j] = dQdbeta[s*(*nSpecies)*(*nDetCovs) + ss*(*nDetCovs) + j] + probZ * term[j];
                    } // j
                    
                    dQdgamma[s*(*nSpecies) + ss] = dQdgamma[s*(*nSpecies) + ss] + probZ * termGamma;
                } // k
            } // i
        } // ss
    
        /*
        // compute derivs of gamma
        for (int ss = 0; ss < (*nSpecies); ss++) {
            if (ss == s) { continue; }
            if (confusion[s*(*nSpecies)+ss] == 0) { continue; }
            for (int i = 0; i < (*nSites); i++) {
                for (int k = 0; k < (*nZs); k++) {
                    int *Z = Zs + k * (*nSpecies);
                    int idx = findIdx(Z,*nSpecies);
                    double probZ = probExpectedOccs[idx*(*nSites) + i];
                                        
                    if (Z[ss] == 0) { continue; }
                    
                    double term = 0;
                    for (int t = 0; t < visits[i]; t++) {
                        double *Xdit = Xd + i*(*nVisits)*(*nDetCovs) + t*(*nDetCovs);
                        
                        // compute prob det
                        double probDet = 1;
                        probDet = probDet * (1-beta0s[s]);
                        for (int sss = 0; sss < *nSpecies; sss++) {
                            double *beta = betas+s*(*nSpecies)*(*nDetCovs)+sss*(*nDetCovs);
                            probDet = probDet * pow(1-logistic(innerProduct(Xdit,beta,*nDetCovs)), Z[sss] * confusion[s*(*nSpecies)+sss]);
                        } // sss
                        probDet = 1 - probDet;
                        
                        if (probDet == 0) { probDet = EPSILON; }
                        if (probDet == 1) { probDet = 1-EPSILON; }
                        
                        double *beta = betas+s*(*nSpecies)*(*nDetCovs)+ss*(*nDetCovs);
                        int Yits = Y[i*(*nVisits)*(*nSpecies) + t*(*nSpecies) + s];
                        double probDetss = logistic(innerProduct(Xdit,beta,*nDetCovs));
                        
                        if (probDetss == 0) { probDetss = EPSILON; }
                        if (probDetss == 1) { probDetss = 1-EPSILON; }
                        term = term + (1 - 1/probDet * Yits) * log(1-probDetss);
                    } // t
                    
                    dQdgamma[s*(*nSpecies) + ss] = dQdgamma[s*(*nSpecies) + ss] + probZ * term;
                } // k
            } // i            
        } // ss
        */
    } // s
    
    // put derivs into an array
    int count = 0;
    for (int s = 0; s < (*nSpecies); s++) {
        for (int ss = 0; ss < (*nSpecies); ss++) {
            if (ss == s) { continue; } 
            derivsOfEJLL[count] = dQdgamma[s*(*nSpecies) + ss];
            count = count + 1;
        }
    }
    
    if (*fixedLeak == 0) {
        for (int s = 0; s < (*nSpecies); s++) {
            derivsOfEJLL[count] = dQdbeta0[s];
            count = count + 1;
        }
    }
    
    for (int s = 0; s < (*nSpecies); s++) {
        for (int j = 0; j < (*nOccCovs); j++) {
            derivsOfEJLL[count+j] = dQdalpha[s*(*nOccCovs)+j];
        }
        count = count + (*nOccCovs);
        
        for (int ss = 0; ss < (*nSpecies); ss++) {
            for (int j = 0; j < (*nDetCovs); j++) {
                derivsOfEJLL[count+j] = dQdbeta[s*(*nSpecies)*(*nDetCovs)+ss*(*nDetCovs)+j];    
            }
            count = count + (*nDetCovs);
        } // ss
    } // s
    
    /*
    for (int j = 0; j < count; j++) {
        Rprintf("%f,",derivsOfEJLL[j]);
    }
    Rprintf("\n");
     */
} // ComputeDerivsOfEJLL


// Predict Detection
void PredictDetNoisyORL1(int *nSites, int *nVisits, int *nOccCovs, int *nDetCovs, int *nSpecies, double *alphas, double *betas, double *beta0s, double* confusion, double *Xo, double *Xd, int *visits, int *Zs, int *nZs, double *predictedDet)
{
    for (int i = 0; i < (*nSites); i++) {
        double *Xoi = Xo+i*(*nOccCovs);
        
        for (int t = 0; t < visits[i]; t++) {
            double *Xdit = Xd + i*(*nVisits)*(*nDetCovs) + t*(*nDetCovs);
                        
            predictedDet[i*(*nVisits)+t] = 0;
            for (int k = 0; k < (*nZs); k++) {
                int *Z = Zs + k * (*nSpecies);
                
                // compute occ prob
                double probZOcc = 1;
                for (int s = 0; s < (*nSpecies); s++) {
                    double *alpha = alphas + s * (*nOccCovs);
                    double probOcc = logistic(innerProduct(Xoi,alpha,*nOccCovs));
                    probZOcc = probZOcc * pow(probOcc, Z[s]) * pow(1-probOcc, 1-Z[s]);
                } // s
                
                // compute det prob
                double probZDet = 1;
                probZDet = probZDet * (1-beta0s[0]);
                for (int s = 0; s < (*nSpecies); s++) {
                    double *beta = betas+s*(*nDetCovs);
                    probZDet = probZDet * pow(1-logistic(innerProduct(Xdit,beta,*nDetCovs)), Z[s]*confusion[s]);
                } // s
                probZDet = 1 - probZDet;
                
                if (probZDet == 0) { probZDet = EPSILON; }
                if (probZDet == 1) { probZDet = 1-EPSILON; }
                
                predictedDet[i*(*nVisits)+t] = predictedDet[i*(*nVisits)+t] + probZOcc*probZDet;
            } // k
        } // t
    } // i

} // PredictDet


