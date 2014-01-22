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

double F(double x) 
{
    return log(1-exp(-x));
}

// update variational parameters
void  UpdateVariationalParams(int *nSites, 
                              int *nVisits, 
                              int *nOccCovs, 
                              int *nDetCovs, 
                              int *nSpecies, 
                              double *alphas, 
                              double *betas, 
                              double *beta0s, 
                              double *Z, 
                              int *Y, 
                              double *Xo, 
                              double *Xd, 
                              int *visits, 
                              double *q)
{
    double theta0[100];    
    for (int s = 0; s < (*nSpecies); s++) {
        theta0[s] = - log(1-beta0s[s]);
        if (theta0[s] == 0) { theta0[s] = EPSILON; }
    }

    for (int i = 0; i < (*nSites); i++) {
        for (int t = 0; t < visits[i]; t++) {
            double *Xdit = Xd + i*(*nVisits)*(*nDetCovs) + t*(*nDetCovs);

            for (int s = 0; s < (*nSpecies); s++) {
                double denom = 0.0;
                for (int r = 0; r < (*nSpecies); r++) {
                    double Zir = *(Z+i*(*nSpecies)+r);
                    double *beta = betas+s*(*nSpecies)*(*nDetCovs)+r*(*nDetCovs);
                    double ditrs = logistic(innerProduct(Xdit,beta,*nDetCovs));
                    if (ditrs == 1) { ditrs = 1-EPSILON; }
                    if (ditrs == 0) { ditrs = EPSILON; }
                    double theta = - log(1-ditrs);
                    
                    //Rprintf("i:%d, t:%d, s:%d, r:%d, d:%f \n",(i+1),(t+1),(s+1),(r+1),ditrs);
                    
                    double qitrs = 0.5;
                    int k = 0; 
                    while (k < 50) {
                        double A = exp(-theta0[s]-theta/qitrs);
                        if (A == 0) { A = EPSILON; } 
                        if (A == 1) { A = 1-EPSILON; } 
                        double F0 = F(theta0[s]);
                        
                        //if (i == 0) { Rprintf("q is %f\n",qitrs); }
                        
                        qitrs = - Zir / F0 * (qitrs*(F(theta0[s]+theta/qitrs) - F0)-A/(1-A)*theta);
                        
                        if (qitrs <= 0) {
                            qitrs = EPSILON;
                            break;
                        }
                        
                        if (isnan(qitrs)) { 
                            Rprintf("it is NAN...\n"); 
                            Rprintf("F0 is %f, theta is %f, A is %f. \n",F0,theta,A);
                        }
                                                
                        k++;
                    } // k
                                                            
                    q[i*(*nVisits)*(*nSpecies)*(*nSpecies)+t*(*nSpecies)*(*nSpecies)+r*(*nSpecies)+s] = qitrs;
                    denom += qitrs;
                    
                    //Rprintf("i:%d, t:%d, s:%d, r:%d, q:%f \n",(i+1),(t+1),(s+1),(r+1),qitrs);
                } // r   
                
                // normalize w.r.t. r
                if (denom == 0) { denom = EPSILON; }
                for (int r = 0; r < (*nSpecies); r++) {
                    q[i*(*nVisits)*(*nSpecies)*(*nSpecies)+t*(*nSpecies)*(*nSpecies)+r*(*nSpecies)+s] /= denom;
                } // r
            } // s
        } // t
    } // i 
}

// 

// compute expected joint log-likelihood
void ComputeEJLLMSODVL(int *nSites, 
                       int *nVisits, 
                       int *nOccCovs, 
                       int *nDetCovs, 
                       int *nSpecies, 
                       double *alphas, 
                       double *betas, 
                       double *beta0s, 
                       double *q, 
                       double *Z, 
                       int *Y, 
                       double *Xo, 
                       double *Xd, 
                       int *visits, 
                       double *EJLL)
{
    double theta0[1024];
    for (int s = 0; s < (*nSpecies); s++) {
        theta0[s] = - log(1-beta0s[s]);
    }

    *EJLL = 0;
    for (int i = 0; i < *nSites; i++) {
        double *Xoi = Xo+i*(*nOccCovs);
        
        // compute occ prob for all species
        double O[1024];
        for (int s = 0; s < (*nSpecies); s++) {
            double *alpha = alphas + s * (*nOccCovs);
            O[s] = logistic(innerProduct(Xoi,alpha,*nOccCovs));
            if (O[s] == 0) { O[s] = EPSILON; }
            if (O[s] == 1) { O[s] = 1-EPSILON; }            
        } // s
        
        // compute EJLL for site i
        double iEJLL = 0;

        // term on occupancy
        for (int r = 0; r < (*nSpecies); r++) {
            double Zir = *(Z+i*(*nSpecies)+r);
            iEJLL += Zir * log(O[r]/(1-O[r])) +log(1-O[r]);
        }
        
        // term on detection
        for (int t = 0; t < visits[i]; t++) {
            double *Xdit = Xd + i*(*nVisits)*(*nDetCovs) + t*(*nDetCovs);
            for (int s = 0; s < (*nSpecies); s++) {
                int Yits = *(Y + i*(*nVisits)*(*nSpecies) + t*(*nSpecies) + s);
                
                iEJLL += (-theta0[s]) * (1-Yits);
                for (int r = 0; r < (*nSpecies); r++) {
                    double Zir = *(Z+i*(*nSpecies)+r); 
                    double qitrs = *(q+i*(*nVisits)*(*nSpecies)*(*nSpecies)+t*(*nSpecies)*(*nSpecies)+r*(*nSpecies)+s);
                    double *beta = betas+s*(*nSpecies)*(*nDetCovs)+r*(*nDetCovs);
                    double ditrs = logistic(innerProduct(Xdit,beta,*nDetCovs));
                    if (ditrs == 1) { ditrs = 1-EPSILON; }
                    double theta = - log(1-ditrs);
                                        
                    iEJLL += (-Zir*theta)*(1-Yits);
                    iEJLL += Zir * Yits * qitrs * (F(theta0[s] + theta/qitrs) - F(theta0[s])) + Yits * qitrs * F(theta0[s]);
                } // r
            } // s
        } // t
        
        *EJLL = *EJLL + iEJLL;
    } // i
} // ComputeEJLLMSODVL

// compute derivs of EJLL
void ComputeDerivsOfEJLLMSODVL(int *nSites, 
                               int *nVisits, 
                               int *nOccCovs, 
                               int *nDetCovs, 
                               int *nSpecies, 
                               double *alphas, 
                               double *betas, 
                               double *beta0s, 
                               double *q, 
                               double *Z, 
                               int *Y, 
                               double *Xo, 
                               double *Xd, 
                               int *visits, 
                               double *derivsOfEJLL)
{
    double dQdalpha[1024];
    for (int i = 0; i < 1024; i++) {
        dQdalpha[i] = 0;
    }
    
    double dQdbeta[4096];
    for (int i = 0; i < 4096; i++) {
        dQdbeta[i] = 0;
    }
    
    double theta0[1024];
    for (int s = 0; s < (*nSpecies); s++) {
        theta0[s] = - log(1-beta0s[s]);
    }
    
    for (int s = 0; s < (*nSpecies); s++) {
        double *alpha = alphas + s*(*nOccCovs);
        
        // compute derivs of alpha
        for (int i = 0; i < (*nSites); i++) {
            double *Xoi = Xo+i*(*nOccCovs);
            double Ois = logistic(innerProduct(Xoi,alpha,*nOccCovs));
            double Zis = *(Z+i*(*nSpecies)+s);
            
            for (int j = 0; j < (*nOccCovs); j++) {
                dQdalpha[s*(*nOccCovs)+j] += (Zis-Ois)*Xoi[j];
            } // j
        } // i
        
        // compute derivs of beta 
        for (int r = 0; r < (*nSpecies); r++) {
            for (int i = 0; i < (*nSites); i++) {
                double Zir = *(Z+i*(*nSpecies)+r); 
            
                for (int t = 0; t < visits[i]; t++) {
                    double *Xdit = Xd+i*(*nVisits)*(*nDetCovs)+t*(*nDetCovs);
                    double *beta = betas+s*(*nSpecies)*(*nDetCovs)+r*(*nDetCovs);
                    int Yits = *(Y + i*(*nVisits)*(*nSpecies) + t*(*nSpecies) + s);
                    double ditrs = logistic(innerProduct(Xdit,beta,*nDetCovs));
                    if (ditrs == 0) { ditrs = EPSILON; }
                    if (ditrs == 1) { ditrs = 1-EPSILON; }
                    double theta = - log(1-ditrs);
                    double qitrs = *(q+i*(*nVisits)*(*nSpecies)*(*nSpecies)+t*(*nSpecies)*(*nSpecies)+r*(*nSpecies)+s);
                    double A = exp(- theta0[s] - theta/qitrs);
                    
                    for (int j = 0; j < (*nDetCovs); j++) {
                        dQdbeta[s*(*nSpecies)*(*nDetCovs)+r*(*nDetCovs)+j] += Zir*(Yits/(1-A)-1)*ditrs*Xdit[j];
                    } // j
                } // t
            } // i
        } // ss
    } // s
    
    // put derivs into an array
    int count = 0;
    for (int s = 0; s < (*nSpecies); s++) {
        for (int j = 0; j < (*nOccCovs); j++) {
            derivsOfEJLL[count+j] = dQdalpha[s*(*nOccCovs)+j];
        } // j
        count = count + (*nOccCovs);
        
        for (int r = 0; r < (*nSpecies); r++) {
            for (int j = 0; j < (*nDetCovs); j++) {
                derivsOfEJLL[count+j] = dQdbeta[s*(*nSpecies)*(*nDetCovs)+r*(*nDetCovs)+j];    
            }
            count = count + (*nDetCovs);
        } // s
    } // r
} // ComputeDerivsOfEJLL


// compute expected occupancy
void ComputeExpectedOccMSODVL(int *nSites, 
                               int *nVisits, 
                               int *nOccCovs, 
                               int *nDetCovs, 
                               int *nSpecies, 
                               double *alphas, 
                               double *betas, 
                               double *beta0s, 
                               double *q, 
                               int *Y, 
                               double *Xo, 
                               double *Xd, 
                               int *visits, 
                               double *Z)
{
    double theta0[100];    
    for (int s = 0; s < (*nSpecies); s++) {
        theta0[s] = - log(1-beta0s[s]);
    }

    for (int i = 0; i < (*nSites); i++) {
        double *Xoi = Xo+i*(*nOccCovs);
        
        for (int r = 0; r < (*nSpecies); r++) {
            double *alpha = alphas + r*(*nOccCovs);
            double Oir = logistic(innerProduct(Xoi,alpha,*nOccCovs));
            
            double Zir1 = Oir;
            double Zir0 = 1-Oir;
            for (int t = 0; t < visits[i]; t++) {
                for (int s = 0; s < (*nSpecies); s++) {
                    double *Xdit = Xd+i*(*nVisits)*(*nDetCovs)+t*(*nDetCovs);
                    double *beta = betas+s*(*nSpecies)*(*nDetCovs)+r*(*nDetCovs);
                    int Yits = *(Y + i*(*nVisits)*(*nSpecies) + t*(*nSpecies) + s);
                    double ditrs = logistic(innerProduct(Xdit,beta,*nDetCovs));
                    if (ditrs == 0) { ditrs = EPSILON; }
                    if (ditrs == 1) { ditrs = 1-EPSILON; }
                    double theta = - log(1-ditrs);
                    double qitrs = *(q+i*(*nVisits)*(*nSpecies)*(*nSpecies)+t*(*nSpecies)*(*nSpecies)+r*(*nSpecies)+s);
                    
                    double probDet1 = 0.0;
                    double probDet0 = 0.0;
                    if (Yits == 1) {
                        probDet1 = exp(qitrs*(F(theta0[s]+theta/qitrs) - F(theta0[s])) + qitrs*F(theta0[s]));
                        probDet0 = exp(qitrs*F(theta0[s]));
                    } else {
                        probDet1 = exp(-theta);
                        probDet0 = 1.0;
                    }
                    
                    Zir1 *= probDet1;
                    Zir0 *= probDet0;
                } // s
            } // t
            
            Z[i*(*nSpecies)+r] = Zir1/(Zir1+Zir0);
        }  // r
    } // i 
}

// Predict Detection
void PredictDetNoisyOR(int *nSites, 
                       int *nVisits, 
                       int *nOccCovs, 
                       int *nDetCovs, 
                       int *nSpecies, 
                       double *alphas, 
                       double *betas, 
                       double *beta0s, 
                       int* confusion, 
                       double *Xo, 
                       double *Xd, 
                       int *visits, 
                       int *Zs, 
                       int *nZs, 
                       double *predictedDet)
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

