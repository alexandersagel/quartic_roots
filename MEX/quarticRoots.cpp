#include "quarticRoots.h"
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Load input variables: alpha_0*x^4+alpha_1*x^3+alpha_2*x^2+alpha_3*x+alpha_4
    double *alpha;
    alpha=mxGetPr(*prhs);
    double a, b,c,d,e;
    a=alpha[0];
    b=alpha[1]/a;
    c=alpha[2]/a;
    d=alpha[3]/a;
    e=alpha[4]/a;
    
    double p=(8.0*c-3.0*pow(b,2))/8.0;
    double q=pow(b,3)/8.0-b*c/2.0+d;
    
    //For special cases of the equation
    bool specialCase=false;
    
    std::complex<double> IMG;
    IMG=-1.0; IMG=sqrt(IMG);
    
    //the root variables. set them to i, so that unassigned values will be ignored
    std::complex<double> x1,x2,x3,x4;
    x1=IMG;
    x2=IMG;
    x3=IMG;
    x4=IMG;
    std::complex<double> Delta_0=pow(c,2)-3.0*b*d+12.0*e;
    std::complex<double> Delta_1=2.0*pow(c,3)-9.0*b*c*d+27.0*pow(b,2)*e+27.0*pow(d,2)-72.0*c*e;
    std::complex<double> Q=pow((Delta_1+sqrt(pow(Delta_1,2)-pow(Delta_0,3)*4.0))/2.0,1.0/3.0);
    std::complex<double> S;
    
    
    if (q==0) {
       //Special case depressed quartic is biquadratic
       specialCase=true;
       double r= (-3.0*pow(b,4)+256.0*e-64.0*b*d+16.0*pow(b,2)*c)/256.0;
       std::complex<double> y1, y2;
       std::complex<double> p2pow2=pow(p/2.0,2);
       y1=-p/2.0+sqrt(p2pow2-r);
       x1=sqrt(y1)-b/4.0;
       if (std::abs(y1)!=0) {
           x2=-sqrt(y1)-b/4.0;
       }
       if (std::abs(p2pow2-r)!=0) {
           y2=-p/2.0-sqrt(p2pow2-r);
           x3=sqrt(y2)-b/4.0;
           if (std::abs(y2)!=0) {
               x4=-sqrt(y2)-b/4.0;
           }
       }
    }
    if (std::abs(Delta_0)==0&&!specialCase) {
        if (std::abs(Delta_1)!=0) {
            //Make sure Q is not 0
            Q=pow(Delta_1,1.0/3.0);
        } else {
            //Special Case: triple root
            // x_1^4+b/2*x_1^2+c/6=0
            // x_2=-b-3*x_1
            specialCase=true;
            x1=-b/4.0+sqrt(pow(b/4,2)-c/6.0);
            if (std::abs((((x1+b)*x1+c)*x1+d)*x1+e)>1e-6) {
                x1=-b/4.0-sqrt(pow(b/4,2)-c/6.0);
            }
            if (x1!=-4.0*b) {
                x2=-b-3.0*x1;
            }
        }
    }
    
    if (!specialCase) {
        S=sqrt(-(2.0*p-(Q+Delta_0/Q))/3.0)/2.0;
        for (unsigned int i=0; std::abs(S)==0; i++) {
            //If S is 0, choose a different root for Q
            Q=std::polar(std::abs(Q),(std::arg(Q)+i*2*3.14159265359)/3.0);
            if (i>2) {
                //Special case: quadruple root: x=b/(4a)
                specialCase=true;
                x1=b/(4.0);
                break;
            }
            S=sqrt(-(2.0*p-(Q+Delta_0/Q))/3.0)/2.0;
        }
    }
    
    //Application of the General formula
    if (!specialCase) {
        x1=-b/4.0-S+sqrt(-4.0*pow(S,2)-2.0*p+q/S)/2.0;
        x2=-b/4.0-S-sqrt(-4.0*pow(S,2)-2.0*p+q/S)/2.0;
        x3=-b/4.0+S+sqrt(-4.0*pow(S,2)-2.0*p-q/S)/2.0;
        x4=-b/4.0+S-sqrt(-4.0*pow(S,2)-2.0*p-q/S)/2.0;
    }
    
    //only return real-valued roots
    unsigned int l=4;
    double tmpRoots[4];
    unsigned int cntr=0;
    if (std::abs(std::imag(x1))<=(1e-6)*std::abs(std::real(x1))) {
        tmpRoots[cntr]=std::real(x1);
        cntr++;
    } else {
        l--;
    }
    if (std::abs(std::imag(x2))<=1e-6*std::abs(std::real(x2))) {
        tmpRoots[cntr]=std::real(x2);
        cntr++;
    } else {
        l--;
    }
    if (std::abs(std::imag(x3))<=1e-6*std::abs(std::real(x3))) {
        tmpRoots[cntr]=std::real(x3);
        cntr++;
    } else {
        l--;
    }
    if (std::abs(std::imag(x4))<=1e-6*std::abs(std::real(x4))) {
        tmpRoots[cntr]=std::real(x4);
        cntr++;
    } else {
        l--;
    }
    
    plhs[0] = mxCreateDoubleMatrix(l,1,mxREAL);
    
    double *roots;
    roots=mxGetPr(*plhs);
    double denom;
    
    for (unsigned int i=0; i<l; i++) {
        //One Newton iteration for numerical stability
        denom=((((4*tmpRoots[i]+3*b)*tmpRoots[i]+2*c)*tmpRoots[i])+d);
        if (denom !=0) {
            tmpRoots[i]=tmpRoots[i]-((((tmpRoots[i]+b)*tmpRoots[i]+c)*tmpRoots[i]+d)*tmpRoots[i]+e)/denom;    
        }
        // }
        roots[i]=tmpRoots[i];
        
    }
    
    
}
