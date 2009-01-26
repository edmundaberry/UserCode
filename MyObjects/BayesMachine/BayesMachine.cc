#include "Riostream.h"
#include "TROOT.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include "TArrow.h"
#include "TBox.h"
#include "TVirtualPad.h"
#include "TF1.h"
#include "TH1.h"
#include "TVector.h"
#include "TVectorD.h"
#include "TClass.h"
#include "BayesMachine.h"

ClassImp(BayesMachine)

/* 

   Method and code taken directly from
   Marc Paterno (paterno@fnal.gov)
   FNAL/CD
   
   This is a class for calculating Bayesian errors on an efficiency calculation, 
   given the number of events that pass (nPass),
   out of the total number of events (nTotal),
   and the confidence level (confLevel, generally 0.683).
   
   Usage is: 
   BayesMachine *b = new BayesMachine();
   
   double effi, errLow, errHigh;
   double confidenceLevel

   b -> Efficiency(nPass,nTotal,confLevel, &effi, &errLow, &errHigh);

   cout << "efficiency is " << effi 
   << "(+ " << errHigh 
   << "), (-" << errLow << ")" << endl;

 */

namespace {
   unsigned long GLOBAL_k;   // used to pass k[i] into equations
   unsigned long GLOBAL_N;   // used to pass N[i] into equations
   double        CONFLEVEL;  // confidence level for the interval
}

BayesMachine::BayesMachine(){}
BayesMachine::~BayesMachine(){}

void BayesMachine::Speak(){
   cout << "Hello World" << endl;
}

double BayesMachine::Ibetai(double a, double b, double x) const
{
   // Calculates the incomplete beta function  I_x(a,b); this is
   // the incomplete beta function divided by the complete beta function

   double bt;
   if (x < 0.0 || x > 1.0) {
      Error("Ibetai","Illegal x in routine Ibetai: x = %g",x);
      return 0;
   }
   if (x == 0.0 || x == 1.0)
      bt=0.0;
   else
      bt=TMath::Exp(TMath::LnGamma(a+b)-TMath::LnGamma(a)-TMath::LnGamma(b)+a*log(x)+b*log(1.0-x));

   if (x < (a+1.0)/(a+b+2.0))
      return bt*TMath::BetaCf(x,a,b)/a;
   else
      return 1.0-bt*TMath::BetaCf(1-x,b,a)/b;
}

double BayesMachine::Interval(double low) const
{
   // Return the length of the interval starting at low
   // that contains CONFLEVEL of the x^GLOBAL_k*(1-x)^(GLOBAL_N-GLOBAL_k)
   // distribution.
   // If there is no sufficient interval starting at low, we return 2.0

   double high = SearchUpper(low, GLOBAL_k, GLOBAL_N, CONFLEVEL);
   if (high == -1.0) return 2.0; //  so that this won't be the shortest interval
   return (high - low);
}

//double BayesMachine::Beta_ab(double a, double b, int k, int N) const
double BayesMachine::Beta_ab(double a, double b, double k, double N) const
{
   // Calculates the fraction of the area under the
   // curve x^k*(1-x)^(N-k) between x=a and x=b

   if (a == b) return 0;    // don't bother integrating over zero range
   //int c1 = k+1;
   //int c2 = N-k+1;

   double c1 = k+1;
   double c2 = N-k+1;
   return Ibetai(c1,c2,b)-Ibetai(c1,c2,a);
}

//double BayesMachine::SearchUpper(double low, int k, int N, double c) const
double BayesMachine::SearchUpper(double low, double k, double N, double c) const
{
   // Integrates the binomial distribution with
   // parameters k,N, and determines what is the upper edge of the
   // integration region which starts at low which contains probability
   // content c. If an upper limit is found, the value is returned. If no
   // solution is found, -1 is returned.
   // check to see if there is any solution by verifying that the integral up
   // to the maximum upper limit (1) is greater than c

   double integral = Beta_ab(low, 1.0, k, N);
   if (integral == c) return 1.0;    // lucky -- this is the solution
   if (integral < c) return -1.0;    // no solution exists
   double too_high = 1.0;            // upper edge estimate
   double too_low = low;
   double test;

   // use a bracket-and-bisect search
   // LM: looping 20 times might be not enough to get an accurate precision.
   // see for example bug https://savannah.cern.ch/bugs/?30246
   // now break loop when difference is less than 1E-15
   // t.b.d: use directly the beta distribution quantile

   for (int loop=0; loop<50; loop++) {
      test = 0.5*(too_low + too_high);
      integral = Beta_ab(low, test, k, N);
      if (integral > c)  too_high = test;
      else too_low = test;
      if ( TMath::Abs(integral - c) <= 1.E-15) break;
   }
   return test;
}

//double BayesMachine::SearchLower(double high, int k, int N, double c) const
double BayesMachine::SearchLower(double high, double k, double N, double c) const
{
   // Integrates the binomial distribution with
   // parameters k,N, and determines what is the lower edge of the
   // integration region which ends at high, and which contains
   // probability content c. If a lower limit is found, the value is
   // returned. If no solution is found, the -1 is returned.
   // check to see if there is any solution by verifying that the integral down
   // to the minimum lower limit (0) is greater than c

   double integral = Beta_ab(0.0, high, k, N);
   if (integral == c) return 0.0;      // lucky -- this is the solution
   if (integral < c) return -1.0;      // no solution exists
   double too_low = 0.0;               // lower edge estimate
   double too_high = high;
   double test;

   // use a bracket-and-bisect search
   // LM: looping 20 times might be not enough to get an accurate precision.
   // see for example bug https://savannah.cern.ch/bugs/?30246
   // now break loop when difference is less than 1E-15
   // t.b.d: use directly the beta distribution quantile

   for (int loop=0; loop<50; loop++) {
      test = 0.5*(too_high + too_low);
      integral = Beta_ab(test, high, k, N);
      if (integral > c)  too_low = test;
      else too_high = test;
      if ( TMath::Abs(integral - c) <= 1.E-15) break;
   }
   return test;
}

double BayesMachine::Brent(double ax, double bx, double cx, double tol, double *xmin) const
{
   // Implementation file for the numerical equation solver library.
   // This includes root finding and minimum finding algorithms.
   // Adapted from Numerical Recipes in C, 2nd edition.
   // Translated to C++ by Marc Paterno

   const int    kITMAX = 100;
   const double kCGOLD = 0.3819660;
   const double kZEPS  = 1.0e-10;

   int iter;
   double a,b,d=0.,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
   double e=0.0;

   a=(ax < cx ? ax : cx);
   b=(ax > cx ? ax : cx);
   x=w=v=bx;
   fw=fv=fx=Interval(x);
   for (iter=1;iter<=kITMAX;iter++) {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*TMath::Abs(x)+kZEPS);
      if (TMath::Abs(x-xm) <= (tol2-0.5*(b-a))) {
         *xmin=x;
         return fx;
      }
      if (TMath::Abs(e) > tol1) {
         r=(x-w)*(fx-fv);
         q=(x-v)*(fx-fw);
         p=(x-v)*q-(x-w)*r;
         q=2.0*(q-r);
         if (q > 0.0) p = -p;
         q=TMath::Abs(q);
         etemp=e;
         e=d;
         if (TMath::Abs(p) >= TMath::Abs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)) d=kCGOLD*(e=(x >= xm ? a-x : b-x));
         else {
            d=p/q;
            u=x+d;
            if (u-a < tol2 || b-u < tol2) d=TMath::Sign(tol1,xm-x);
         }
      } else {
         d=kCGOLD*(e=(x >= xm ? a-x : b-x));
      }
      u=(TMath::Abs(d) >= tol1 ? x+d : x+TMath::Sign(tol1,d));
      fu=Interval(u);
      if (fu <= fx) {
         if (u >= x) a=x; else b=x;
         v  = w;
         w  = x;
         x  = u;
         fv = fw;
         fw = fx;
         fx = fu;
      } else {
         if (u < x) a=u; else b=u;
         if (fu <= fw || w == x) {
            v=w;
            w=u;
            fv=fw;
            fw=fu;
         } else if (fu <= fv || v == x || v == w) {
            v=u;
            fv=fu;
         }
      }
   }
   Error("Brent","Too many interations");
   *xmin=x;
   return fx;
}


//void BayesMachine::Efficiency(int k, int N, double conflevel,
void BayesMachine::Efficiency(double k, double N, double conflevel,
			      double& mode, double& low, double& high) const {
  // Calculate the shortest central confidence interval containing the required
   // probability content.
   // Interval(low) returns the length of the interval starting at low
   // that contains CONFLEVEL probability. We use Brent's method,
   // except in two special cases: when k=0, or when k=N
   // Main driver routine
   // Author: Marc Paterno

   //If there are no entries, then we know nothing, thus return the prior...
   if (0==N) {
      mode = .5; low = 0.0; high = 1.0;
      return;
   }
 
   // Calculate the most probable value for the posterior cross section.
   // This is easy, 'cause it is just k/N
   double efficiency = (double)k/N;

   double low_edge;
   double high_edge;

   if (k == 0) {
      low_edge = 0.0;
      high_edge = SearchUpper(low_edge, k, N, conflevel);
   } else if (k == N) {
      high_edge = 1.0;
      low_edge = SearchLower(high_edge, k, N, conflevel);
   } else {
      GLOBAL_k = k;
      GLOBAL_N = N;
      CONFLEVEL = conflevel;
      Brent(0.0, 0.5, 1.0, 1.0e-9, &low_edge);
      high_edge = low_edge + Interval(low_edge);
   }

   // return output
   mode = efficiency;
   low  = efficiency - low_edge;
   high = high_edge - efficiency;
} 

