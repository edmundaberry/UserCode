#ifndef ROOT_BayesMachine
#define ROOT_BayesMachine

using namespace std;

class BayesMachine {
 public:  
  BayesMachine();
  virtual ~BayesMachine();
  //virtual void    Efficiency(int k, int N, double conflevel,
  virtual void    Efficiency(double k, double N, double conflevel,
			     double& mode, double& low, double& high) const;
  void Speak();
 private:
  Double_t        Beta_ab    (double a, double b, double k, double N) const;
  Double_t        Ibetai     (double a, double b, double x) const;
  Double_t        Betai      (double a, double b, double x) const;
  Double_t        SearchLower(double high, double k, double N, double c) const;
  Double_t        SearchUpper(double low, double k, double N, double c) const;
  Double_t        Brent(double ax, double bx, double cx, double tol, double *xmin) const;
  Double_t        Interval(double low) const;
};
#endif
  
