#include <iostream>
#include <functional>
#include <complex>

using namespace std;


complex<double> betaFunction(complex<double> mu, complex<double> a) {
  double beta1 = 9.0/2.0;
  return -1.0/mu*beta1*pow(a, 2);
}


complex<double> RungeKutta(complex<double> t, complex<double> y, const complex<double> &h, const function<complex<double>(complex<double>, complex<double>)> &f, const uint &n) {
  complex<double> k1, k2, k3, k4;
  for (uint i = 1; i <= n; i++) {
    k1 = f(t, y);
    k2 = f(t + h/2.0, y + h*k1/2.0);
    k3 = f(t + h/2.0, y + h*k2/2.0);
    k4 = f(t + h, y + h*k3);
    t = t + h;
    y = y + h/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
  }

  return y;
}

complex<double> runAlphaStrong(complex<double> mu2, complex<double> mu1, complex<double> a1) {
  uint stepCount = 1000;
  complex<double> h = (mu2-mu1)/(double)stepCount;
  return RungeKutta(mu1, a1, h, betaFunction, stepCount);
}

int main() {
  complex<double> a1(0.1, 0.0);
  complex<double> mu1(1.0, 1.0);
  complex<double> mu2(2.15, 3.22);

  std::cout << runAlphaStrong(mu2, mu1, a1) << std::endl;

   return 0;
}
