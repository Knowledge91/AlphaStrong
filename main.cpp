#include <iostream>
#include <functional>
#include <complex>

using namespace std;


class RunAlphaStrong {
public:
  RunAlphaStrong() {
    nf_ = 3;
  }

  complex<double> run(complex<double> mu2, complex<double> mu1, complex<double> a1) {
    uint stepCount = 1000;
    complex<double> h = (mu2-mu1)/(double)stepCount;
    function<complex<double>(complex<double>, complex<double>)> rhs = [&](complex<double> mu, complex<double> a) -> complex<double> {
      return -1.0/mu*betaFunction(a);
    };
    return RungeKutta(mu1, a1, h, rhs, stepCount);
  }


private:
  complex<double> betaFunction(complex<double> a) {
    double zeta_3 = 1.2020569031595942;

    double beta_1 = 11.0/2.0 - 1.0/3.0*nf_;
    double beta_2 = 51.0/4.0 - 19.0/12.0*nf_;
    double beta_3 = 2857.0/64.0 - 5033.0/576.0*nf_ + 325.0/1728.0*pow(nf_, 2);
    double beta_4 = 149753.0/768.0 + 891.0/32.0*zeta_3 - (1078361.0/20736.0 + 1627.0/864.0*zeta_3)*nf_
      + (50065.0/20736.0 + 809.0/1296.0*zeta_3)*pow(nf_, 2) + 1093.0/93312.0*pow(nf_, 3);
    return beta_1*pow(a, 2) + beta_2*pow(a, 3) + beta_3*pow(a, 4) + beta_4*pow(a, 5);
  };

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

  uint nf_;
};

// complex<double> betaFunction(complex<double> mu, complex<double> a) {
//   double beta1 = 9.0/2.0;
//   return -1.0/mu*beta1*pow(a, 2);
// }



// complex<double> runAlphaStrong(complex<double> mu2, complex<double> mu1, complex<double> a1) {
//   uint stepCount = 1000;
//   complex<double> h = (mu2-mu1)/(double)stepCount;
//   return RungeKutta(mu1, a1, h, betaFunction, stepCount);
// }

int main() {
  complex<double> a1(0.101627, 0.0);
  complex<double> mu1(3.1570893124000001, 0.0);
  complex<double> mu2(3.0, 0.0);

  RunAlphaStrong runner;
  std::cout << runner.run(mu2, mu1, a1) << std::endl;

   return 0;
}
