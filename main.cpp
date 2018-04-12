#include <iostream>
#include <iomanip>
#include <functional>
#include <complex>
#include <algorithm>
#include <stdexcept>

using namespace std;


class RunAlphaStrong {
public:
  RunAlphaStrong() {
    nf_ = 3;
  }

  complex<double> run(complex<double> mu2, complex<double> mu1, complex<double> a1) {
    function<complex<double>(complex<double>, complex<double>)> rhs = [&](complex<double> mu, complex<double> a) -> complex<double> {
      return -1.0/mu*betaFunction(a);
    };
    complex<double> h = (mu2-mu1)/1000.0;
    complex<double> y = a1, t = mu1;
    int i = 0;
    while(t != mu2) {
      i++;
      // check if overstep
      if (abs(mu2 - (t+h)) > abs(mu2 - t)) {
        h = mu2 - t;
      }
      std::cout << i << std::endl;
      RungeKutta(t, y, h, rhs);
    }
    return y;
  }


private:
  complex<double> betaFunction(complex<double> a) {
    double zeta_3 = 1.2020569031595942;
    double zeta_4 = 1.0823232337111382;
    double zeta_5 = 1.0369277551433699;

    double beta_1 = 11.0/2.0 - 1.0/3.0*nf_;
    double beta_2 = 51.0/4.0 - 19.0/12.0*nf_;
    double beta_3 = 2857.0/64.0 - 5033.0/576.0*nf_ + 325.0/1728.0*pow(nf_, 2);
    double beta_4 = 149753.0/768.0 + 891.0/32.0*zeta_3 - (1078361.0/20736.0 + 1627.0/864.0*zeta_3)*nf_
      + (50065.0/20736.0 + 809.0/1296.0*zeta_3)*pow(nf_, 2) + 1093.0/93312.0*pow(nf_, 3);
    double beta_5 = 2.0/pow(4.0, 5)*(8157455.0/16.0 + 621885.0/2.0*zeta_3 - 88209.0/2.0*zeta_4 - 288090.0*zeta_5
                                     - (336460813.0/1944.0 + 4811164.0/81.0*zeta_3 - 33935.0/6.0*zeta_4 - 1358995.0/27.0*zeta_5)*nf_
                                     + (25960913.0/1944.0 + 698531.0/81.0*zeta_3 - 10526.0/9.0*zeta_4 - 381760.0/81.0*zeta_5)*pow(nf_, 2)
                                     - (630559.0/5832.0 + 48722.0/243.0*zeta_3 - 1618.0/27.0*zeta_4 - 460.0/9.0*zeta_5)*pow(nf_, 3)
                                     + (1205.0/2916.0 - 152.0/81.0*zeta_3)*pow(nf_, 4));
    return beta_1*pow(a, 2) + beta_2*pow(a, 3) + beta_3*pow(a, 4) + beta_4*pow(a, 5) + beta_5*pow(a, 6);
  };

  // Adaptive Explicit Runge-Kutta 4th order (Cash-Karp method)
  void RungeKutta(complex<double> &t, complex<double> &y, complex<double> &h, const function<complex<double>(complex<double>, complex<double>)> &f) {
    double eps = 1e-12;
    complex<double> k1, k2, k3, k4, k5, k6;
    complex<double> yError = 0., yTemp;
    complex<double> hTemp = h;

    // Cash-Karp coefficients
    double c2 = 1.0/5.0, c3 = 3.0/10.0, c4 = 3.0/5.0, c5 = 1.0, c6 = 7.0/8.0;
    double b1 = 37.0/378.0, b2 = 0, b3 = 250.0/621.0, b4 = 125.0/594.0, b5 = 0.0, b6 = 512.0/1771.0;
    double B1 = 2825.0/27648.0, B2 = 0.0, B3 = 18575.0/48384.0, B4 = 13525.0/55296.0, B5 = 277.0/14336.0, B6 = 1.0/4.0;
    double a21 = 1.0/5.0;
    double a31 = 3.0/40.0, a32 = 9.0/40.0;
    double a41 = 3.0/10.0, a42 = -9.0/10.0, a43 = 6.0/5.0;
    double a51 = -11.0/54.0, a52 = 5.0/2.0, a53 = -70.0/27.0, a54 = 35.0/27.0;
    double a61 = 1631.0/55296.0, a62 = 175.0/512.0, a63 = 575.0/13824.0, a64 = 44275.0/110592.0, a65 = 253.0/4096.0;

    while(true) {
      k1 = f(t, y);
      k2 = f(t + c2*hTemp, y + hTemp*(a21*k1));
      k3 = f(t + c3*hTemp, y + hTemp*(a31*k1 + a32*k2));
      k4 = f(t + c4*hTemp, y + hTemp*(a41*k1 + a42*k2 + a43*k3));
      k5 = f(t + c5*hTemp, y + hTemp*(a51*k1 + a52*k2 + a53*k3 + a54*k4));
      k6 = f(t + c6*hTemp, y + hTemp*(a61*k1 + a62*k2 + a63*k4 + a64*k4 + a65*k5));

      yTemp = y + hTemp*(b1*k1 + b2*k2 + b3*k3 + b4*k4 + b5*k5 + b6*k6);
      yError = (b1-B1)*k1 + (b2-B2)*k2 + (b3-B3)*k3 + (b4-B4)*k4 + (b5-B5)*k5 + (b6-B6)*k6;

      double err = abs(yError)/eps;
      if(err < 1) {
        t = t + hTemp;
        y = yTemp;
        if(err>1.89e-4) {
          h = 0.9*hTemp*pow(abs(yError),-0.25);
        } else {
          h = 5.0*hTemp;
        }
        break;
      } else {
        hTemp = 0.1*hTemp;
        std::cout << "hTemp, t \t" << hTemp << "\t" << t << "\t" << yError << std::endl;
        if ( t+hTemp == t ) {
          throw std::runtime_error("RK: stepsize too small");
        }
        continue;
      }
    }
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
  std::cout << std::setprecision(17);
  complex<double> a1(0.101627, 0.0);
  complex<double> mu1(3.1570893124000001, 0.0);
  complex<double> mu2(3.2, 0.0);

  RunAlphaStrong runner;
  std::cout << runner.run(mu2, mu1, a1) << std::endl;

   return 0;
}
