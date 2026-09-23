#pragma once
#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>
namespace q {
using cx = std::complex<double>;
struct State {
  std::size_t nq{0};
  std::vector<cx> amp;
  static State zeros(std::size_t n) {
    State s; s.nq = n; s.amp.assign(1ull << n, cx{0,0}); s.amp[0] = 1; return s;
  }
  double norm2() const { double a=0; for (auto z: amp) a += std::norm(z); return a; }
  void u1(std::size_t qb, cx a, cx b, cx c, cx d) {
    const std::size_t bit = 1ull << qb;
    auto out = amp;
    for (std::size_t i = 0; i < amp.size(); ++i) {
      if (i & bit) continue;
      auto j = i | bit;
      auto x0 = amp[i], x1 = amp[j];
      out[i] = a*x0 + b*x1;
      out[j] = c*x0 + d*x1;
    }
    amp.swap(out);
  }
  void H(std::size_t q) { double s=1.0/std::sqrt(2.0); u1(q,s,s,s,-s); }
  void X(std::size_t q) { u1(q,0,1,1,0); }
  void Z(std::size_t q) { u1(q,1,0,0,-1); }
  void CNOT(std::size_t ctrl, std::size_t tgt) {
    const std::size_t cb = 1ull<<ctrl, tb = 1ull<<tgt;
    auto out = amp;
    for (std::size_t i = 0; i < amp.size(); ++i) {
      if ((i & cb) && !(i & tb)) { auto j = i | tb; out[i]=amp[j]; out[j]=amp[i]; }
    }
    amp.swap(out);
  }
  double prob(std::size_t i) const { return std::norm(amp[i]); }
  double expectZ(std::size_t qb) const {
    double e=0; const std::size_t bit=1ull<<qb;
    for (std::size_t i=0;i<amp.size();++i) { double p=std::norm(amp[i]); e += (i&bit)?-p:p; }
    return e;
  }
};
inline void grover2(State& s) {
  s.H(0); s.H(1);
  s.Z(0); s.Z(1); s.CNOT(0,1); s.Z(1); s.CNOT(0,1); s.Z(0); s.Z(1);
  s.H(0); s.H(1); s.X(0); s.X(1); s.H(1); s.CNOT(0,1); s.H(1); s.X(0); s.X(1); s.H(0); s.H(1);
}
}
