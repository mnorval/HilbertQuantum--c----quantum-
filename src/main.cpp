#include "hilbert.hpp"
#include <cstdio>
int main() {
  auto bell = q::State::zeros(2);
  bell.H(0); bell.CNOT(0,1);
  std::printf("Bell |Phi+>  ||psi||^2=%.6f\n", bell.norm2());
  for (int i=0;i<4;++i)
    std::printf("  |%d%d> p=%.4f amp=(%+.3f,%+.3f)\n", (i>>1)&1, i&1, bell.prob(i),
      bell.amp[i].real(), bell.amp[i].imag());
  std::printf("  <Z0>=%.3f <Z1>=%.3f\n", bell.expectZ(0), bell.expectZ(1));
  auto g = q::State::zeros(2);
  q::grover2(g);
  std::printf("Grover-2 marked |11>\n");
  for (int i=0;i<4;++i)
    std::printf("  |%d%d> p=%.4f\n", (i>>1)&1, i&1, g.prob(i));
}
