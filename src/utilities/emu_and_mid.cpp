#include "utilities/emu_and_mid.h"

#include <algorithm>
#include <tuple>


namespace khnum {
// convolution
Mid operator*(const Mid &lhs, const Mid &rhs) {
    Mid convolve_result(lhs.size() + rhs.size() - 1, 0.0);
    for (size_t mass_shift = 0; mass_shift < convolve_result.size(); ++mass_shift) {
        for (size_t lhs_mass_shift = 0; lhs_mass_shift < lhs.size(); ++lhs_mass_shift) {
            if (mass_shift >= lhs_mass_shift) {
                size_t rhs_mass_shift = mass_shift - lhs_mass_shift;
                if (rhs_mass_shift < rhs.size()) {
                    convolve_result[mass_shift] += lhs[lhs_mass_shift] * rhs[rhs_mass_shift];
                }
            }
        }
    }

    return convolve_result;
}


bool operator==(const EmuAndMid &lhs, const EmuAndMid &rhs) {
    return std::tie(lhs.emu, lhs.mid) == std::tie(rhs.emu, rhs.mid);
}


bool operator<(const EmuAndMid &lhs, const EmuAndMid &rhs) {
    return std::tie(lhs.emu, lhs.mid) < std::tie(rhs.emu, rhs.mid);
}

Mid Normalize(Mid mid) {
    double sum = 0.0;
    for (double part : mid) {
        sum += part;
    }

    for (double& part : mid) {
        part /= sum;
    }
    return mid;
}
} //namespace khnum