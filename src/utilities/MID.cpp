#include "MID.h"

#include <algorithm>

MID operator*(const MID &lhs, const MID &rhs) {
    MID convolve_result(lhs.size() + rhs.size() - 1, 0.0);
    for (int mass_shift = 0; mass_shift < convolve_result.size(); ++mass_shift) {
        for (int lhs_mass_shift = 0; lhs_mass_shift < lhs.size(); ++lhs_mass_shift) {
            int rhs_mass_shift = mass_shift - lhs_mass_shift;
            if (rhs_mass_shift >= 0 || rhs_mass_shift < rhs.size()) {
                convolve_result[mass_shift] += lhs[lhs_mass_shift] * rhs[rhs_mass_shift];
            }
        }
    }

    return convolve_result;
}