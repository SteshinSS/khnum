#include "modeller/check_model.h"

#include <vector>
#include <stdexcept>

#include "utilities/measurement.h"

namespace khnum {
namespace modelling_utills {

void CheckMeasurementsMID (std::vector<Measurement> measurements) {
    const double epsilon = 0.0001;
    for (const Measurement& measurement : measurements) {
        double total_mass = 0.0;
        for (double mass : measurement.mid) {
            total_mass += mass;
        }
        bool is_ok = (1.0 - epsilon < total_mass) && (total_mass < 1.0 + epsilon);
        if (!is_ok) {
            throw std::runtime_error("There is measurements with MID is not sum to 1");
        }
    }
}

}
}