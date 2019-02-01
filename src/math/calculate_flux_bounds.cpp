#include "calculate_flux_bounds.h"

#include "reaction_struct.h"

#include <vector>
#include <cmath>


std::vector<Reaction> CalculateFluxBounds(std::vector<Reaction> reactions) {
    for (Reaction &reaction : reactions) {
        if (std::isnan(reaction.basis)) {
            reaction.computed_upper_bound = reaction.setted_upper_bound;
            reaction.computed_lower_bound = reaction.setted_lower_bound;
        } else {
            if (std::isnan(reaction.deviation)) {
                reaction.computed_upper_bound = reaction.basis;
                reaction.computed_lower_bound = reaction.basis;
            } else {
                reaction.computed_upper_bound = reaction.basis + reaction.deviation;
                reaction.computed_lower_bound = reaction.basis - reaction.deviation;
            }
        }
    }

    return reactions;
}