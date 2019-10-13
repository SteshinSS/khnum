#pragma once

#include "alglib/ap.h"
#include "alglib/dataanalysis.h"

#include <vector>


namespace khnum {
class Clasterizer {
public:
    Clasterizer(const std::vector<alglib::real_1d_array> &allSolution);
    void Start();

private:

    alglib::real_2d_array FillSolutionArray(const std::vector<alglib::real_1d_array> &allSolutions);

    int FindSet(int v);

    void InitilizeParents();

private:
    const double max_cluster_size = 1.0;

    const std::vector<alglib::real_1d_array> &allSolutions_;
    std::vector<int> parents_;
    alglib::ahcreport report_;
};
} // namespace khnum