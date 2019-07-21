#include "clusterizer.h"

#include <unordered_map>
#include <iostream>

#include "alglib/dataanalysis.h"

Clasterizer::Clasterizer(const std::vector<alglib::real_1d_array> &allSolutions) :
        allSolutions_{allSolutions} {

}


void Clasterizer::Start() {
    std::cout << "Start Clusterization..." << std::endl;

    alglib::real_2d_array solutionsToClustering = FillSolutionArray(allSolutions_);

    alglib::clusterizerstate state;
    alglib::clusterizercreate(state);
    alglib::clusterizersetpoints(state, solutionsToClustering, 2); // 2 -- means Euclidian metrics
    alglib::clusterizerrunahc(state, report_);


    // Create Disjoint Set Union structure
    InitilizeParents();
    for (int i = 0; i < report_.z.rows(); ++i) {
        if (report_.mergedist[i] > max_cluster_size) {
            break;
        }

        int a = find_set(report_.z(i, 0));
        int b = find_set(report_.z(i, 1));

        parents_[b] = a;
    }

    std::unordered_map<int, int> clasters_size;
    for (int cluster : parents_) {
        ++clasters_size[find_set(cluster)];
    }

    for (auto [cluster, size] : clasters_size) {
        std::cout << allSolutions_[cluster].tostring(3) << ", size: " << size << std::endl;
    }
}

alglib::real_2d_array Clasterizer::FillSolutionArray(const std::vector<alglib::real_1d_array>& allSolutions) {
    alglib::real_2d_array solutionsToClustering;
    solutionsToClustering.setlength(allSolutions.size(), allSolutions.at(0).length());

    for (int solutionNum = 0; solutionNum < allSolutions.size(); ++solutionNum) {
        for (int fluxNum = 0; fluxNum < allSolutions[solutionNum].length(); ++fluxNum) {
            solutionsToClustering[solutionNum][fluxNum] = allSolutions[solutionNum][fluxNum];
        }
    }
    return solutionsToClustering;
}

void Clasterizer::InitilizeParents() {
    parents_.resize(allSolutions_.size());
    for (int i = 0; i < allSolutions_.size(); ++i) {
        parents_[i] = i;
    }
}

int Clasterizer::find_set(int v) {
    while (v >= allSolutions_.size()) {
        v = report_.z(v - allSolutions_.size(), 0);
    }

    if (v == parents_[v]) {
        return v;
    } else {
        return parents_[v] = find_set(parents_[v]);
    }
}
