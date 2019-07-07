#include "clusterizer.h"

#include "stdafx.h"
#include "dataanalysis.h"

Clasterizer::Clasterizer(std::vector<alglib::real_1d_array> &allSolutions) {
    alglib::real_2d_array solutionsToClustering;
    solutionsToClustering.setlength(allSolutions.size(), allSolutions.at(0).length());

    std::cout << "Start Clusterization" << std::endl;

    for (int solutionNum = 0; solutionNum < allSolutions.size(); ++solutionNum) {
        std::cout << solutionNum << " ";
        for (int fluxNum = 0; fluxNum < allSolutions[solutionNum].length(); ++fluxNum) {
            solutionsToClustering[solutionNum][fluxNum] = allSolutions[solutionNum][fluxNum];
            std::cout << allSolutions[solutionNum][fluxNum] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    alglib::clusterizerstate state;
    alglib::ahcreport report;
    alglib::clusterizercreate(state);
    alglib::clusterizersetpoints(state, solutionsToClustering, 2); // 2 -- means Euclidian metrics
    alglib::clusterizerrunahc(state, report);

    std::cout << report.z.tostring().c_str();

}