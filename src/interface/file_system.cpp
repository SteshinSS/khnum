#include "file_system.h"

#include <string>

FileSystem::FileSystem(const std::string &root) : _Root(root) {}

std::string FileSystem::getModelFile() {
    return _Root + "/model.csv";
}

std::string FileSystem::getMeasuredIsotopesFile() {
    return _Root + "/measured_isotopes.txt";
}

std::string FileSystem::getMeasurementsFile() {
    return _Root + "/measurements.csv";
}

std::string FileSystem::getSubstrateInputFile() {
    return _Root + "/substrate_input.csv";
}

std::string FileSystem::getExcludedMetabolitesFile() {
    return _Root + "/excluded_metabolites.txt";
}