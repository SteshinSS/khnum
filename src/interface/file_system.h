#ifndef CFLEX_FILE_SYSTEM_H
#define CFLEX_FILE_SYSTEM_H

#include <string>

class FileSystem {
public:
    FileSystem (const std::string& root);

    std::string getModelFile();
    std::string getMeasuredIsotopesFile();
    std::string getMeasurementsFile();
    std::string getSubstrateInputFile();
    std::string getExcludedMetabolitesFile();

private:
    const std::string _Root;
};

#endif //CFLEX_FILE_SYSTEM_H
