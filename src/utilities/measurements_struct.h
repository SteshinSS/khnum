#ifndef CFLEX_MEASUREMENTS_STRUCT_H
#define CFLEX_MEASUREMENTS_STRUCT_H

using MeasurementValue = double;
using MeasurementError = double;

struct Measurement {
  MeasurementValue value;
  MeasurementError error;
};

#endif //CFLEX_MEASUREMENTS_STRUCT_H
