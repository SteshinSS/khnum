#ifndef CFLEX_LINEAR_PROBLEM_H
#define CFLEX_LINEAR_PROBLEM_H

#include "glpk.h"

// glpk wrapper for automatic memory management
// should remake
class LinearProblem {
 public:
    explicit LinearProblem(int size) {
      _linear_problem = glp_create_prob();

      // see glpk Reference Manual pg. 11
      ia = new int [size + 1];
      ja = new int [size + 1];
      ar = new double [size + 1];
  }

  ~LinearProblem() {
      glp_delete_prob(_linear_problem);

      delete [] ia;
      delete [] ja;
      delete [] ar;
  }

  LinearProblem& operator=(const LinearProblem&) = delete;
  LinearProblem(const LinearProblem& item) = delete;
  LinearProblem(LinearProblem&&) = delete;

  operator glp_prob *() {
      return _linear_problem;
  }

  int *ia;
  int *ja;
  double *ar;

 private:
  glp_prob *_linear_problem;
};

#endif //CFLEX_LINEAR_PROBLEM_H
