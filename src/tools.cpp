#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  unsigned int est_size = estimations.size();

  if (est_size != ground_truth.size() || est_size == 0) {
    std::cout << "error, unequal size or no data";
    return rmse;
  }

  for (unsigned int i=0; i<est_size; i++) {
    VectorXd diff = estimations[i] - ground_truth[i];
    VectorXd diff_sq = diff.array()*diff.array();
    rmse += diff_sq;
  }

  rmse /= est_size;
  rmse = rmse.array().sqrt();
  return rmse;
}
