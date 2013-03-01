#####################################################
#
#  20 October 2009
#  Bryan A. Toth
#  University of California, San Diego
#  btoth@physics.ucsd.edu
# 
#  This script writes the program file for a C++ IPOPT 
#  program defined by the vector field in the file
#  equations.txt and written by the script makecode.py.
#  This file creates an instance of the non-linear file
#  defined by makehpp.py and makecode.py, and interfaces
#  with the IPOPT libraries necessary for solving the
#  optimization problem.
#
#  This script has been developed as part of a suite of 
#  python scripts to define a dynamic parameter estimation
#  problem using the optimization software IPOPT, but is 
#  generally applicable to any application needing
#  discretized derivatives of a vector field.
#
######################################################

import discretize

src = """
// %(probl)s_main.cpp
// Main file for use with IPOPT

#include "IpIpoptApplication.hpp"
#include "%(probl)s_nlp.hpp"
#include <cstdio>

using namespace Ipopt;

int main(int argv, char* argc[])
{

  // Create a new instance of your nlp
  //  (use a SmartPtr, not raw)
  SmartPtr<TNLP> mynlp = new %(probu)s_NLP();

  // Create a new instance of IpoptApplication
  //  (use a SmartPtr, not raw)
  SmartPtr<IpoptApplication> app = new IpoptApplication();

  // Change some options
  // Note: The following choices are only examples, they might not be
  //       suitable for your optimization problem.
  app->Options()->SetNumericValue("tol", 1e-12);
  app->Options()->SetStringValue("mu_strategy", "adaptive");
  app->Options()->SetStringValue("output_file", "ipopt.out");
  // The following overwrites the default name (ipopt.opt) of the
  // options file
  app->Options()->SetStringValue("option_file_name", "%(probl)s.opt");

  // Intialize the IpoptApplication and process the options
  ApplicationReturnStatus status;
  status = app->Initialize();
  if (status != Solve_Succeeded) {
    printf("\\n\\n *** Error during initialization!\\n");
    return (int) status;
  }

  // Ask Ipopt to solve the problem
  status = app->OptimizeTNLP(mynlp);

  if (status == Solve_Succeeded) {
    printf("\\n\\n*** The problem solved!\\n");
  }
  else {
    printf("\\n\\n*** The problem FAILED!\\n");
  }

  // As the SmartPtrs go out of scope, the reference count
  // will be decremented and the objects will automatically
  // be deleted.

  return (int) status;
}
"""

prob = discretize.Problem

open(prob.lower() + '_main.cpp','wt').write(src % dict(probl = prob.lower(),
                                                       probu = prob.upper()))
