#####################################################
#
#  20 October 2009
#  Bryan A. Toth
#  University of California, San Diego
#  btoth@physics.ucsd.edu
# 
#  This script writes the class header file for a C++ IPOPT 
#  program defined by the vector field in the file
#  equations.txt and written by the script makecode.py.
#  The file written by this script defines a class that is 
#  further described in the file written by makecode.py.
#
#  This script has been developed as part of a suite of 
#  python scripts to define a dynamic parameter estimation
#  problem using the optimization software IPOPT, but is 
#  generally applicable to any application needing
#  discretized derivatives of a vector field.
#
######################################################

import discretize
prob = discretize.Problem

template = """
// %(probl)s.hpp
// Header file for %(probl)s.cpp

#ifndef __%(probu)s_NLP_HPP__
#define __%(probu)s_NLP_HPP__

#include "IpTNLP.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstring>

using namespace std;
using namespace Ipopt;

class %(probu)s_NLP : public TNLP {
public:
  /** default constructor */
  %(probu)s_NLP();

  /** default destructor */
  virtual ~%(probu)s_NLP();

  /**@name Overloaded from TNLP */
  //@{
  /** Method to return some info about the nlp */
  virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag, IndexStyleEnum& index_style);

  /** Method to return the bounds for my problem */
  virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                               Index m, Number* g_l, Number* g_u);

  /** Method to return the starting point for the algorithm */
  virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                  bool init_z, Number* z_L, Number* z_U,
                                  Index m, bool init_lambda,
                                  Number* lambda);

  /** Method to return the objective value */
  virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

  /** Method to return the gradient of the objective */
  virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

  /** Method to return the constraint residuals */
  virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

  /** Method to return:
   *   1) The structure of the jacobian (if "values" is NULL)
   *   2) The values of the jacobian (if "values" is not NULL)
   */
  virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                          Index m, Index nele_jac, Index* iRow, Index *jCol,
                          Number* values);

  /** Method to return:
   *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
   *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
   */
  virtual bool eval_h(Index n, const Number* x, bool new_x,
                      Number obj_factor, Index m, const Number* lambda,
                      bool new_lambda, Index nele_hess, Index* iRow,
                      Index* jCol, Number* values);

  //@}

  /** @name Solution Methods */
  //@{
  /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
  virtual void finalize_solution(SolverReturn status,
                                 Index n, const Number* x, const Number* z_L, const Number* z_U,
                                 Index m, const Number* g, const Number* lambda,
                                 Number obj_value,
				 const IpoptData* ip_data,
				 IpoptCalculatedQuantities* ip_cq);
  //@}

private:
  /**@name Methods to block default compiler methods.
   * The compiler automatically generates the following three methods.
   *  Since the default compiler implementation is generally not what
   *  you want (for all but the most simple classes), we usually 
   *  put the declarations of these methods in the private section
   *  and never implement them. This prevents the compiler from
   *  implementing an incorrect "default" behavior without us
   *  knowing. (See Scott Meyers book, "Effective C++")
   *  
   */
  //@{
  //  %(probu)s_NLP();

  double %(variable_declarations)s;

  int nU;
  int nP;
  int nY;
  int nI;
  int skip;
  double* K11val;
  double* K11val2;
  double* K11valp1;
  double* dK11val;
  double* dK11val2;
  double* dK11valp1;
  double* Xdval;
  double* Xdval2;
  double* Xdvalp1;
  double* Xval;
  double* Xval2;
  double* Xvalp1;
  double* Pval;
  double* Ival;
  double* Ival2;
  double* Ivalp1;
  double* Rvalue;
  int Time;
  double hstep;
  string buffer;
  string* specs;
  double** bounds;
  %(probu)s_NLP(const %(probu)s_NLP&);
  %(probu)s_NLP& operator=(const %(probu)s_NLP&);
  //@}
};

#endif
"""

variable_decl = ", ".join("* %s, * %sdummy" % (x,x) for x in discretize.Ldata + discretize.Lstimuli)

open(prob.lower() + "_nlp.hpp",'wt').write(template % dict(probu = prob.upper(), probl = prob.lower(),
                                                    variable_declarations = variable_decl))
