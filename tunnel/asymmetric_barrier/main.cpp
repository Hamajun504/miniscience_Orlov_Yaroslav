#include <dolfin.h>
#include "Shrodinger.h"
#include "cmath"
#include <iostream>

using namespace dolfin;

class OutBoundary : public SubDomain
{
  double w;
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return on_boundary;
  }
public:
  OutBoundary(double w): SubDomain(), w(w) {}
};


class OutFunctionRight : public Expression
{
  double w, k;
public:
  OutFunctionRight(double w, double k) : Expression(2), w(w), k(k) {}

  void eval(Array<double>& values, const Array<double>& x) const
  {
    if (x[0] < 0) {
      values[0] = 0;
      values[1] = 0.5;
    }
    else{
        values[0] = 1;
        values[1] = 0;
    }
  }

  std::size_t value_rank()const override
  {
    return 1;
  }

  std::vector<std::size_t> value_shape()const override
  {
    return {2};
  }
};

class OutFunctionLeft : public Expression
{
  double w, k;
public:
  OutFunctionLeft(double w, double k) : Expression(2), w(w), k(k) {}

  void eval(Array<double>& values, const Array<double>& x) const
  {
    if (x[0] > 0) {
      values[0] = 0;
      values[1] = 0.5;
    }
    else{
        values[0] = 1;
        values[1] = 0;
    }
  }

  std::size_t value_rank()const override
  {
    return 1;
  }

  std::vector<std::size_t> value_shape()const override
  {
    return {2};
  }
};

class Potential : public Expression
{
  double U, w;
public:

  Potential(double U, double w) : Expression(0), U(U), w(w) {}
  //rectangle
  /* void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = x[0] > 1 && x[0] <= 1 + w ? U : 0;
  } */
  //triangle
  void eval(Array<double>& values, const Array<double>& x) const
  {
    if (x[0] < -0.01) values[0] = 0;
    else if (x[0] < -0.05) values[0] = U / 0.05 * (x[0] + 0.1);
    else if (x[0] < 0.05) values[0] = -U / 0.15 * (x[1] - 0.1);
    else values[0] = 0;
  }

  std::size_t value_rank()const override
  {
    return 0;
  }

  std::vector<std::size_t> value_shape()const override
  {
    return {1};
  }
};



int main()
{
    double Ed = 49, U = 50;
    double width = 0.2;
    std::string out = "triangle";
    std::cout << "E = " << Ed << std::endl;

    auto E = std::make_shared<Constant>(Ed);

    auto mesh = std::make_shared<IntervalMesh>(100000, -1 - width / 2, 1 + width / 2);

    auto V = std::make_shared<Shrodinger::FunctionSpace>(mesh);

    auto out_boundary = std::make_shared<OutBoundary>(width);
    auto zero = std::make_shared<Constant>(0);
    auto one = std::make_shared<Constant>(1, 0);
    auto potential = std::make_shared<Potential>(U, width);
    auto out_function_right = std::make_shared<OutFunctionRight>(width, sqrt(Ed));
    auto out_function_left = std::make_shared<OutFunctionLeft>(width, sqrt(Ed));

    DirichletBC bc_right(V, out_function_right, out_boundary);
    DirichletBC bc_left(V, out_function_left, out_boundary);
    
    auto psi = std::make_shared<Function>(V);
    
    Shrodinger::BilinearForm a(V, V);
    Shrodinger::LinearForm L(V);

    a.E = E;
    a.U = potential;
    L.f = zero;
 
    solve(a==L, *psi, bc_right);
    File file_right("results/" + out + "_right.raw");
    file_right << *psi;

    solve(a==L, *psi, bc_left);
    File file_left("results/" + out + "_left.raw");
    file_left << *psi;

    return 0;
}
