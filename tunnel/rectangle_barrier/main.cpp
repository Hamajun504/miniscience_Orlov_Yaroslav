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
    return (x[0] > 2 + w - DOLFIN_EPS);
  }
public:
  OutBoundary(double w): SubDomain(), w(w) {}
};

class OutFunction : public Expression
{
  double w, k;
public:
  OutFunction(double w, double k) : Expression(2), w(w), k(k) {}

  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = cos(k*(x[0] - 1 - w));
    values[1] = sin(k*(x[0] - 1 - w));
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

  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = x[0] > 1 && x[0] <= 1 + w ? U : 0;
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
    double Ed, U;
    std::string width;
    std::cin >> Ed >> U >> width;
    std::cout << "E = " << Ed << std::endl;

    auto E = std::make_shared<Constant>(Ed);

    auto tentative_mesh = std::make_shared<IntervalMesh>(100000, 0, 2 + stod(width) * 0.1);

    auto mesh = tentative_mesh;

    auto V = std::make_shared<Shrodinger::FunctionSpace>(mesh);

    auto out_boundary = std::make_shared<OutBoundary>(stod(width) * 0.1);
    auto zero = std::make_shared<Constant>(0);
    auto one = std::make_shared<Constant>(1, 0);
    auto potential = std::make_shared<Potential>(U, stod(width) * 0.1);
    auto out_function = std::make_shared<OutFunction>(stod(width) * 0.1, sqrt(Ed));

    DirichletBC bc_in(V, out_function, out_boundary);
    std::vector<const DirichletBC*> bcs{&bc_in};
    
    auto psi = std::make_shared<Function>(V);
    
    Shrodinger::BilinearForm a(V, V);
    Shrodinger::LinearForm L(V);

    a.E = E;
    a.U = potential;
    L.f = zero;

    
    solve(a==L, *psi, bcs);
    std::string out;
    {
      out = "w";
      out += width;
      out += "/U";
      out += std::to_string(U);
      out += "E";
      out += std::to_string(Ed);
	    out += ".raw";
    }
    File file(out);
    file << *psi;

    return 0;
}
