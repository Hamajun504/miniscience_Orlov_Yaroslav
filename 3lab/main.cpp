// Hyperelasticity (C++)
// =====================
// 
// 
// Background
// ----------
// 
// See the section :ref:`hyperelasticity` for some mathematical
// background on this demo.
// 
// 
// Implementation
// --------------
// 
// The implementation is split in two files: a form file containing the
// definition of the variational forms expressed in UFL and a C++ file
// containing the actual solver.
// 
// Running this demo requires the files: :download:`main.cpp`,
// :download:`HyperElasticity.ufl` and :download:`CMakeLists.txt`.
// 
// 
// UFL form file
// ^^^^^^^^^^^^^
// 
// The UFL file is implemented in :download:`Hyperelasticity.ufl`, and
// the explanation of the UFL file can be found at :doc:`here
// <Hyperelasticity.ufl>`.
// 
// 
// C++ program
// ^^^^^^^^^^^
// 
// The main solver is implemented in the :download:`main.cpp` file.
// 
// At the top, we include the DOLFIN header file and the generated header
// file "HyperElasticity.h" containing the variational forms and function
// spaces.  For convenience we also include the DOLFIN namespace.
// 
// .. code-block:: cpp

#include <dolfin.h>
#include "HyperElasticity.h"

using namespace dolfin;

// We begin by defining two classes, deriving from :cpp:class:`SubDomain`
// for later use when specifying domains for the boundary conditions.
// 
// .. code-block:: cpp

const double R = 1;
const double r = 0.35; 
class Quater : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return ((x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + R * R - r * r) *
            (x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + R * R - r * r) - 
            4 * R * R * (x[0] * x[0] + x[1] * x[1]) < DOLFIN_EPS) && x[0] <= 0 && x[1] <= 0 && on_boundary;
  }
};

// We also define two classes, deriving from :cpp:class:`Expression`, for
// later use when specifying values for the boundary conditions.
// 
// .. code-block:: cpp

// Dirichlet boundary condition for clamp at left end
class Clamp : public Expression
{
public:

  Clamp() : Expression(3) {}

  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = 0.0;
    values[1] = 0.0;
    values[2] = 0.0;
  }

};

// Next:
// 
// .. code-block:: cpp

int main()
{

// Inside the ``main`` function, we begin by defining a tetrahedral mesh
// of the domain and the function space on this mesh. Here, we choose to
// create a unit cube mesh with 25 ( = 24 + 1) verices in one direction
// and 17 ( = 16 + 1) vertices in the other two directions. With this
// mesh, we initialize the (finite element) function space defined by the
// generated code.
// 
// .. code-block:: cpp

  // Create mesh and define function space
  auto mesh = std::make_shared<Mesh>("torus.xml");
  auto V = std::make_shared<HyperElasticity::FunctionSpace>(mesh);

// Now, the Dirichlet boundary conditions can be created using the class
// :cpp:class:`DirichletBC`, the previously initialized
// :cpp:class:`FunctionSpace` ``V`` and instances of the previously
// listed classes ``Left`` (for the left boundary) and ``Right`` (for the
// right boundary), and ``Clamp`` (for the value on the left boundary)
// and ``Rotation`` (for the value on the right boundary).
// 
// .. code-block:: cpp

  // Define Dirichlet boundaries
  auto left_half = std::make_shared<Quater>();

  // Define Dirichlet boundary functions
  auto c = std::make_shared<Clamp>();
  // auto r = std::make_shared<Rotation>();

  // Create Dirichlet boundary conditions
  DirichletBC bcl(V, c, left_half);
  //std::vector<const DirichletBC*> bcs = {{&bcl}};

// The two boundary conditions are collected in the container ``bcs``.
// 
// We use two instances of the class :cpp:class:`Constant` to define the
// source ``B`` and the traction ``T``.
// 
// .. code-block:: cpp

  // Define source and boundary traction functions
  auto B = std::make_shared<Constant>(0.0, 0.0, -0.5);
  auto T = std::make_shared<Constant>(0.0,  0.0, 0.0);

// The solution for the displacement will be an instance of the class
// :cpp:class:`Function`, living in the function space ``V``; we define
// it here:
// 
// .. code-block:: cpp

  // Define solution function
  auto u = std::make_shared<Function>(V);

// Next, we set the material parameters
// 
// .. code-block:: cpp

  // Set material parameters
  const double E  = 10.0;
  const double nu = 0.3;
  auto mu = std::make_shared<Constant>(E/(2*(1 + nu)));
  auto lambda = std::make_shared<Constant>(E*nu/((1 + nu)*(1 - 2*nu)));

// Now, we can initialize the bilinear and linear forms (``a``, ``L``)
// using the previously defined :cpp:class:`FunctionSpace` ``V``. We
// attach the material parameters and previously initialized functions to
// the forms.
// 
// .. code-block:: cpp

  // Create (linear) form defining (nonlinear) variational problem
  HyperElasticity::ResidualForm F(V);
  F.mu = mu; F.lmbda = lambda; F.u = u;
  F.B = B; F.T = T;

  // Create Jacobian dF = F' (for use in nonlinear solver).
  HyperElasticity::JacobianForm J(V, V);
  J.mu = mu; J.lmbda = lambda; J.u = u;

// Now, we have specified the variational forms and can consider the
// solution of the variational problem.
// 
// .. code-block:: cpp

  // Solve nonlinear variational problem F(u; v) = 0
  Parameters params("nonlinear_variational_solver");
  Parameters newton_params("newton_solver");
  newton_params.add("relaxation_parameter", 0.3);
  newton_params.add("maximum_iterations", 500);
  params.add(newton_params);
  solve(F == 0, *u, bcl, J, params);

// Finally, the solution ``u`` is saved to a file named
// ``displacement.pvd`` in VTK format.
// 
// .. code-block:: cpp

  // Save solution in VTK format
  File file("displacement.pvd");
  file << *u;

  return 0;
}
