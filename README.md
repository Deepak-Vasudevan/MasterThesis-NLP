# Master-Thesis
This project is an exploration attempt into the mathematical optimization of a water distribution network under the influence of leakages caused by deterioration within the distribution network.

The optimization process is carried out in two steps

      Step 1: Optimal placement of pressure valves

      Step 2: Optimal control of the valves in the presence of known leaks

The optimization framework addresses the minimisation of the average network pressure in an extended time setting by imposing the hydraulic equations as nonlinear constraints. The hydraulic components, namely the pressure reduction valves, are modelled as integer variables resulting in a non-convex and non-linear optimization problem that falls under the
class of optimization problems known as mixed-integer nonlinear programming (MINLP).

The project implements two reformulation methods that solves the
MINLP problem as a sequence of regular nonlinear programs (NLPs) - The Penalty method and Reformulation Method.

While there is sufficient research on water network optimization using various mathematical methods, this project tries to combine a leakage model within the optimization framework.

The project also includes a generation of various pseudo random demand patters for simulating a real distribution network simulated under appropriate demand conditions, though explicit details of the consumer data points are not included with the code files.

Keywords: Nonconvex optimization, NLP, MINLP, BONMIN, IPOPT, Demand management, Pressure reduction, Leakage minimisation.

Coding tool : MATLAB

Simulation Tool : Mike Urban
