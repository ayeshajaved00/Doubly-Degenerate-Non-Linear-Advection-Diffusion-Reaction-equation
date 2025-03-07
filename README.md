## Double-Degenerate-Non-Linear-Advection-Diffusion-Reaction-equation
This repository contains numerical test cases for different nonlinear equations solved using three iterative schemes: **Newton Scheme**, **M Scheme**, and **Adaptive M Scheme**. The equations included are

1. **Porous Medium Equation**
2. **Biofilm Equation**
3. **Richards Equation with Advection**
4. **Double Degenerate Toy Model**

Each equation has been implemented for both **1D** and **2D** cases. The initial condition for all cases is taken as the **Barenblatt solution**, with the parameter  **m**  (Barenblatt parameter) set to **6**.

## Parameters

The following parameters are used throughout the numerical simulations:

- **T**: Total time (final time for the simulation)
- **dt**: Time step (the size of each time step during the simulation)
- **N**: Total number of time steps
- **x**: Space grid (the spatial grid used for the discretization)
- **ee**: Source term (external source or forcing term in the equation)
- **M (MB,Mb)**: M scheme parameter (used in the M scheme method)
- **n_values**: Number of mesh points
- **L**: Left domain boundary (the left boundary for the spatial domain)
- **R**: Right domain boundary (the right boundary for the spatial domain)

## Mathematical Formulas
**Contraction Rate**

Contraction Rate (α_i) = error_i/ error_(i-1)  ` ∀ i ∈ N`

**Convergence Order**

Convergence Order = log( α_i / α_(i-1))      ` ∀ i ∈ N`

## 1D Numerical Tests
For each equation in the 1D case, the following numerical experiments have been conducted:
- **Convergence Order**
- **Contraction Rate**
- **Error vs. Iteration Number**
- **Average Number of Iterations vs. Mesh Size**

## 2D Numerical Tests
In 2D cases, only the **Average Number of Iterations vs. Mesh Size** has been analyzed. The simulations were executed on the VSC supercomputer.

## Repository Structure
```
📂 ## Doubly-Degenerate-Non-Linear-Advection-Diffusion-Reaction-equation
│-- 📂 1D_Biofilm
│   ├── 📂 Convergence_order
│   │   ├── order_convergence.m              # MATLAB script for Convergence order
│   ├── 📂 Error_vs_Iter
│   │   ├── ErrorVsiterationNumber.m         # MATLAB script for Error vs Iteration Number
│   ├── 📂 avg_iter_vs_mesh_size
│   │   ├── Simulation.m                     # MATLAB script for Average number of iterations Vs mesh size (h)
│   ├──  📂 Contraction_Rate
│   │   ├── M_scheme.m                        # MATLAB script for M scheme
│   │   ├── MAdap_scheme.m                    # MATLAB script for adaptive M scheme
│   │   ├── Newton_scheme.m                   # MATLAB script for Newton scheme
│   │   ├── Biofilm_contraction_Arithmetic.m  # The MATLAB script that runs all the results of the schemes and plots them
│-- 📂 1D_PME
│   ├── 📂 Convergence_order                  
│   │   ├── order_convergence.m               # MATLAB script for Convergence order
│   ├──  📂 Contraction_Rate
│   │   ├── M_scheme.m                        # MATLAB script for M scheme
│   │   ├── AdapM_scheme.m                    # MATLAB script for adaptive M scheme
│   │   ├── Newton_scheme.m                   # MATLAB script for Newton scheme
│   │   ├── Contraction_rate.m                # The MATLAB script that runs all the results of the schemes and plots them
│   ├── 📂 avg_iter_vs_mesh_size
│   │   ├── Simulation.m                       # MATLAB script for Average number of iterations Vs mesh size (h)
│   ├── 📂 error_vs_iter
│   │   ├── ErrVsIter.m                        # MATLAB script for Error vs Iteration Number
│-- 📂 1D_Richards
│   ├── 📂 Avg-Iter_Vs_mesh-size
│   │   ├── Simulation.m                        # MATLAB script for Average number of iterations Vs mesh size (h)
│   ├── 📂 Contraction_Rate
│   │   ├── MScheme.m                            # MATLAB script for M scheme
│   │   ├── Richards_Newton.m                    # MATLAB script for Newton scheme
│   │   ├── Adap_Richards.m                      # MATLAB script for adaptive M scheme
│   │   ├── contraction_Arithmetic.m             # The MATLAB script that runs all the results of the schemes and plots them
│   ├── 📂 Convergence_Order
│   │   ├── order_convergenec.m                   # MATLAB script for Convergence order
│   ├── 📂 Error_Vs_Iter
│   │   ├── ErrorVsIterationNumber.m              # MATLAB script for Error vs Iteration Number
│-- 📂 1D_Toymodel
│   ├── 📂 Contraction_Rate
│   │   ├── DD_mixed.m                            # MATLAB script for M scheme
│   │   ├── DD_Newton.m                           # MATLAB script for Newton scheme
│   │   ├── DD_AdapM.m                            # MATLAB script for adaptive M scheme
│   │   ├── Contraction_plot_Arithmetic.m         # MATLAB script that runs all the results of the schemes and plots them
│   ├── 📂 Convergence_Order
│   │   ├── order_convergence.m                   # MATLAB script for convergence order
│   ├── 📂 error_vs_iter
│   │   ├── ErrorVsIterationNumber.m              # MATLAB script for Error vs Iteration Number
│   ├── 📂 average-iter_Vs_mesh-size
│   │   ├── Simulation.m                           # MATLAB script for Average number of iterations Vs mesh size (h)
│-- 📂 2D_Biofilm/Average_iter_Vs_mesh
│   ├── BioSimulation.m
│-- 📂 2D_PME/Average_iter_Vs_mesh                 
│   ├── PSimulation.m                               # MATLAB script for Average number of iterations Vs mesh size (h)
│-- 📂 2D_Richards/Average_iter_Vs_mesh
│   ├── RSimulation.m                               # MATLAB script for Average number of iterations Vs mesh size (h)
│-- 📂 2D_Toymodel/Average_iter_Vs_mesh
│   ├── DDSimulation.m                              # MATLAB script for Average number of iterations Vs mesh size (h)



