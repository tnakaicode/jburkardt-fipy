# Fenics]

## dpg_laplace_adapt

- "Error:   Unable to successfully call PETSc function 'KSPSolve'."
  - PETSc error code is: 76 (Error in external library).
  - <https://fenicsproject.org/qa/4177/reason-petsc-error-code-is-76/>
- "Calling FFC just-in-time (JIT) compiler, this may take some time."
  - <https://fenicsproject.discourse.group/t/calling-ffc-just-in-time-jit-compiler-this-may-take-some-time/2275>

| mesh                                  | sol                                |
| ------------------------------------- | ---------------------------------- |
| ![pic](./dpg_laplace_adapt/msh.png)   | ![pic](./dpg_laplace_adapt/u.png)  |
| ![pic](./dpg_laplace_adapt/mesh0.png) | ![pic](./dpg_laplace_adapt/u0.png) |
| ![pic](./dpg_laplace_adapt/mesh1.png) | ![pic](./dpg_laplace_adapt/u1.png) |
| ![pic](./dpg_laplace_adapt/mesh2.png) | ![pic](./dpg_laplace_adapt/u2.png) |
| ![pic](./dpg_laplace_adapt/mesh3.png) | ![pic](./dpg_laplace_adapt/u3.png) |
| ![pic](./dpg_laplace_adapt/mesh4.png) | ![pic](./dpg_laplace_adapt/u4.png) |
| ![pic](./dpg_laplace_adapt/mesh5.png) | ![pic](./dpg_laplace_adapt/u5.png) |
| ![pic](./dpg_laplace_adapt/mesh6.png) | ![pic](./dpg_laplace_adapt/u6.png) |
| ![pic](./dpg_laplace_adapt/mesh7.png) | ![pic](./dpg_laplace_adapt/u7.png) |
| ![pic](./dpg_laplace_adapt/mesh8.png) | ![pic](./dpg_laplace_adapt/u8.png) |

## heat_explicit

![pic](./heat_explicit/heat_explicit_mesh.png)

|      num | sol                                                     |          |                                                         |
| -------: | ------------------------------------------------------- | -------: | ------------------------------------------------------- |
|    num=0 | ![pic](./heat_explicit/heat_explicit_solution_0.png)    |  num=500 | ![pic](./heat_explicit/heat_explicit_solution_500.png)  |
| num= 100 | ![pic](./heat_explicit/heat_explicit_solution_100.png)  |  num=700 | ![pic](./heat_explicit/heat_explicit_solution_700.png)  |
| num= 200 | ![pic](./heat_explicit/heat_explicit_solution_200.png)  |  num=800 | ![pic](./heat_explicit/heat_explicit_solution_800.png)  |
| num= 300 | ![pic](./heat_explicit/heat_explicit_solution_300.png)  |  num=900 | ![pic](./heat_explicit/heat_explicit_solution_900.png)  |
| num= 400 | ![pic](./heat_explicit/heat_explicit_solution_400.png ) | num=1000 | ![pic](./heat_explicit/heat_explicit_solution_1000.png) |

## heat_implicit

![pic](./heat_implicit/heat_implicit_mesh.png)

|    num | sol                                                   |    num | sol                                                   |
| -----: | ----------------------------------------------------- | -----: | ----------------------------------------------------- |
| num= 0 | ![pic](./heat_implicit/heat_implicit_solution_0.png)  | num=12 | ![pic](./heat_implicit/heat_implicit_solution_12.png) |
| num= 2 | ![pic](./heat_implicit/heat_implicit_solution_2.png)  | num=14 | ![pic](./heat_implicit/heat_implicit_solution_14.png) |
| num= 4 | ![pic](./heat_implicit/heat_implicit_solution_4.png)  | num=16 | ![pic](./heat_implicit/heat_implicit_solution_16.png) |
| num= 6 | ![pic](./heat_implicit/heat_implicit_solution_6.png)  | num=18 | ![pic](./heat_implicit/heat_implicit_solution_18.png) |
| num= 8 | ![pic](./heat_implicit/heat_implicit_solution_8.png)  | num=20 | ![pic](./heat_implicit/heat_implicit_solution_20.png) |
| num=10 | ![pic](./heat_implicit/heat_implicit_solution_10.png) |

## Step11 Error Decay for a Sequence of Adaptive Meshes

uses the Discontinuous Petrov Galerkin (DPG) method
to solve a Poisson problem, and repeatedly refines the mesh, guided by DPG error indicators

| case0                                         | case1                                         |
| --------------------------------------------- | --------------------------------------------- |
| ![pic](./step11/step11_error_decay_case0.png) | ![pic](./step11/step11_error_decay_case1.png) |
| ![pic](./step11/step11_case0_mesh_level0.png) | ![pic](./step11/step11_case1_mesh_level0.png) |
| ![pic](./step11/step11_case0_mesh_level1.png) | ![pic](./step11/step11_case1_mesh_level1.png) |
| ![pic](./step11/step11_case0_mesh_level2.png) | ![pic](./step11/step11_case1_mesh_level2.png) |
| ![pic](./step11/step11_case0_mesh_level3.png) | ![pic](./step11/step11_case1_mesh_level3.png) |
| ![pic](./step11/step11_case0_mesh_level4.png) | ![pic](./step11/step11_case1_mesh_level4.png) |
