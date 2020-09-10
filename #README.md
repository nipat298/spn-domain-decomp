# SPN approximation codes

Readme for codes for SPN approximation based modeling of the radiative transfer through a medium. 
Supports non-homogenous media with varying number of inhomogeneities. 

Uses a tetrahedral/ triangular mesh and uses the finite element method for solving the forward problem. Presently supports only cuboidal domains. (to include support for other mesh generators/ adaptive meshing in future).

Supports fluorescence and elastic scattering simulations.
The medium properties are provided by the user through a `SettingsFile` . A sample file has been included for both fluorescence and elastic scattering simulation.

[TOC]

## The `SettingsFile`
This file is the only essential input from the user end. This allow the user to specify the problem as 3D or 2D, elastic or fluorescence, and specify domain size, inclusions, source and detector specifications. The source and detector are specified through structure quantities.  Detailed documentation is in the sample files provided.

Presently support only single type of inhomogeneity. (Note: this can be changed to allow different types of inhomogeneities by using a structure to specify the properties of the inhomogeneity)

Check `sample_fl.m` (for fluorescence) and `sample_el.m` (for elastic scattering)

## Meshing and Assembly data

We use MATLAB's delaunayTriangulation function to generate the triangulation mesh. Once the mesh is generated appropriate element and node numbers are assigned to different entities such as source, detectors and inclusions. All this is encapsulated in the function `MeshGen2D_new`  for 2D and `MeshGen3D` for 3D. While the present version has been tested for 2D, and a 3D framework is available, testing with 3D is incomplete.

For a regular mesh, the mass and stiffness matrices can be precomputed and stored. This is done with the help of the functions `Assemble2D_new` and `Assemble3D` for 2D and 3D respectively. 

## Data generation

The function `DataGen` can be used to generate measurement Data for reconstruction studies. It takes as inputs details of the problem through the [SettingsFile](##The-SettingsFile) and returns the exiting partial current . Noise can be also added to the measurements by specifying a desired SNR. The simulation parameters can be set to choose the desired [Forward Solver](##Forward-Solvers).

## Forward Solvers

The forward solver is based on the SPN approximation to the radiative transfer equation. Three solver modules are presently supported namely the basic SPN approximation [BSPN](###Forward-Solver-for-BSPN) ( as in Klose et al), the SPN approximation with the first collision source [FSPN](###Forward-Solver-for-FSPN) and the SPN approximation with the delta-Eddington approximation to the phase function [DSPN](###Forward-Solver-for-DSPN).  By suitable choice of the input arguments the same code can be used for fluorescence and elastic scattering simulations. The output in each case is the exiting partial current and optionally Fluence and/or supplementary data for the reconstruction routines. 

Each of the forward solver codes uses four main subroutines:  `SP3Formulation` , `SpnFwd`, `GetSource` and `GetPartCurr`. Other subroutines used are listed as [auxiliary subroutines](####Other-auxiliary-subroutines-used).

#### The `SP3Formulation` subroutine

This subroutine uses the information about the optical properties and Source-detector positions to form the elementwise SP3 approximation coefficient matrices. Details can be found in (include reference here)

#### The `SpnFwd` subroutine

This subroutine makes use of the coefficient matrices and the assembly matrices to generate the global FEM assembly matrix which is used subsequently in the linear solver to solve for Fluence. 

#### The `GetSource` subroutine

This Subroutine generates the nodal source vector corresponding to each source. It makes use subroutines `GetSource2D_internal` for 2D (`GetSource3D_internal` for 3D) for internal sources (to be tested) and subroutines `GetSourceX2D` for 2D (`GetSourceX3DN` for 3D) for a boundary source. The 3D codes are yet to be tested. 

#### The `GetPartCurr` subroutine

This subroutine is used to generate the operator matrix that will operate on the fluence to give the exiting partial current values. It makes use of function `GetExitCurrent2D` (2D) or `GetExitCurrent3D` (3D) to obtain the operator matrices. 

#### Other auxiliary subroutines used 

**`CalcReflCoeff`** - function to calculate coefficients dependent on the reflection coefficient.

**`DefineBoundMat`** - function used to define the coefficient matrices for the boundary condition when the external refractive index is different across different interfaces. 

**`CalcDetCoeff`** - function to calculate refractive index dependent coefficients for measurement matrices.

**`DirectedDetection`** - function to calculate refractive index dependent coefficients for measurement matrices when detector is not isotropic.

**`GetAssembledMat`** - generate global FEM assembly matrix for given pair of assembly and coefficient matrix.

**`GetSourceX2D`** - generate the source vector for a point source on the boundary.

**`GetDetFluence`** - obtain detected fluence on boundary detectors of either point type or surface detectors. 

**`get_raytrace_2D`** - this file is used to obtain the trace of the collimated ray from the point of incidence to the point of exit from the medium for collimated light sources only. The elements that the ray crosses are saved in the source structure as well as the coordinates of the point of intersections for each element and the length of the ray segment that cuts through the element. The source is now marked as an internal source instead of a boundary source.
**`get_pfunc_coeff`** - this script is used to choose the phase function. Presently three phase functions are supported, Henyey-Greenstein, Modified-Henyey-Greenstein and Two-Term-Henyey-Greenstein. A flag allows to choose between their conventional or delta-Eddington approximations. 

### Forward Solver for BSPN

Use the function `FwdSolver` to run a simulation using the Basic SPN approximation. Source type can be **collimated** or **isotropic** source. Only **point** sources are presently supported. Support for beam sources is to be added. 

### Forward Solver for FSPN

Use the function `FwdSolver_fcs_4` to run a simulation using the SPN approximation in conjunction with the first collision source. This is helpful to better approximate the anisotropy in radiance when simulating collimated/ directed sources. This simulation type only supports sources of type **collimated**. Presently only point sources are supported, and support for beam-type sources can easily be incorporated. (Not tested yet).

The first collision source is computed using the subroutine `fwdsolver_unscattered_v7`. This routine evaluates the source strength at each element along the ray as well as its spatial derivative which is required while computing the source vector through the routine `GetSource2D_fcs_2`. Details of the derivation of the FSPN approximation are given in report 1 (Insert dated reference to report here). 

### Forward Solver for DSPN

Use the function `FwdSolver_fcs_de_3` to run a simulation using the SPN approximation in conjunction with the first collision source and the delta-Eddington approximation to the phase function. This improves on the **FSPN** by taking into account anisotropy of the phase function. 

The first collision source is computed using the subroutine `fwdsolver_unscattered_de_3`. This routine evaluates the source strength at each element along the ray as well as its spatial derivative which is required while computing the source vector through the routine `GetSource2D_fcs_de_2`. Details of the derivation of the DSPN approximation are given in report 1 (Insert dated reference to report here). 

## Adjoint sensitvity

Frechet derivatives for each forward solver are computed through the routines `FrechetDeriv1`, `FrechetDeriv1_fcs_3`and `FrechetDeriv1_de_3` for BSPN, FSPN and DSPN variants of the solver respectively. Each routine consists of 3 main  sections:

1. **Evaluating the adjoints**: The subroutine `AdjointFwd` is used to generate the coefficient matrices for the adjoint problem and the adjoint fields are evalutated using a linear solver. 
2. **Computing the derivative of the coefficient matrices**: The subroutine `dSP3Formulation` (for BSPN, FSPN and `dSP3Formulation_de` for DSPN) is used to compute the derivatives of the coefficient matrices for the forward problem wrt the parameter of interest i.e. absorption coefficient (for elastic scattering) or fluorophore absorption coefficient (for fluorescence). 
   In addition the subroutine `GetDiffSource2D_fcs_3` and `GetDiffSource2D_de_3` are used to compute the derivative of the source term in the FSPN and DSPN variants respectively.
3. **Assembling the Jacobian**: A vectorized implementation is used to assemble the Jacobian using the derivatives of the coefficient matrices and the forward and adjoint fields. 

`TestJac2.m` is a sample script for comparing the adjoint sensitivity with the finite difference evaluation of the same and plotting these comparisons.

## Reconstruction Algorithms

The reconstruction routine `ResidualMinAlgo` allows the user to choose between four different gradient based residual least squares minimisation algorithm, 

	1. The first order Levenberg-Marquardt. (`folm_step`)
 	2. The second order Levenberg_marquardt, which is a second-degree scheme based on the predictor-corrector method of Hettlich and Rundell (`solm_step_optim`)
 	3. The Levenberg-Marquardt scheme with geodesic acceleration (`galm_step_optim`)
 	4. The Tensor-Newton scheme. (`tnlm_step_FO_optim2_2`)

The routine also allows the user to opt for either a backtracking line search or a quadratic-cubic line search through the subroutine `LineSearch`. The Levenberg-Marquardt parameter in each case is updated using the subroutine `update_lambda` using one of three strategies:

1.  a continuous update strategy by Madsen
2. the 2:3 update rule
3. or the 3:5 update rule

The subroutine `eval_residual` is used to compute the residual wrt measurement data provided and the subroutines `eval_fderiv` and `eval_sderiv_optim` are used to compute the first and second order derivatives of the residual respectively.

## Plotting routine

A routine `Plot_mu2` is provided to plot the desired parameter of interest and a chosen cross section. 

## References

The following is a list of references that are helpful in understanding various sections of the code. 
(to be populated)

## Created by

Nishigandha Patil
Department of Electrical Engineering
Indian Institute of Technology Kanpur

- Primary Email: <nipat@iitk.ac.in>
- Alternate Email: nishi.rpatil@gmail.com

# Self Notes

## Things to add in readme

- Sample usage for subroutines
- References
- Quick start guide, demo script
- Sample problems

## Things to add in the codes

- Support 3D in present version
- Higher order reconstruction schemes 
- Demo script
- Include Intermediate print commands to update status of simulation. 
- Sample problems