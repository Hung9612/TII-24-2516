
---------------------------------------File Tree-----------------------------------------------

--./SSRD/Derivation_ssrd/Dmat  for derivation of d_ij
--./SSRD/Derivation_ssrd/Transform_func_numerical  for derivation of linearized homogenous transformation matrix (HTM) which is based on local product of exponential.
--./SSRD/Derivation_ssrd/best_Tmat_cal  for derivation of Lemmas.
--./SSRD/Derivation_ssrd/encode  for derivation of derivation of linearized HTMs in encoded space.
--./SSRD/Derivation_ssrd/parameterVar  for define parameter variables.
--./SSRD/Derivation_ssrd/simplification  for eliminate the redundant multivariate polynomials.
--./SSRD/Derivation_ssrd/torque  for eliminate the redundant parameters by full rank decomposition.

--./otherTest  for case study test.

--./genfile  for generate the code.

--./FunGen  for function to generate the code.

--./testExample  for an usual example to shown the proposed method to efficient derivation of linear-in-parameter dynamics.


---------------------------------------TEST-----------------------------------------------
%% Cofig Package: install matlab mingw C++ complier and matlab symbolic toolbox


%% Add the path of the following folder and their subfolders into path of matlab workspace 
-- FunGen
-- genfile
-- SSRD
-- otherTest

%% 
-- run "test.m" in the /testExample


