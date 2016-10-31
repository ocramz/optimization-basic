---------------------------------------------------------------------------
-- | Module    : Numeric.Optimization.Basic.LP
-- Copyright   : (c) 2016 Marco Zocca
-- License     : GPL-3
--
-- Maintainer  : zocca.marco gmail
-- Stability   : experimental
-- Portability : portable
--
-- This module implements a general interior-point method for LP, based on the logarithmic barrier penalization 
--
-- * [1] Primal-dual interior point algorithms for linear programming, G. Tzallas-Regas
--
--------------------------------------------------------------------------
module Numeric.Optimization.Basic.LP where


import Data.Sparse -- (reexports Data.Sparse.Common from sparse-linear-algebra)

import Numeric.LinearAlgebra.Sparse


{- |
Penalized LP (PLP):

max_x <c, x> - rho B(x)

s.t. A x = b

where B(x) := sum_j( log(x_j) ) 

1. Write Lagrangian for PLP L(x, rho, y)
2. Write KKT conditions (first order optimality of L(x, rho, y))
3. Write Newton equations for KKT system
-}



-- | Update primal `x`, dual `y` and slack `s` variables by applying Newton's method to the KKT system of the log-penalized LP. This requires one linear solve and a line search for a stepsize `alpha` that preserves feasibility : (x, s) > 0.
-- NB : It assumes the dimensions of vectors are x :: Rn , y :: Rm , s :: Rn 
pdIpLPstep ::
  (SpVector Double ->
   SpVector Double ->
   SpVector Double ->
   SpVector Double ->
   Double) ->
  SpMatrix Double ->
  Double ->
  (SpVector Double, SpVector Double, SpVector Double, Int) ->
  (SpVector Double, SpVector Double, SpVector Double, Int)
pdIpLPstep chooseAlpha aa tau (x, y, s, k) = (xNew, yNew, sNew, k+1) where
  gamma = (x `dot` s) / fromIntegral n  -- duality measure
  rho = tau * gamma                     -- Lagrange multiplier of barrier penalty
  ss = diagonalSM s 
  xx = diagonalSM x
  arow1 = zn -||- transposeSM aa -||- eye n
  arow2 = aa -||- zm             -||- zm
  arow3 = ss -||- zm             -||- xx
  amat = arow1 -=- (arow2 -=- arow3)
  rhs = concatSV zvn $ concatSV zvm $ negated (xx #> s) ^+^ (rho .* onesSV m)
  xysStep = amat <\> rhs
  xStep = takeSV n xysStep
  yStep = takeSV m (dropSV n xysStep)
  sStep = dropSV (m+n) xysStep
  -- xys = concatSV x $ concatSV y s
  xNew = x ^+^ (alpha .* xStep)
  yNew = y ^+^ (alpha .* yStep)
  sNew = s ^+^ (alpha .* sStep)  
  alpha = chooseAlpha x s xStep sStep
  (m, n) = dim aa
  zm = zeroSM m m
  zn = zeroSM n n
  zvn = zeroSV n
  zvm = zeroSV m

