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


import Control.Monad.State


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



-- | Update primal `x`, dual `y` and slack `s` variables by applying Newton's method to the KKT system of the log-penalized LP.
-- This requires one linear solve and a line search for a stepsize `alpha` that preserves feasibility : (x, s) > 0.
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
  xStep = takeSV n xysStep            -- 0 .. n-1            
  yStep = takeSV m (dropSV n xysStep) -- n .. n+m-1
  sStep = dropSV (m + n) xysStep      -- n+m .. 2n+m 
  xNew = x ^+^ (alpha .* xStep)
  yNew = y ^+^ (alpha .* yStep)
  sNew = s ^+^ (alpha .* sStep)  
  alpha = chooseAlpha x s xStep sStep
  (m, n) = dim aa
  zm = zeroSM m m
  zn = zeroSM n n
  zvn = zeroSV n
  zvm = zeroSV m



feasX :: (Foldable t, Ord a, Num a) => t a -> Bool
feasX = all (> 0)



linesearch amax x s dx ds = undefined where



-- * Karmarkar algorithm


  
karmIpLP :: SpMatrix Double -> SpVector Double -> KarmLPData
karmIpLP aa c = execState (untilConverged (\x -> c `dot` _xKarm x) (karmarkarLPstep aa c ar)) kInit
  where
   ar = alpha * r
   kInit = KarmLPData x0 y0 0
   (m, n) = dimSM aa
   n' = fromIntegral n
   r = 1 / sqrt (n' * (n'-1))
   alpha = (n' - 1)/ (3*n')
   x0 = recip (fromIntegral n) .* onesSV n
   y0 = x0

data KarmLPData =
  KarmLPData { _xKarm, _yKarm:: SpVector Double, _niter :: Int } deriving (Eq, Show)
  
  
-- | alphar = alpha * r
karmarkarLPstep ::
  SpMatrix Double -> SpVector Double -> Double -> KarmLPData -> KarmLPData
karmarkarLPstep aa c alphar (KarmLPData x y k) = KarmLPData xnew ynew (k + 1) where
  (m, n) = dim aa
  ddk = diagonalSM x
  pp = (aa #~# ddk) -=- svToSM (onesSV n)
  cdk = c <# ddk
  cp = cdk ^-^ (transposeSM pp #> ((pp ##^ pp) <\> (pp #> cdk)))
  cp0 = normalize 2 cp
  ynew = y ^-^ (alphar .* cp0)
  dy = ddk #> ynew
  xnew = recip (sum dy) .* dy
  -- z = c `dot` xnew



-- linObjConverged c = modifyInspectN 200 (fproj c _xKarm) where
--   fproj c f [x1, x2] = almostZero $ c `dot` f x1 - c `dot` f x2

---

-- untilConverged :: MonadState a m => (a -> SpVector Double) -> (a -> a) -> m a
untilConverged fproj = modifyInspectN 200 (normDiffConverged fproj)

normDiffConverged f [x1,x2] = almostZero $ f x1 - f x2

-- normDiffConverged fp xx = nearZero $ normSq (foldrMap fp (^-^) (zeroSV 0) xx)
