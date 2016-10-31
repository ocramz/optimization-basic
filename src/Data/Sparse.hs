module Data.Sparse (module X,
                    takeSV, dropSV,
                    diagonalSM) where


import qualified Data.IntMap as IM
import Data.Sparse.Common as X -- from sparse-linear-algebra





takeSV, dropSV :: Int -> SpVector a -> SpVector a
takeSV n (SV _ sv) = SV n $ IM.filterWithKey (\i _ -> i < n) sv

dropSV n (SV n0 sv) = SV (n0 - n) $ IM.mapKeys (subtract n) $ IM.filterWithKey (\i _ -> i >= n) sv


-- * Present in sparse-linear-algebra-0.2.1.0 :

-- | Indexed fold over SpVector
ifoldSV :: (Int -> a -> b -> b) -> b -> SpVector a -> b
ifoldSV f e (SV d im) = IM.foldWithKey f e im

-- | Fill the diagonal of a SpMatrix with the components of a SpVector
diagonalSM :: SpVector a -> SpMatrix a
diagonalSM sv = ifoldSV iins (zeroSM n n) sv where
  n = dim sv
  iins i = insertSpMatrix i i






