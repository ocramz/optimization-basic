module Data.Sparse (module X,
                    takeSV, dropSV
                   ) where


import qualified Data.IntMap as IM
import Data.Sparse.Common as X -- from sparse-linear-algebra



takeSV, dropSV :: Int -> SpVector a -> SpVector a
takeSV n (SV _ sv) = SV n $ IM.filterWithKey (\i _ -> i < n) sv

dropSV n (SV n0 sv) = SV (n0 - n) $ IM.mapKeys (subtract n) $ IM.filterWithKey (\i _ -> i >= n) sv









