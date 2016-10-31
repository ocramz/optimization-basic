module LibSpec where

import Test.Hspec
-- import Test.Hspec.QuickCheck

-- import Lib 

import Data.Sparse
import Numeric.Optimization.Basic.LP


main :: IO ()
main = hspec spec

spec :: Spec
spec =
  describe "Lib" $ do
    it "works" $ do
      True `shouldBe` True
    -- prop "ourAdd is commutative" $ \x y ->
    --   ourAdd x y `shouldBe` ourAdd y x




-- | LP

{- | p.1

max x1 + x2

s.t.

[1,  2]    [4 ]
[4,  2] <= [12]
[-1, 1]    [1 ]

x1, x2 >= 0

xOpt = [8/3, 2/3]

-}


aa0 :: SpMatrix Double
aa0 = fromListDenseSM 3 [1,4,-1,2,2,1]

c0, b0, x0, y0, s0, xOpt :: SpVector Double
c0 = onesSV 2 
b0 = fromListDenseSV 3 [4, 12, 1]
x0 = fromListDenseSV 2 [1,1]
y0 = fromListDenseSV 2 [2,2]
s0 = fromListDenseSV 2 [3,3]
-- optimum
xOpt = fromListDenseSV 2 [8/3, 2/3]

-- smoke test
alpha0 = 0.5
tau = 0.3
(x0', y0', s0', _) = pdIpLPstep (\_ _ _ _ -> alpha0) aa0 tau (x0, y0, s0, 1)


testAmat aa x s = amat where
  (m,n) = dim aa
  ss = diagonalSM s 
  xx = diagonalSM x  
  arow1 = zn -||- transposeSM aa -||- eye n
  arow2 = aa -||- zm             -||- zm
  arow3 = ss -||- zm             -||- xx
  amat = arow1 -=- (arow2 -=- arow3)
  zm = zeroSM m m
  zn = zeroSM n n


--


{- | p.2

max c^T x
c = [0, 2, -1]

s.t.

[1,  -2, 1]    [4 ]
[1,  1,  1] <= [12]
[-1, 1]    [1 ]

x1, x2 >= 0

xOpt = [0, 1/3, 2/3]

-}

-- karmTest aa c = karmarkarLPstep aa c alphar (x, y, k) where
--   k = 0
--   alphar = alpha * r
--   (x, y, alpha, r) = karmarkarLPinit aa c


aa1 :: SpMatrix Double
aa1 = fromListDenseSM 2 [1,1,-2,1,1,1]

c1, b1, xOpt1 :: SpVector Double
c1 = fromListDenseSV 3 [0, 2, -1]
b1 = fromListDenseSV 2 [0, 1]

-- optimum
xOpt1 = fromListDenseSV 3 [0, 1/3, 2/3]


(KarmLPData x y niter) = karmIpLP aa1 c1
