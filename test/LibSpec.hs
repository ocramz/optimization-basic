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

-}


aa0 :: SpMatrix Double
aa0 = fromListDenseSM 3 [1,4,-1,2,2,1]

b0, x0, y0, s0 :: SpVector Double
b0 = fromListDenseSV 3 [4, 12, 1]
x0 = fromListDenseSV 2 [1,1]
y0 = fromListDenseSV 2 [2,2]
s0 = fromListDenseSV 2 [3,3]

alpha0 = 1e-33

tau = 0.3

(x1, y1, s1, _) = pdIpLPstep (\_ _ _ _ -> alpha0) aa0 tau (x0, y0, s0, 1)
