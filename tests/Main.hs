module Main where
import           Approx               ( assertApproxEqual )
import           Data.Complex         ( Complex(..) )
import           Math.Eisenstein
                                      ( eisensteinE4,
                                        eisensteinE6 )
import           Test.Tasty           ( defaultMain, testGroup )
import           Test.Tasty.HUnit     ( testCase )

i_ :: Complex Double
i_ = 0.0 :+ 1.0

z :: Complex Double
z = 1.0 :+ 1.0

tau1 :: Complex Double 
tau1 = i_

tau2 :: Complex Double 
tau2 = i_ / 10.0

tau3 :: Complex Double 
tau3 = 2.0 :+ 2.0

main :: IO ()
main = defaultMain $
  testGroup "Tests"
  [ 
    testCase "E4 is modular - condition 1" $ do
      let e4_tau   = eisensteinE4 tau1 
          e4_taup1 = eisensteinE4 (tau1 + 1)
      assertApproxEqual "" 14 e4_tau e4_taup1,

    testCase "E4 is modular - condition 2" $ do
      let e4  = eisensteinE4 (-1 / tau2) 
          e4' = tau2**4 * eisensteinE4 tau2
      assertApproxEqual "" 14 e4 e4',

    testCase "E6 is modular - condition 1" $ do
      let e6_tau   = eisensteinE6 tau2 
          e6_taup1 = eisensteinE6 (tau2 + 1)
      assertApproxEqual "" 14 e6_tau e6_taup1,

    testCase "E6 is modular - condition 2" $ do
      let e6  = eisensteinE6 (-1 / tau3) 
          e6' = tau3**6 * eisensteinE6 tau3
      assertApproxEqual "" 14 e6 e6'

  ]
