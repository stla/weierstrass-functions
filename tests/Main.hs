module Main where
import           Approx               ( assertApproxEqual )
import           Data.Complex         ( Complex(..) )
import           Math.Eisenstein
                                      ( eisensteinE4,
                                        eisensteinE6,
                                        kleinJ,
                                        agm,
                                        kleinJinv )
import           Math.Gamma           ( gamma )
import           Test.Tasty           ( defaultMain, testGroup )
import           Test.Tasty.HUnit     ( testCase )
import           Math.Weierstrass     ( halfPeriods, ellipticInvariants )

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
      assertApproxEqual "" 12 e4_tau e4_taup1,

    testCase "E4 is modular - condition 2" $ do
      let e4  = eisensteinE4 (-1 / tau2) 
          e4' = tau2**4 * eisensteinE4 tau2
      assertApproxEqual "" 12 e4 e4',

    testCase "E6 is modular - condition 1" $ do
      let e6_tau   = eisensteinE6 tau2 
          e6_taup1 = eisensteinE6 (tau2 + 1)
      assertApproxEqual "" 7 e6_tau e6_taup1,

    testCase "E6 is modular - condition 2" $ do
      let e6  = eisensteinE6 (-1 / tau3) 
          e6' = tau3**6 * eisensteinE6 tau3
      assertApproxEqual "" 10 e6 e6',

    testCase "a value of Klein J-function" $ do
      let expected = 66**3
          obtained = kleinJ (2 * i_)
      assertApproxEqual "" 7 expected obtained,

    testCase "a value of agm" $ do
      let expected = 2 * pi ** 1.5 * sqrt 2 / gamma 0.25 ** 2
          obtained = agm 1 (sqrt 2)
      assertApproxEqual "" 14 expected obtained,

    testCase "kleinJ o kleinJinv = id" $ do
      let expected =  0.2 :+ 0.2
          obtained = kleinJ (kleinJinv (0.2 :+ 0.2))
      assertApproxEqual "" 12 expected obtained,

    testCase "Elliptic invariants - 1/2" $ do
      let g2 = (-7) :+ 9
          g3 = 5 :+ 3
          (omega1, omega2) = halfPeriods g2 g3
          (g2', _) = ellipticInvariants omega1 omega2
      assertApproxEqual "" 10 g2 g2',

    testCase "Elliptic invariants - 2/2" $ do
      let g2 = (-7) :+ 9
          g3 = 5 :+ 3
          (omega1, omega2) = halfPeriods g2 g3
          (_, g3') = ellipticInvariants omega1 omega2
      assertApproxEqual "" 10 g3 g3'

  ]
