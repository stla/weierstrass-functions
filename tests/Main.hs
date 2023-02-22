module Main where
import           Approx               ( assertApproxEqual )
import           Data.Complex         ( Complex(..) )
import           Math.Eisenstein      ( eisensteinE4,
                                        eisensteinE6,
                                        kleinJ,
                                        agm,
                                        kleinJinv, 
                                        etaDedekind,
                                        lambda )
import           Math.Gamma           ( gamma )
import           Test.Tasty           ( defaultMain, testGroup )
import           Test.Tasty.HUnit     ( testCase )
import           Math.Weierstrass     ( halfPeriods, 
                                        ellipticInvariants,
                                        weierstrassP,
                                        weierstrassPdash,
                                        weierstrassPinv )

i_ :: Complex Double
i_ = 0.0 :+ 1.0

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
      assertApproxEqual "" 12 g2 g2',

    testCase "Elliptic invariants - 2/2" $ do
      let g2 = (-7) :+ 9
          g3 = 5 :+ 3
          (omega1, omega2) = halfPeriods g2 g3
          (_, g3') = ellipticInvariants omega1 omega2
      assertApproxEqual "" 12 g3 g3',

    testCase "a value of weierstrassP" $ do
      let z = 0.1 :+ 0.1
          g2 = 2 :+ 1
          g3 = 2 :+ (-1)
          obtained = weierstrassP z g2 g3
          expected = (-0.0010285443715) :+ (-49.9979857342848)
      assertApproxEqual "" 11 expected obtained,

    testCase "Equianharmonic case" $ do
      let omega2 = gamma (1/3) ** 3 / 4 / pi
          z0 = omega2 * (1 :+ (1 / sqrt 3))
          obtained = weierstrassP z0 0 1
          expected = 0
      assertApproxEqual "" 13 obtained expected,

    testCase "Differential equation" $ do
      let z = 1 :+ 1
          g2 = 2 :+ 1
          g3 = 2 :+ (-1)
          w = weierstrassP z g2 g3
          wdash = weierstrassPdash z g2 g3
          left = wdash ** 2
          right = 4 * w ** 3 - g2 * w - g3
      assertApproxEqual "" 11 left right,

    testCase "weierstrassPinv works" $ do
      let w = 0.1 :+ 1
          g2 = 2 :+ 2
          g3 = 0 :+ 3
          z = weierstrassPinv w g2 g3
          obtained = weierstrassP z g2 g3
          expected = w
      assertApproxEqual "" 13 expected obtained,

    testCase "a value of Dedekind eta" $ do
      let expected = gamma 0.25 / 2 ** (11/8) / pi ** 0.75
          obtained = etaDedekind (2 * i_)
      assertApproxEqual "" 14 expected obtained,

    testCase "lambda modular identity" $ do
      let x = sqrt 2
          expected = 1
          obtained = lambda (i_ * x) + lambda (i_ / x)
      assertApproxEqual "" 14 expected obtained

  ]
