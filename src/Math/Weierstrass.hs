module Math.Weierstrass
    ( halfPeriods,
      ellipticInvariants,
      weierstrassP
    ) where
import           Data.Complex     ( Complex(..) )
import           Internal         ( (%^%) )
import           Math.Eisenstein  ( eisensteinE4, eisensteinE6, kleinJinv) 
import           Math.JacobiTheta ( jtheta2, jtheta3, jtheta1, jtheta4 )
import           Math.Gamma       ( gamma )


i_ :: Complex Double
i_ = 0.0 :+ 1.0

eisensteinG4 :: Complex Double -> Complex Double
eisensteinG4 tau = pi * pi * pi * pi / 45 * eisensteinE4 tau

eisensteinG6_over_eisensteinG4 :: Complex Double -> Complex Double
eisensteinG6_over_eisensteinG4 tau = 
  2 * pi * pi / 21 * eisensteinE6 tau / eisensteinE4 tau

omega1_and_tau :: 
  Complex Double -> Complex Double -> (Complex Double, Complex Double)
omega1_and_tau g2 g3 = (omega1, tau)
  where
    (omega1, tau) 
      | g2 == 0 = 
        (
          gamma (1/3) %^% 3 / (4 * pi * g3 ** (1/6)),
          0.5 :+ (sqrt 3 / 2)
        )
      | g3 == 0 = 
        (
          i_ * sqrt(sqrt(3.75 * eisensteinG4 tau' / g2)),
          tau'
        )
      | otherwise = 
        (
          sqrt(7 * g2 * eisensteinG6_over_eisensteinG4 tau' / (12 * g3)),
          tau'
        )
      where 
        g2cube = g2 * g2 * g2
        j = 1728 * g2cube / (g2cube - 27 * g3 * g3)
        tau' = kleinJinv j

-- | Half-periods from elliptic invariants.
halfPeriods :: 
    Complex Double -- ^ g2
 -> Complex Double -- ^ g3
 -> (Complex Double, Complex Double) -- ^ omega1, omega2
halfPeriods g2 g3 = (omega1, tau * omega1)
  where
    (omega1, tau) = omega1_and_tau g2 g3

g_from_omega1_and_tau :: 
  Complex Double -> Complex Double -> (Complex Double, Complex Double)
g_from_omega1_and_tau omega1 tau = (g2, g3)
  where
    q = exp (i_ *pi * tau)
    j2 = jtheta2 0 q
    j3 = jtheta3 0 q
    j2pow4  = j2 %^% 4
    j2pow8  = j2pow4 * j2pow4
    j2pow12 = j2pow4 * j2pow8
    j3pow4  = j3 %^% 4
    j3pow8  = j3pow4 * j3pow4
    j3pow12 = j3pow4 * j3pow8
    g2 = 4/3 * (pi / 2 / omega1) %^% 4 * (j2pow8 - j2pow4 * j3pow4 + j3pow8)
    g3 = 8/27 * (pi / 2 / omega1) %^% 6 *
      (j2pow12 - ((1.5 * j2pow8 * j3pow4) + (1.5 * j2pow4 * j3pow8)) + j3pow12)

-- | Elliptic invariants from half-periods.
ellipticInvariants :: 
    Complex Double -- ^ omega1
 -> Complex Double -- ^ omega2
 -> (Complex Double, Complex Double) -- ^ g2, g3
ellipticInvariants omega1 omega2 = 
  g_from_omega1_and_tau omega1 (omega2 / omega1)

weierstrassP_from_tau :: Complex Double -> Complex Double -> Complex Double
weierstrassP_from_tau z tau = 
  (pi * j2 * j3 * j4 / j1) %^% 2 - pi*pi * (j2 %^% 4 + j3 %^% 4) / 3
  where
    q = exp (i_ *pi * tau)
    j2 = jtheta2 0 q
    j3 = jtheta3 0 q
    j1 = jtheta1 z q
    j4 = jtheta4 z q

weierstrassP_from_omega :: 
  Complex Double -> Complex Double -> Complex Double -> Complex Double
weierstrassP_from_omega z omega1 omega2 = 
  weierstrassP_from_tau (z/omega1/2) (omega2/omega1) / (4 * omega1 * omega1)

-- | Weierstrass p-function
weierstrassP ::
    Complex Double -- ^ z
 -> Complex Double -- ^ elliptic invariant g2
 -> Complex Double -- ^ elliptic invariant g3
 -> Complex Double
weierstrassP z g2 g3 = weierstrassP_from_omega z omega1 omega2
  where
    (omega1, omega2) = halfPeriods g2 g3

