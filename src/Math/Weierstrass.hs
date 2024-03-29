{-|
Module      : Math.Eisenstein
Description : Some Weierstrass functions.
Copyright   : (c) Stéphane Laurent, 2023
License     : BSD3
Maintainer  : laurent_step@outlook.fr

Provides some Weierstrass functions and related functions.
-}
module Math.Weierstrass
    ( halfPeriods,
      ellipticInvariants,
      weierstrassP,
      weierstrassP',
      weierstrassPdash,
      weierstrassPdash',
      weierstrassPinv,
      weierstrassPinv',
      weierstrassSigma,
      weierstrassSigma',
      weierstrassZeta,
      weierstrassZeta'
    ) where
import           Data.Complex           ( Complex(..) )
import           Internal               ( (%^%) )
import           Math.Eisenstein        ( eisensteinE4, 
                                          eisensteinE6, 
                                          kleinJinv, 
                                          jtheta1DashDashDash0 ) 
import           Math.JacobiTheta       ( jtheta1', 
                                          jtheta2', 
                                          jtheta3', 
                                          jtheta4',
                                          jtheta1Dash0,
                                          jtheta1Dash )
import           Math.Gamma             ( gamma )
import           Math.EllipticIntegrals ( carlsonRF' )



i_ :: Complex Double
i_ = 0.0 :+ 1.0

eisensteinG4 :: Complex Double -> Complex Double
eisensteinG4 tau = pi %^% 4 / 45 * eisensteinE4 tau

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
    j2 = jtheta2' 0 tau
    j3 = jtheta3' 0 tau
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
  (pi * j2 * j3 * j4 / j1) %^% 2 - pi * pi * (j2 %^% 4 + j3 %^% 4) / 3
  where
    j2 = jtheta2' 0 tau
    j3 = jtheta3' 0 tau
    z' = pi * z
    j1 = jtheta1' z' tau
    j4 = jtheta4' z' tau

-- | Weierstrass p-function given the half-periods
weierstrassP ::
    Complex Double -- ^ z
 -> Complex Double -- ^ half-period omega1
 -> Complex Double -- ^ half-period omega2
 -> Complex Double
weierstrassP z omega1 omega2 = 
  weierstrassP_from_tau 
    (z / omega1 / 2) (omega2 / omega1) / (4 * omega1 * omega1)

-- | Weierstrass p-function given the elliptic invariants
weierstrassP' ::
    Complex Double -- ^ z
 -> Complex Double -- ^ elliptic invariant g2
 -> Complex Double -- ^ elliptic invariant g3
 -> Complex Double
weierstrassP' z g2 g3 = weierstrassP z omega1 omega2
  where
    (omega1, omega2) = halfPeriods g2 g3

-- | Derivative of Weierstrass p-function given the half-periods
weierstrassPdash ::
    Complex Double -- ^ z
 -> Complex Double -- ^ half-period omega1
 -> Complex Double -- ^ half-period omega2
 -> Complex Double
weierstrassPdash z omega1 omega2 = 2 / (w1 %^% 3) * j2 * j3 * j4 * f
  where
    w1 = 2 * omega1 / pi
    tau = omega2 / omega1
    q = exp (i_ * pi * tau)
    z' = z / w1 
    j1 = jtheta1' z' tau
    j2 = jtheta2' z' tau
    j3 = jtheta3' z' tau
    j4 = jtheta4' z' tau
    j1dash = jtheta1Dash0 q
    j2zero = jtheta2' 0 tau
    j3zero = jtheta3' 0 tau
    j4zero = jtheta4' 0 tau
    f = j1dash %^% 3 / (j1 %^% 3 * j2zero * j3zero * j4zero)

-- | Derivative of Weierstrass p-function given the elliptic invariants
weierstrassPdash' ::
    Complex Double -- ^ z
 -> Complex Double -- ^ elliptic invariant g2
 -> Complex Double -- ^ elliptic invariant g3
 -> Complex Double
weierstrassPdash' z g2 g3 = weierstrassPdash z omega1 omega2
  where
    (omega1, omega2) = halfPeriods g2 g3

-- | Inverse of Weierstrass p-function given the half-periods
weierstrassPinv ::
    Complex Double -- ^ w
 -> Complex Double -- ^ half-period omega1
 -> Complex Double -- ^ half-period omega2
 -> Complex Double
weierstrassPinv w omega1 omega2 = carlsonRF' 1e-14 (w - e1) (w - e2) (w - e3)
  where
    e1 = weierstrassP omega1 omega1 omega2
    e2 = weierstrassP omega2 omega1 omega2
    e3 = weierstrassP (-omega1 - omega2) omega1 omega2

-- | Inverse of Weierstrass p-function given the elliptic invariants
weierstrassPinv' ::
    Complex Double -- ^ z
 -> Complex Double -- ^ elliptic invariant g2
 -> Complex Double -- ^ elliptic invariant g3
 -> Complex Double
weierstrassPinv' z g2 g3 = weierstrassPinv z omega1 omega2
  where
    (omega1, omega2) = halfPeriods g2 g3

-- | Weierstrass sigma function given the half-periods
weierstrassSigma ::
    Complex Double -- ^ z
 -> Complex Double -- ^ half-period omega1
 -> Complex Double -- ^ half-period omega2
 -> Complex Double
weierstrassSigma z omega1 omega2 = w1 * exp (h * z * z1 / pi) * j1 / j1dash
  where
    tau = omega2 / omega1
    q = exp (i_ * pi * tau)
    w1 = -2 * omega1 / pi
    z1 = z / w1
    j1 = jtheta1' z1 tau
    j1dash = jtheta1Dash0 q
    h = - pi / (6 * w1) * jtheta1DashDashDash0 tau / j1dash

-- | Weierstrass sigma function given the elliptic invariants
weierstrassSigma' ::
    Complex Double -- ^ z
 -> Complex Double -- ^ elliptic invariant g2
 -> Complex Double -- ^ elliptic invariant g3
 -> Complex Double
weierstrassSigma' z g2 g3 = weierstrassSigma z omega1 omega2
  where
    (omega1, omega2) = halfPeriods g2 g3

-- | Weierstrass zeta function given the half-periods
weierstrassZeta ::
    Complex Double -- ^ z
 -> Complex Double -- ^ half-period omega1
 -> Complex Double -- ^ half-period omega2
 -> Complex Double
weierstrassZeta z omega1 omega2 = - eta1 * z + p * lj1dash
  where
    tau = omega2 / omega1
    q = exp (i_ * pi * tau)
    w1 = - omega1 / pi
    p = 0.5 / w1
    j1dash = jtheta1Dash0 q
    eta1 = p * jtheta1DashDashDash0 tau / (6 * w1 * j1dash)
    pz = p * z
    lj1dash = jtheta1Dash pz q / jtheta1' pz tau

-- | Weierstrass zeta function given the elliptic invariants
weierstrassZeta' ::
    Complex Double -- ^ z
 -> Complex Double -- ^ elliptic invariant g2
 -> Complex Double -- ^ elliptic invariant g3
 -> Complex Double
weierstrassZeta' z g2 g3 = weierstrassZeta z omega1 omega2
  where
    (omega1, omega2) = halfPeriods g2 g3
