module Math.Eisenstein
    ( lambda,
      eisensteinE2,
      eisensteinE4,
      eisensteinE6,
      kleinJ,
      kleinJinv,
      modularDiscriminant,
      agm,
      etaDedekind,
      jtheta1DashDashDash0
    ) where
import           Data.Complex           ( Complex(..) )
import           Internal               ( (%^%) )
import           Math.EllipticIntegrals ( ellipticF', ellipticE' )
import           Math.JacobiTheta       ( jtheta2, jtheta3, jtheta4 )


i_ :: Complex Double
i_ = 0.0 :+ 1.0

-- | Lambda modular function (square of elliptic modulus)
lambda :: 
    Complex Double -- ^ tau
 -> Complex Double
lambda tau = (j2 / j3) %^% 4
    where
      q = exp (i_ * pi * tau)
      j2 = jtheta2 0 q
      j3 = jtheta3 0 q

-- | Eisenstein series of weight 2
eisensteinE2 :: 
    Complex Double -- ^ tau
 -> Complex Double
eisensteinE2 tau = 
  6 / pi * ellE * j3 - j3 * j3 - j4
    where
      q = exp (i_ * pi * tau)
      j3 = jtheta3 0 q %^% 2
      j4 = jtheta3 0 q %^% 4
      ellE = ellipticE' 1e-14 (pi/2) (lambda tau)

-- | Eisenstein series of weight 4
eisensteinE4 :: 
    Complex Double -- ^ tau
 -> Complex Double
eisensteinE4 tau = 
  (jtheta2 0 q %^% 8 + jtheta3 0 q %^% 8 + jtheta4 0 q %^% 8) / 2 
    where
      q = exp (i_ * pi * tau)

-- | Eisenstein series of weight 6
eisensteinE6 :: 
    Complex Double -- ^ tau
 -> Complex Double
eisensteinE6 tau = 
  (jtheta3 0 q %^% 12 + jtheta4 0 q %^% 12 - 3 * jtheta2 0 q %^% 8 
    * (jtheta3 0 q %^% 4 + jtheta4 0 q %^% 4)) / 2
    where
      q = exp (i_ * pi * tau)

-- | Modular discriminant
modularDiscriminant ::
    Complex Double -- ^ tau
 -> Complex Double
modularDiscriminant tau = 
  (eisensteinE4 tau %^% 3 - eisensteinE6 tau %^% 2) / 1728

-- | Klein J-function
kleinJ :: 
    Complex Double -- ^ tau
 -> Complex Double
kleinJ tau = 
  eisensteinE4 tau %^% 3 / modularDiscriminant tau

-- | Arithmetic-geometric mean
agm :: 
    Complex Double -- ^ x 
 -> Complex Double -- ^ y
 -> Complex Double
agm x y = 
  if x == 0 || y == 0 || x + y == 0
    then 0
    else pi/4 * (x + y) / ellipticF' 1e-15 (pi/2) (((x - y) / (x + y)) %^% 2)

-- | Inverse Klein-J function
kleinJinv :: 
    Complex Double
 -> Complex Double
kleinJinv j = 
  if j == 0
    then 0.5 :+ (sqrt 3 / 2)
    else i_ * agm 1 (sqrt(1 - lbd)) / agm 1 (sqrt lbd)
  where
    j2 = j * j
    j3 = j2 * j
    t = -j3 + 2304 * j2 + 12288 * sqrt(3 * (1728 * j2 - j3)) - 884736 * j
    u = t ** (1/3)
    x = 4 + (u - j) / 192 - (1536 * j - j2) / (192 * u)
    lbd = -(-1 - sqrt(1 - x)) / 2

-- | Dedekind eta function
etaDedekind ::
    Complex Double -- ^ tau
 -> Complex Double
etaDedekind tau = exp (ipitau / 12) * j3
  where
    ipitau = i_ * pi * tau
    q = exp (3 * ipitau)
    j3 = jtheta3 (pi / 2 * (tau + 1)) q

-- | Third derivative at 0 of the first Jacobi theta function
jtheta1DashDashDash0 :: 
    Complex Double -- ^ tau
 -> Complex Double
jtheta1DashDashDash0 tau = -2 * etaDedekind tau %^% 3 * eisensteinE2 tau 


