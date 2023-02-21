module Math.Eisenstein
    ( eisensteinE4,
      eisensteinE6,
      kleinJ,
      modularDiscriminant,
      agm
    ) where
import           Data.Complex           ( Complex(..) )
import           Math.EllipticIntegrals ( ellipticF' )
import           Math.JacobiTheta       ( jtheta2, jtheta3, jtheta4 )


i_ :: Complex Double
i_ = 0.0 :+ 1.0

-- | Eisenstein series of weight 4
eisensteinE4 :: 
    Complex Double -- ^ tau
 -> Complex Double
eisensteinE4 tau = 
  (jtheta2 0 q ** 8 + jtheta3 0 q ** 8 + jtheta4 0 q ** 8) / 2 
    where
      q = exp (i_ * pi * tau)

-- | Eisenstein series of weight 6
eisensteinE6 :: 
    Complex Double -- ^ tau
 -> Complex Double
eisensteinE6 tau = 
  (jtheta3 0 q ** 12 + jtheta4 0 q ** 12 - 3 * jtheta2 0 q ** 8 
    * (jtheta3 0 q ** 4 + jtheta4 0 q ** 4)) / 2
    where
      q = exp (i_ * pi * tau)

-- | Modular discriminant
modularDiscriminant ::
    Complex Double -- ^ tau
 -> Complex Double
modularDiscriminant tau = 
  (eisensteinE4 tau ** 3 - eisensteinE6 tau ** 2) / 1728

-- | Klein J-function
kleinJ :: 
    Complex Double -- ^ tau
 -> Complex Double
kleinJ tau = 
  eisensteinE4 tau ** 3 / modularDiscriminant tau

-- | Arithmetic-geometric mean
agm :: 
    Complex Double -- ^ x 
 -> Complex Double -- ^ y
 -> Complex Double
agm x y = 
  if x == 0 || y == 0 || x + y == 0
    then 0
    else pi/4 * (x + y) / ellipticF' 1e-15 (pi/2) (((x-y)/(x+y))**2)

