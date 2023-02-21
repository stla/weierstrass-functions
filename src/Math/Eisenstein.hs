module Math.Eisenstein
    ( eisensteinE4,
      eisensteinE6
    ) where
import           Data.Complex         ( Complex(..) )
import           Math.JacobiTheta     ( jtheta2, jtheta3, jtheta4 )

i_ :: Complex Double
i_ = 0.0 :+ 1.0

-- | Eisenstein series of weight 4
eisensteinE4 :: 
    Complex Double -- ^ tau
 -> Complex Double
eisensteinE4 tau = 
  jtheta2 0 q ** 8 + jtheta3 0 q ** 8 + jtheta4 0 q ** 8 
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

