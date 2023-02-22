module Internal ((%^%)) where
import Data.Complex ( Complex (..) )

(%^%) :: Complex Double -> Int -> Complex Double
(%^%) z p = z ^^ p
