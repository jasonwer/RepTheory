{-# OPTIONS_GHC -Wall                    #-} 
{-# LANGUAGE OverloadedStrings           #-}

-- |
-- Module      : RepTheory.SuperWeights
-- Copyright   : (c) Jason Werry 2016-2017
-- License     : BSD-style
--
-- Maintainer  : Jason Werry <jason.werry@uqconnect.edu.au>
-- Stability   : experimental
-- Portability : POSIX 
--
-- A library for Lie superalgebra gl(m|n) weights and generation of (super)
-- <https://arxiv.org/abs/math/0211289 Gelfand-Tsetlin> (GT) bases. The core of
-- this library is @getBasis@ which generates a GT basis given a highest gl(m|n) weight.
-- The basis may be exported in TeX format while a Python script to display the
-- 'crystal' structure can be generated.

-- For general information regarding representation theory, weights, roots, bases etc. 
-- see [HUM]. For the super-case see [GIW1] and references therein. 
--
-- References

-- [GIW1] M. D. Gould, P. S. Isaac, and J. L. Werry. Matrix elements for type 1 
-- unitary irreducible representations of the Lie superalgebra gl(m|n). 
-- J. Math. Phys., 55, 011703 (2014).
-- https://arxiv.org/abs/1311.4246
--
-- [HUM] J. E. Humphreys. Introduction to Lie algebras and representation theory
-- volume 9 of Graduate Texts in Mathematics. Springer, New York, 1972.

module RepTheory.SuperWeights (
  
  -- * Construction
  SuperWeight,
  SuperGTPattern,
  SuperBasis,  
{-
  -- * SuperWeight property checking
  zeroWeight,
  isWeaklyDecreasing,
  isStrictDecreasing,
  isDominant,
  hasNonNegLabels,
  
  -- * Algorithms
  getBasis,
  basisWeights,
  weightFromGT,

  -- * Utility functions
  dimensionB,
  dimensionW,
  rhoTimes2,
  getLabel,
  getRootsGL,
  getPositiveRootsGL,
  sumWeights,
  getBetweenness,

  -- * Data exporting

  -- ** TeX
  writeBasis,
  renderGTPattern,
  renderBasis,

  -- ** Python
  writeCrystal,
-}
) where 

import Control.Lens
import Data.Ix
import Data.Text() 
import Data.Text.IO()

-- | A finite dimensional gl(m) highest weight. i.e. A weakly decreasing series of non-negative integers of length m.
type Weight    = [Int]

-- | A gl(m|n) highest weight. 
-- i.e. A pair of weakly decreasing series of non-negative integers of length m and n respecively.

type SuperWeight  = (Weight,Weight)

-- | Gelfand-Tsetlin pattern obeying usual betweenness conditions.
type SuperGTPattern = [SuperWeight] 

-- | A Gelfand-Tsetlin basis of a gl(m|n) module.
type SuperBasis     = [SuperGTPattern]

-- | Vector addition
(+.) :: Weight -> Weight -> Weight
(+.) = zipWith (+)

-- | Even Bosonic part
bose :: SuperWeight -> Weight
bose = fst

-- | Odd Fermionic part
fermi :: SuperWeight -> Weight
fermi = snd

-- | Super vector addition
(+..) :: SuperWeight -> SuperWeight -> SuperWeight 
(+..) a b = ( (fermi a) +. (fermi b) , (bose a) +. (bose b) )

-- | The 1-dimensional zero weight.
zeroWeight :: Int -> Int -> SuperWeight
zeroWeight m n = (replicate m (0::Int) , replicate n (0::Int) ) :: SuperWeight

-- | Component-wise sum of weights
sumWeights :: [SuperWeight] -> SuperWeight
--sumWeights [] = []
sumWeights w = foldl (+..) (zeroWeight m n) w where
  m = length $ fst $ head w
  n = length $ snd $ head w

-- | Test part weight for dominance.
isWeaklyDecreasingPart :: (SuperWeight -> Weight) -> SuperWeight -> Bool
isWeaklyDecreasingPart p s = and $ zipWith (>=) w (tail w) where
  w = p s

-- | Test weight for dominance.
isWeaklyDecreasing :: SuperWeight ->  Bool
isWeaklyDecreasing s = (isWeaklyDecreasingPart bose s) && (isWeaklyDecreasingPart fermi s)

-- | Test weight for strictly decreasing property.
isStrictDecreasingPart :: (SuperWeight -> Weight) -> SuperWeight -> Bool
isStrictDecreasingPart p s = and $ zipWith (>) w (tail w) where
  w = p s

-- | Aliased function to match RepTheory terminology (=isWeaklyDecreasing).
isDominant :: SuperWeight -> Bool
isDominant = isWeaklyDecreasing  

-- | Test the set of weight labels for non-negativity.
hasNonNegLabels :: SuperWeight -> Bool
hasNonNegLabels s = (all (>= 0) $ bose s) && (all (>= 0) $ fermi s)

-- | Generate a list of possible values that obey the gl(m|n)
-- betweenness conditions. See [GIW1], Theorem 7.
-- e.g All numbers between 3 and 0
-- getBetweenness [3,0]
-- [[3,2,1,0]]
-- e.g All numbers between 5 and 4, 4 and 2, 2 and 0
-- [[5,4],[4,3,2],[2,1,0]]
getBetweenness :: SuperWeight -> [Weight]
getBetweenness w = evens ++ odds where
  evens = map ( \(x,y) -> (reverse $ range(y,x) ) ) $ zip t $  (map (subtract 1) t)
  odds = map ( \(x,y) -> (reverse $ range(y,x) ) ) $ zip s $ tail s
  t = bose w
  s = fermi w


-- | Generate all gl(m|n-1) highest SuperWeights given a gl(m|n) highest SuperWeight. 
-- Numerically, we generate all combinations of the possibilities obtained 
-- from the betweenness conditions. See [GIW1], Theorem 7.
getSubModulesHW :: SuperWeight -> [ SuperWeight ]
getSubModulesHW w = map (\t -> (splitAt m t)  ) $ sequence $ getBetweenness w where
  m = length $ bose w

-- References

-- [GIW1] M. D. Gould, P. S. Isaac, and J. L. Werry. Matrix elements for type 1 
-- unitary irreducible representations of the Lie superalgebra gl(m|n). 
-- J. Math. Phys., 55, 011703 (2014).
-- https://arxiv.org/abs/1311.4246

-- [WGI1] J. L. Werry, M. D. Gould, and P. S. Isaac. Matrix elements and duality for type 2 
-- unitary irreducible representations of the Lie superalgebra gl(m|n). 
-- J. Math. Phys., 56, 121703 (2015).
-- https://arxiv.org/abs/1506.07593

-- [HUM] J. E. Humphreys. Introduction to Lie algebras and representation theory
-- volume 9 of Graduate Texts in Mathematics. Springer, New York, 1972.


---------------------------------------------

