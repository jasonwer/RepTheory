{-# OPTIONS_GHC -Wall                    #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing #-}
{-# OPTIONS_GHC -fno-warn-orphans        #-}
{-# OPTIONS_GHC -fdefer-type-errors      #-} 

{-# LANGUAGE OverloadedStrings           #-}

-- |
-- Module      : RepTheory.Weights
-- Copyright   : (c) Jason Werry 2016-2017
-- License     : BSD-style
--
-- Maintainer  : Jason Werry <jason.werry@uqconnect.edu.au>
-- Stability   : experimental
-- Portability : POSIX 
--
-- A library for gl(m) weights (in the representation theoretic sense) and generation of 
-- <https://arxiv.org/abs/math/0211289 Gelfand-Tsetlin> (GT) bases. The core of
-- this library is @getBasis@ which generates a GT basis given a highest weight.
-- The basis may be exported in TeX format while a Python script to display the
-- 'crystal' structure can be generated.
--
-- References
-- [Hum] J. E. Humphreys. Introduction to Lie algebras and representation theory
-- volume 9 of Graduate Texts in Mathematics. Springer, New York, 1972.

module RepTheory.Weights (
  
  -- * Construction
  Weight,
  GTPattern,
  Basis,  

  -- * Weight property checking
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
) where 

import Control.Lens
import Data.Ix
import Data.Text() 
import Data.Text.IO()

-- | A finite dimensional gl(m) highest weight. i.e. A weakly decreasing series of non-negative integers of length m.
type Weight    = [Int]

-- | Gelfand-Tsetlin pattern obeying usual betweenness conditions.
type GTPattern = [Weight] 

-- | A Gelfand-Tsetlin basis of a gl(m) module.
type Basis     = [GTPattern]

--(*.) :: Weight -> Weight -> Int
--(*.) = (sum .) . zipWith (*)

-- | Vector addition
(+.) :: Weight -> Weight -> Weight 
(+.) = zipWith (+)

-- | The 1-dimensional zero weight.
zeroWeight :: Int -> Weight
zeroWeight m = replicate m (0::Int) :: Weight

-- | Component-wise sum of weights
sumWeights :: [Weight] -> Weight
sumWeights [] = []
sumWeights w = foldl (+.) (zeroWeight m) w where
  m = length $ head w

-- | Test weight for dominance.
isWeaklyDecreasing :: Weight -> Bool
isWeaklyDecreasing s = and $ zipWith (>=) s (tail s)

-- | Test weight for strictly decreasing property.
isStrictDecreasing :: Weight -> Bool
isStrictDecreasing s = and $ zipWith (>) s (tail s)

-- | Aliased function to match RepTheory terminology (=isWeaklyDecreasing).
isDominant :: Weight -> Bool
isDominant = isWeaklyDecreasing  

-- | Test the set of weight labels for non-negativity.
hasNonNegLabels :: Weight -> Bool
hasNonNegLabels = all (>= 0) 

------------------------------

-- | Generate a list of possible values that obey the standard gl(m)
-- betweenness conditions.
-- e.g All numbers between 3 and 0
-- getBetweenness [3,0]
-- [[3,2,1,0]]
-- e.g All numbers between 5 and 4, 4 and 2, 2 and 0
-- [[5,4],[4,3,2],[2,1,0]]
getBetweenness :: Weight -> [Weight]
getBetweenness s = map ( \(x,y) -> (reverse $ range(y,x) ) :: Weight) $ zip s $ tail s

-- | Generate all gl(m-1) highest weights given a gl(m) highest weight. 
-- Numerically, we generate all combinations of the possibilities obtained 
-- from the betweenness conditions.
getSubModulesHW :: Weight -> [ Weight ]
getSubModulesHW = sequence.getBetweenness

-- | Given a partial GT pattern (containing rows p to p-k > 1) generate the
-- set of GT patterns with the (p-k)'th row obeying the betweenness conditions.
generateSubModules :: GTPattern -> [GTPattern]
generateSubModules g = map (\x -> g ++ [x] :: GTPattern) s
  where s = getSubModulesHW $ last g

-- | Recursively generate GT pattern rows from the top down.
generateBasis :: Basis -> Basis
generateBasis b 
 | length (last $ head b) == 1 = b
 | otherwise = generateBasis $ foldl1 (++) $ map generateSubModules b

-- | Generate a GT basis (a set of GT patterns) given a highest weight. 
getBasis :: Weight -> Basis
getBasis wt = generateBasis b
  where b   = [gt1] 
        gt1 = [wt] 

-- | Obtain i'th weight label from p'th level
--   where the 1st level is the top row of the GT pattern
--   and the labels of each row are indexed from 0.
getLabel :: GTPattern -> Int -> Int -> Int
getLabel g p i = g !! (length g - p ) !! (i-1)
------------------------------------- 
-- | Dimension of a gl(m) highest weight module
dimensionW :: Weight -> Int
dimensionW w = dimensionB $ getBasis w

-- | Dimension of a gl(m) GT basis
dimensionB :: Basis -> Int
dimensionB = length 

-- | Weight of a GT pattern (set of neighbouring row-sum differences)
weightFromGT :: GTPattern -> Weight
weightFromGT gt = reverse (v2 ++ [last v]) -- so that zipWith (-) [4,3,0] [3,0] gives [1,3,0] (not [1,3])
  where v2 = zipWith (-) v (tail v)  -- subtract the sum of each GT row from the one above (left)
        v  = map sum gt            -- sum each row of the GT pattern

-- | Get the spectrum of weights from the basis 
basisWeights :: Basis -> [Weight]
basisWeights = map weightFromGT

-- | 2*rho = sum of the positive roots (See [HUM]) 
rhoTimes2 :: Int -> Weight
rhoTimes2 m = sumWeights $ getPositiveRootsGL m

-- | Get all roots of gl(m) - see [HUM]
getRootsGL :: Int -> [Weight]
getRootsGL m = map (\(x,y) -> zeroWeight m & element x .~ 1 & element y .~ (-1) ) pairs where
  pairs = [(i::Int,j::Int) | i <- [0..m-1] , j <- [0..m-1], i /= j ]

-- | Get all positive roots of gl(m) - see [HUM]. 
getPositiveRootsGL :: Int -> [Weight]
getPositiveRootsGL m = map (\(x,y) -> zeroWeight m & element x .~ 1 & element y .~ (-1) ) pairs where
  pairs = [(i,j) | i <- [0..m-1] , j <- [0..m-1], i < j ]
-------------------------------------------------------------

weightToStrings :: Weight -> [String]
weightToStrings = map show

padWeight :: [String] -> Int -> [String]
padWeight w l = w ++ replicate (l - length w) " "

columnDelimit :: [String] -> Int -> String
columnDelimit w l = concat (zipWith (++) (padWeight w l) (replicate l "&")) ++ "\\\\ \n"

---------------- TeX document file rendering--------------------

-- | Writes a TeX file containing the GT basis patterns.
writeBasis :: Basis 
             -> String -- ^ FilePath 
             -> IO()
writeBasis b f = writeFile f (texHeader ++ renderBasis b ++ texFooter)

-- e.g input of [ [4,2,1,0],[4,1,1],[3,1],[3] ]
-- returns the TeX 
{-
\left\lvert
\begin{matrix}
 4&2&1&0& \\
 4&1&1& & \\
 3&1& & & \\ 
 3& & & & \\
\end{matrix}
\right)
-}
renderGTPattern :: GTPattern -> String
renderGTPattern g = startDelim ++ startMat ++ patRows ++ endMat ++ endDelim where
  startDelim = "\\left\\lvert"
  startMat = "\\begin{matrix}\n"
  patRows = concatMap (\x -> columnDelimit (weightToStrings x) padWidth) g
  endMat = "\\end{matrix}"
  endDelim = "\\right),\\\\ \n"
  padWidth = length $ head g

renderBasis :: Basis -> String
renderBasis b = "\\begin{align}\n" ++ concatMap renderGTPattern b ++ "\\end{align}\n"

texHeader :: String
texHeader = "\\documentclass[12pt]{article} \n\\usepackage{amsmath}\n\n\\begin{document}\n"

texFooter :: String
texFooter = "\n\\end{document}"

---------------- Python graphics file rendering--------------------

-- | Writes a Python script that generates a graphical representation of the weight spectrum.
writeCrystal :: Basis
             -> String -- ^ FilePath
             -> IO()  
writeCrystal b f = writeFile f (pyHeader ++ renderNodes b ++ renderEdges ++ renderPositions ++ pyFooter) where 
  renderEdges = "\n todo - renderEdges \n"
  renderPositions = "\n todo - renderPositions \n "                

renderNodes :: Basis -> String
renderNodes = concatMap (\x -> "G.add_node('" ++ show x ++ "') \n" )  

--renderEdges :: Basis -> String
--renderPositions :: Basis -> String

pyHeader :: String
pyHeader = "# -*- coding: utf-8 -*- \n\
\import networkx as nx \n\
\import matplotlib.pyplot as plt \n\
\import math \n\
\ \n\
\G = nx.DiGraph() \n\
\ \n"

pyFooter :: String
pyFooter = "\nfixed_nodes = fixed_positions.keys() \n\
\pos = nx.spring_layout(G,pos=fixed_positions, fixed = fixed_nodes, iterations=5)\n\
\    \n\
\fig = plt.gcf() \n\
\fig.set_size_inches(12.5, 8.5, forward=True) \n\
\fig.savefig('crystal.png', dpi=100) \n\
\\n\
\nx.draw_networkx(G,pos,alpha=0.5,font_size=12) \n"


-- References

-- [Hum] J. E. Humphreys. Introduction to Lie algebras and representation theory
-- volume 9 of Graduate Texts in Mathematics. Springer, New York, 1972.


---------------------------------------------



