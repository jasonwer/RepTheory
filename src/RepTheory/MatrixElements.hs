{-# OPTIONS_GHC -Wall                    #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing #-}
{-# OPTIONS_GHC -fno-warn-orphans        #-}
{-# OPTIONS_GHC -fdefer-type-errors      #-} 

{-# LANGUAGE OverloadedStrings           #-}

-- |
-- Module      : RepTheory.MatrixElements
-- Copyright   : (c) Jason Werry 2016-2017
-- License     : BSD-style
--
-- Maintainer  : Jason Werry <jason.werry@uqconnect.edu.au>
-- Stability   : experimental
-- Portability : POSIX 
--
-- A library for calculating matrix elements and Wigner coefficients of gl(m) modules.

module RepTheory.MatrixElements (
  getMatrixElement,
  diffLabels,    -- :: GTPattern -> Level1 -> Level2 -> Offset -> [Int]
) where 

import RepTheory.Weights

-- |Differences of basis labels.
--Returns a list of terms (evaluated to integers) that give
--the differences of the basis labels added to the difference
--of the indexes themselves plus an offset.
--
--Example: For a GTPattern [ [4,2,1,0], [4,2,0] ,[3,2], [2] ]
--the first sublist [4,2,1,0] has subalgebra level 4, and so
--on to [2] which has subalgebra level 1. 

diffLabels :: GTPattern -- ^ a complete GT pattern containing subalgebra levels m, m-1, ..., 1
           -> Int       -- ^ The subalgebra level 'm' of the fixed label
           -> Int       -- ^ The fixed label 'r' (1-based)
           -> Int       -- ^ The 2nd subalgebra level (if ==m then the i=j term is skipped)
           -> Int       -- ^ Offset of each term (usually -1,0, or 1)
           -> [Int]     -- ^ Returned list of label differences
diffLabels g m r p offset =  filter (/=0) $ map (\i -> fixed - row!!(i-1) + i - r + offset) [1..p]
  where row = g!!(length g - p)
        fixed = g!!(length g - m)!!(r-1)

-- | Returns the matrix element N^m_r given a subalgebra level m and position r of the weight to be shifted.
getMatrixElement :: GTPattern
           -> Int         -- ^ The subalgebra level 'm' of the fixed label
           -> Int         -- ^ The fixed label 'r' (1-based)
           -> Float       -- ^ Returned matrix element
getMatrixElement g m r = fromIntegral ((product p1) * (product p2)) /  fromIntegral ((product p3) * (product p4))
  where p1 = diffLabels g m r (m+1) 0
        p2 = diffLabels g m r (m-1) 1
        p3 = diffLabels g m r m     0
        p4 = diffLabels g m r m     1 
        










