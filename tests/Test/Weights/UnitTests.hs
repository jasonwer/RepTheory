module Test.Weights.UnitTests where

  import RepTheory.Weights
  import Test.QuickCheck

  runTests :: IO()
  runTests = do
    putStrLn "Running tests for module RepTheory.Weights"
    putStrLn "  Weight functions"
    quickCheck (zeroWeight 5 == [0,0,0,0,0])
    quickCheck (sumWeights [[1,2,3,4],[4,3,2,1],[1,1,1,1]] == [6,6,6,6])

    putStrLn "  Betweenness"
    quickCheck (getBetweenness [5,4,2,0] == [[5,4],[4,3,2],[2,1,0]])

    putStrLn "  GT Basis of (x,0,0,..)"
    quickCheck (getBasis [0] == [[[0]]])

    quickCheck (getBasis [0,0,0] == [[[0,0,0],[0,0],[0]]])

    quickCheck (getBasis [1,0,0] == [[[1,0,0],[1,0],[1]],[[1,0,0],[1,0],[0]],[[1,0,0],[0,0],[0]]])

    quickCheck (getBasis [2,0,0] == [[[2,0,0],[2,0],[2]],[[2,0,0],[2,0],[1]],[[2,0,0],[2,0],[0]],[[2,0,0],[1,0],[1]],[[2,0,0],[1,0],[0]],[[2,0,0],[0,0],[0]]])

    putStrLn "  GT Basis of (3,1,0)"
    quickCheck (getBasis [3,1,0] == [[[3,1,0],[3,1],[3]],[[3,1,0],[3,1],[2]],[[3,1,0],[3,1],[1]],[[3,1,0],[3,0],[3]],[[3,1,0],[3,0],[2]],[[3,1,0],[3,0],[1]],[[3,1,0],[3,0],[0]],[[3,1,0],[2,1],[2]],[[3,1,0],[2,1],[1]],[[3,1,0],[2,0],[2]],[[3,1,0],[2,0],[1]],[[3,1,0],[2,0],[0]],[[3,1,0],[1,1],[1]],[[3,1,0],[1,0],[1]],[[3,1,0],[1,0],[0]]])

    putStrLn "  GT Basis of other highest weights"
    quickCheck (weightFromGT [[2,1,0],[2,1],[2]] == [2,1,0])
    quickCheck (weightFromGT [[2,1,0],[2,1],[0]] == [0,3,0])
    quickCheck (weightFromGT [[4,2,1,0],[3,2,0],[2,1],[1]] == [1,2,2,2])

    putStrLn "  Calculate half-sum of postive roots"
    quickCheck (rhoTimes2 6 == [5,3,1,-1,-3,-5])

    putStrLn "  Dimensions of trivial reps"
    quickCheck (dimensionW [0] == 1)
    quickCheck (dimensionW [0,0,0,0,0] == 1)
    quickCheck (dimensionW [5,5,5,5,5] == 1)

    putStrLn "  Dimensions of vector reps"
    quickCheck (dimensionW [1,0,0] == 3)
    quickCheck (dimensionW [1,0,0,0] == 4)
    quickCheck (dimensionW [1,0,0,0,0] ==5)

    putStrLn "  Dimensions of tensor reps"
    quickCheck (dimensionW [2,0,0] == 6)
    quickCheck (dimensionW [3,0,0] == 10)
    quickCheck (dimensionW [4,0,0] == 15)

    putStrLn "  Large dimensional reps"
    quickCheck (dimensionW [4,3,2,1,0] == 1024)
    quickCheck (dimensionW [5,3,3,1,0] == 1890)

