module Algs.Test

open NUnit.Framework
open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.LinearAlgebra.Double
open Plotly.NET

open Algs.Core


[<TestFixture>]
module TestDynamicModeDecomposition =

  let inline splitData (X: Matrix<'a>) =
    let m = X.ColumnCount
    X[*, 0 .. m - 2], X[*, 1 .. m - 1] // Zero indices

  let matrixExp (X: Matrix<float>) =

    let evd = X.ToComplex().Evd()
    let Q = evd.EigenValues.PointwiseExp() |> CreateMatrix.DenseOfDiagonalVector
    let Z = evd.EigenVectors * Q * evd.EigenVectors.Inverse()
    Z.Real()

  let plot ts (t: Matrix<float>) =
    seq { for row in t.EnumerateRows() -> Chart.Line(x = ts, y = row) }
    |> Chart.combine

  /// A test to show the Dmd model is captured from the data.
  /// Show the model can be used to reconstruct the signal within tolerance.
  [<TestCase>]
  let ``Test DMD algorithm Simple`` () =

    // Create stable spring mass damper system
    let numPoints = 10
    let dt = 0.01<Dmd.sec>

    let m = 1.0
    let k = 1.0
    let c = 1.0
    let A = [ [ 0.0; 1.0 ]; [ -k / m; -c / m ] ] |> matrix

    // Initial condition
    let u0 = [ 0.0; 5.0 ] |> vector

    // Propogate u0 through system A
    let folder u0 (t: float) =
      u0 @ [ (matrixExp (A * t)) * (List.head u0) ]

    let ts = [ 0.0 .. float dt .. float numPoints ]
    let trajectory = ts |> List.fold folder [ u0 ] |> CreateMatrix.DenseOfColumnVectors

    // Setup the problem and compute model
    let (X0, X1) = splitData trajectory
    let model = Dmd.exactDmd dt X0 X1

    let folder' acc _ =
      acc @ [ Dmd.predictNext dt (List.last acc) model.Phi model.Omega ]

    let trajectory' =
      ts |> List.fold folder' [ u0 ] |> CreateMatrix.DenseOfColumnVectors

    // Assert the reconstructed trajectory closely matches the baseline.
    let err = trajectory - trajectory'

    for row in err.EnumerateColumns() do
      Assert.That(row.L2Norm(), Is.LessThan(1e-1))

  [<TestCase(0, 10, 0.1, 0.5, 0.050, ExpectedResult = true)>]
  [<TestCase(1, 10, 0.1, 0.5, 0.050, ExpectedResult = true)>]
  [<TestCase(3, 10, 0.1, 0.5, 0.050, ExpectedResult = true)>]
  [<TestCase(4, 10, 0.1, 0.5, 0.075, ExpectedResult = true)>]
  [<TestCase(5, 10, 0.1, 0.5, 0.075, ExpectedResult = true)>]
  let ``Test Sequential Threshold Least Squares Norm Error`` (seed, numIters, (gain: float), lambda, tolerance) =
    let x = [ 1.23; 0.0; 4.56; 0.0; 7.89 ] |> vector // 3 active terms
    let noise = gain * DenseVector.randomSeed 5 seed
    let A = DenseMatrix.randomSeed 1000 5 seed |> DenseMatrix.OfMatrix

    let bNoisy = (A * (x + noise)).ToColumnMatrix() |> DenseMatrix.OfMatrix

    let xTilde = Solvers.stlsq numIters lambda A bNoisy
    let errors = xTilde - x.ToColumnMatrix()
    let errorNorm = errors.TransposeThisAndMultiply errors
    errorNorm[0, 0] < tolerance

  [<TestCase(0, 10, 0.10, 0.5, ExpectedResult = 3)>]
  [<TestCase(1, 10, 0.20, 0.5, ExpectedResult = 3)>]
  [<TestCase(3, 10, 0.25, 0.5, ExpectedResult = 3)>]
  [<TestCase(4, 10, 0.20, 0.5, ExpectedResult = 3)>]
  [<TestCase(5, 10, 0.10, 0.5, ExpectedResult = 3)>]
  let ``Test Sequential Threshold Least Squares Term Count`` (seed, numIters, (gain: float), lambda) =
    let x = [ 1.23; 0.0; 4.56; 0.0; 7.89 ] |> vector // 3 active terms
    let noise = gain * DenseVector.randomSeed 5 seed
    let A = DenseMatrix.randomSeed 1000 5 seed |> DenseMatrix.OfMatrix

    let bNoisy = (A * (x + noise)).ToColumnMatrix() |> DenseMatrix.OfMatrix

    let xTilde = Solvers.stlsq numIters lambda A bNoisy
    let numTerms = xTilde[*, 0].Enumerate() |> Seq.countBy (fun x -> x > 0) |> Map.ofSeq

    numTerms[true] // Should be the number of active terms in the vector x

  [<TestCase(0, 100, 0.2, 0.100, 0.50, 0.25, 1e-6, 1e-3, 0.50, ExpectedResult = true)>]
  [<TestCase(1, 100, 0.2, 0.100, 0.50, 0.25, 1e-6, 1e-3, 0.50, ExpectedResult = true)>]
  [<TestCase(3, 100, 0.2, 0.075, 0.05, 0.25, 1e-6, 1e-3, 0.25, ExpectedResult = true)>]
  [<TestCase(4, 100, 0.2, 0.075, 0.05, 0.25, 1e-6, 1e-3, 0.25, ExpectedResult = true)>]
  [<TestCase(5, 100, 0.2, 0.075, 0.05, 0.25, 1e-6, 1e-3, 0.25, ExpectedResult = true)>]
  let ``Test Least Absolute Shrinkage and Selection Operator Norm Error``
    (seed, numIters, (gain: float), lambda, rho, alpha, absTol, relTol, tolerance)
    =
    let x = [ 1.23; 0.0; 4.56; 0.0; 7.89 ] |> vector // 3 active terms
    let noise = gain * DenseVector.randomSeed 5 seed
    let A = DenseMatrix.randomSeed 1000 5 seed |> DenseMatrix.OfMatrix

    let bNoisy = (A * (x + noise))

    let lasso = Solvers.lasso lambda rho alpha absTol relTol numIters

    let (xTilde, _, _, _) = lasso A bNoisy
    let errors = xTilde - x
    let errorNorm = errors * errors
    errorNorm < tolerance

  [<TestCase(0, 100, 0.2, 0.150, 0.50, 0.50, 1e-6, 1e-3, ExpectedResult = 5)>]
  [<TestCase(1, 100, 0.2, 0.100, 0.50, 0.25, 1e-6, 1e-3, ExpectedResult = 4)>]
  [<TestCase(3, 100, 0.2, 0.075, 0.05, 0.25, 1e-6, 1e-3, ExpectedResult = 4)>]
  [<TestCase(4, 100, 0.2, 0.075, 0.05, 0.25, 1e-6, 1e-3, ExpectedResult = 3)>]
  [<TestCase(5, 100, 0.2, 0.075, 0.05, 0.25, 1e-6, 1e-3, ExpectedResult = 4)>]
  let ``Test Least Absolute Shrinkage and Selection Operator Term Count``
    (seed, numIters, (gain: float), lambda, rho, alpha, absTol, relTol)
    =
    (*
            Note - LASSO algorithm consistently performs poorer than STLS in finding the correct parsimony.
            There is consistently more active terms detected in the system than actually present.

            The general conclusion is not prefer STLS over LASSO.
        *)

    let x = [ 1.23; 0.0; 4.56; 0.0; 7.89 ] |> vector // 3 active terms
    let noise = gain * DenseVector.randomSeed 5 seed
    let A = DenseMatrix.randomSeed 1000 5 seed |> DenseMatrix.OfMatrix

    let bNoisy = (A * (x + noise))

    let lasso = Solvers.lasso lambda rho alpha absTol relTol numIters

    let (xTilde, _, _, _) = lasso A bNoisy
    let numTerms = xTilde.Enumerate() |> Seq.countBy (fun x -> x > 0) |> Map.ofSeq

    numTerms[true] // Should be the number of active terms in the vector x