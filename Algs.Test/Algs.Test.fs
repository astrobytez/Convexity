module Algs.Test

open NUnit.Framework
open MathNet.Numerics.LinearAlgebra
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
    [<Test>]
    let TestDmdSimple () =

        // Create stable spring mass damper system
        let numPoints = 10
        let dt = 0.01<Dmd.sec>

        let m = 1.0
        let k = 1.0
        let c = 1.0
        let A = [ [ 0.0; 1.0 ]; [ -k/m; -c/m ] ] |> matrix

        // Initial condition
        let u0 = [0.0; 5.0] |> vector

        // Propogate u0 through system A
        let folder u0 (t: float) =
            u0 @ [ (matrixExp (A * t)) * (List.head u0) ]

        let ts = [ 0.0 .. float dt .. float numPoints ]
        let trajectory =
            ts
            |> List.fold folder [ u0 ]
            |> CreateMatrix.DenseOfColumnVectors

        // Setup the problem and compute model
        let (X0, X1) = splitData trajectory
        let model = Dmd.exactDmd dt X0 X1

        let folder' acc _ =
            acc @ [ Dmd.predictNext dt (List.last acc) model.Phi model.Omega ]

        let trajectory' =
            ts
            |> List.fold folder' [ u0 ]
            |> CreateMatrix.DenseOfColumnVectors

        // Chart.show <| plot ts trajectory
        // Chart.show <| plot ts trajectory'

        // Assert the reconstructed trajectory closely matches the baseline.
        let err = trajectory - trajectory'
        for row in err.EnumerateColumns() do
            Assert.That(row.L2Norm(), Is.LessThan(1e-1))
