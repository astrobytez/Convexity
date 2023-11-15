namespace Algs.Core

open System.Numerics
open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.LinearAlgebra.Complex

module internal Helpers =

    let (|IsEven|IsOdd|) n =
        match n % 2 with
        | 0 -> IsEven n
        | _ -> IsOdd n

    // Compute the median value from a set of values.
    let median v =
        let sorted = v |> Array.sort
        let length: int = sorted.Length

        match length with
        | x when x < 1 -> 0.0
        | _ ->
            match length with
            | IsEven length ->
                let next = length / 2
                let prev = length / 2 - 1
                (sorted.[prev] + sorted.[next]) / 2.0
            | IsOdd length ->
                let index = int <| ceil (float length) / 2.0
                sorted.[index]

    /// Compute the optimal singular value threshold for the system.
    /// Adapted from https://arxiv.org/pdf/1305.5870.pdf
    let singularValueThresholding (M: Matrix<float>) =
        let rc = (float M.RowCount)
        let cc = (float M.ColumnCount)
        let beta = cc / rc

        let y = M.Diagonal().ToArray() |> median

        let omega = 0.56 * beta ** 3 - 0.95 * beta ** 2 + 1.82 * beta + 1.43

        y * omega // Compute the optimal threshold

[<RequireQualifiedAccessAttribute>]
module Dmd =
    open Helpers

    [<Measure>]
    type sec // A type to represent time units in seconds


    type DmdModel =
        { Phi: DenseMatrix; Omega: DenseVector }


    /// Compute the Dmd Modes for the X1 = A.X0 dynamical system.
    /// Uses the SVD approach with truncation of singular values to optimal set
    /// based on the singular value thresholding algorithm above.
    /// Adapted from https://cwrowley.princeton.edu/theses/tu.pdf
    ///
    /// seriesDeltaTime: A number of time spanning each snapshot in seconds.
    /// X0: matrix containing evenly spaced time series snapshots 0 .. m-1
    /// X1: matrix containing evenly spaced time series snapshots 1 .. m-0
    ///
    /// New papers show interesting results for improvements to the algorithm:
    /// https://doi.org/10.3390/computation10120210
    let exactDmd (seriesDeltaTime: float<sec>) (X0: Matrix<float>) (X1: Matrix<float>) =

        let minModes = 2
        let maxModes = 20

        let (U, S, VT) =
            let svd = X0.Svd(true)
            svd.U, svd.S, svd.VT

        let r =
            let tau =
                // Compute truncated singular values r from SVT algorithm
                DenseMatrix.ofDiag S |> singularValueThresholding

            min
                maxModes
                (seq {
                    for s in S do
                        if s > tau then
                            yield s
                 }
                 |> Seq.length)
            |> max minModes

        let UrT = U[0 .. (r - 1), *]
        let Sr = DenseMatrix.ofDiag S[0 .. (r - 1)]
        let Vr = VT[0 .. (r - 1), *]
        let SinV = Vr.ConjugateTranspose() * Sr.Inverse() // Review - this is expensive with inverse and mutliply.
        let X1SinV = X1 * SinV
        let Atilde = UrT * X1SinV

        let X1SinV = X1SinV.ToComplex()
        let Atilde = Atilde.ToComplex()

        /// Compute Eigen Decomposition for A
        let Ev = Atilde.Evd()
        let Phi = X1SinV * Ev.EigenVectors |> DenseMatrix.OfMatrix // Review - this copy is inefficient

        let dt = float seriesDeltaTime

        let Omega =
            Ev.EigenValues.PointwiseLog().Multiply(Complex(1.0 / dt, 0.0))
            |> DenseVector.OfVector // Review - this copy is inefficient

        { Phi = Phi; Omega = Omega }

    /// Take the DMD eigen system and project forward one discrete step.
    /// Extrapolate a given Dmd model.
    ///
    /// seriesDeltaTime: A number of time spanning each snapshot in seconds.
    /// u0: A initial condition vector to extrapolate the model forward in time from.
    /// Phi: The Dmd modes matrix.
    /// The Dmd eigen values vector.
    ///
    /// Uses a naive approach to compute the modes by solving directly.
    /// There are more efficient methods to compute this as well as methods
    /// which better extract the dominant modes, e.g. sparsity promoting.
    let predictNext (seriesDeltaTime: float<sec>) (u0: Vector<float>) (phi: Matrix<Complex>) (omega: Vector<Complex>) =
        let dt = seriesDeltaTime
        let b = phi.Solve(u0.ToComplex()) // Todo - this is inefficient there are better solutions for computing the modes.

        let A =
            omega
            |> Vector.map (fun x -> exp (x * Complex(float dt, 0.0)))
            |> DenseMatrix.ofDiag

        (phi * A * b).Real()

    /// Projects a given Dmd model a number of times.
    /// Returns the last predicted vector using the intermediates as
    /// recursive initial conditions for the subsequent call to predict next.
    ///
    /// seriesDeltaTime: A number of time spanning each snapshot in seconds.
    /// numPredictions: Number of snapshot windows forward in time to project the model.
    /// u0: A initial condition vector to extrapolate the model forward in time from.
    /// model: A provided Dmd model solution to use in the projection.
    let extrapolateDmdModel (seriesDeltaTime: float<sec>) (numPredictions: int) (u0: Vector<float>) (model: DmdModel) =

        let folder u _ =
            predictNext seriesDeltaTime u model.Phi model.Omega

        [| 1..numPredictions |] |> Array.fold folder u0
