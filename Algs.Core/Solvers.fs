module Solvers

open MathNet.Numerics
open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.LinearAlgebra.Double

open Matrix.Extensions

/// The Sequential Threshold Least Squares algorithm.
///
/// numIters: Number of times to regress terms.
/// lambda: Sparsity promoting threshold used to truncate small terms.
/// A: System matrix
/// B: Regressor matrix
///
/// Algorithm adapted from:
/// https://arxiv.org/pdf/1509.03580.pdf.
let stlsq numIters lambda (A: DenseMatrix) (B: DenseMatrix) =
    let allRows = Seq.replicate A.RowCount 1 |> vector
    let mutable (Xi: DenseMatrix) = DenseMatrix.OfMatrix(A.Solve(B))

    do
        for _ in 1..numIters do

            // Truncate small values toward zero
            let mutable smallIndices = Xi.PointwiseAbsLessThan lambda
            Xi[smallIndices] <- 0.0

            // Regress dynamics onto remaining terms to find sparse Xi
            let activeIndices = DenseMatrix.PointwiseNot smallIndices

            for (colIndex, bigInds) in activeIndices.EnumerateColumns() |> Seq.indexed do
                let Sq = A[allRows, bigInds]
                let neSol = Sq.Solve(B[*, colIndex])
                Xi[bigInds, colIndex] <- neSol

    Xi

/// An implmentation of the LASSO algorithm using
/// Alternating Direction Method of Multipliers.
///
/// A: the system matrix
/// y: the regressor vector to solve.
/// lambda: resgularisation parameter
/// rho: augmented lagrangian parameter
/// alpha: relaxation parameter
/// absTol: Aboslute tolerance
/// relTol: Relative tolerance
/// matIters: The maximum number of iterations
///
/// Adapted from:
/// https://arxiv.org/ftp/arxiv/papers/2208/2208.11544.pdf
let lasso (A: Matrix<float>) (y: Vector<float>) (lambda: float) (rho: float) (alpha: float) absTol relTol maxIters =
    let n = A.RowCount
    let m = A.ColumnCount
    let I = CreateMatrix.DenseIdentity m

    let mutable x: Vector<float> = DenseVector.Create(m, 0.0)
    let mutable z: Vector<float> = DenseVector.Create(m, 0.0)
    let mutable u: Vector<float> = DenseVector.Create(m, 0.0)

    let mutable primResidual = 1e1
    let mutable dualResidual = 1e1

    let shrink eps (s: Vector<float>) =
        s.Enumerate()
        |> Seq.map (fun x -> if abs (x) < eps then 0.0 else x)
        |> DenseVector.OfEnumerable

    let mutable primTol = 1e0
    let mutable dualTol = 1e0

    let rec loop k =
        if k < maxIters && (primResidual > primTol || dualResidual > dualTol) then

            // Update x
            let x1 = A.TransposeThisAndMultiply(A) + rho * I
            let x' = x1.Solve(A.TransposeThisAndMultiply(y) + rho * (z - u))
            x <- alpha * x' + (1.0 - alpha) * z

            // Update z
            let z' = z
            z <- shrink (lambda / rho) (x + u)

            // Update u
            u <- u + x - z

            // Compute residuals
            primResidual <- (x - z).L2Norm()
            dualResidual <- (-rho * (z - z')).L2Norm()

            primTol <- sqrt (float m) * absTol + relTol * max (x.L2Norm()) (z.L2Norm())
            dualTol <- sqrt (float m) * absTol + relTol * (rho * u).L2Norm()

            loop (k + 1)
        else
            x, primResidual, dualResidual, k

    loop 0
