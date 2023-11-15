module Matrix.Extensions

open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.LinearAlgebra.Double

/// The following code extendes the MathNet Matrix type to
/// allow logical row and column indexing.
///
/// No work has been done to benchmark performance.
type DenseMatrix with

    /// Support for logical indexing using ints
    member this.Item
        with set (logical: Matrix<int>) value =

            if logical.RowCount <> this.RowCount || logical.ColumnCount <> this.ColumnCount then
                failwith "Logical Indexes must have same shape"

            do
                for (rowIndex, row) in logical.EnumerateRows() |> Seq.indexed do
                    for (colIndex, elem) in row |> Seq.indexed do
                        if elem <> 0 then
                            this.Item(rowIndex, colIndex) <- value

    member this.Item
        with get ((logicalRows: Vector<int>), (logicalCols: Vector<int>)) =
            if logicalRows.Count <> this.RowCount then
                failwith "Error logical rows must have same count"

            if logicalCols.Count <> this.ColumnCount then
                failwith "Error logical columns must have same count"

            let rowIndices = { 0 .. this.RowCount - 1 }
            let colIndices = { 0 .. this.ColumnCount - 1 }

            // These next operations may not be efficient
            let slicedRows =
                Seq.zip rowIndices logicalRows
                |> Seq.filter (fun (_, selected) -> selected = 1)
                |> Seq.map (fun (index, _) -> this[index, *])
                |> DenseMatrix.ofRowSeq

            Seq.zip colIndices logicalCols
            |> Seq.filter (fun (_, selected) -> selected = 1)
            |> Seq.map (fun (index, _) -> slicedRows[*, index])
            |> DenseMatrix.ofColumnSeq

    member this.Item
        with set ((logicalRows: Vector<int>), colIndex) values =

            if logicalRows.Count <> this.RowCount then
                failwith "Error logical rows must have same count"

            if logicalRows.Count <> (values |> Seq.length) then
                failwith "Error row count must equal values count"

            let rowIndices = { 0 .. this.RowCount - 1 }

            let selectedRows =
                Seq.zip rowIndices logicalRows
                |> Seq.filter (fun (_, selected) -> selected = 1)
                |> Seq.map (fun (index, _) -> index)

            Seq.zip selectedRows values
            |> Seq.iter (fun (rowIndex, value) -> this[rowIndex, colIndex] <- value)

    member this.PointwiseAbsLessThan value =
        seq {
            for row in this.EnumerateRows() do
                seq { for elem in row -> if abs (elem) < value then 0 else 1 }
        }
        |> matrix

    static member PointwiseNot(m: Matrix<int>) =
        seq {
            for row in m.EnumerateRows() do
                seq { for elem in row -> if elem = 1 then 0 else 1 }
        }
        |> matrix
