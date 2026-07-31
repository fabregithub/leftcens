#' Groundwater trace-metal concentrations (interval layout)
#'
#' A small simulated example of the kind of table an analytical laboratory
#' exports, in the **three-tier `.left` / `.right`** layout: one row per sample
#' and a lower/upper interval-bound pair for each metal. Concentrations are in
#' micrograms per litre. Each metal shows a mix of non-detect `(0, MDL)`,
#' detected-not-quantified `(MDL, LCMRL)`, and quantified `(value, value)`
#' results. Feed it straight to [as_cens_list()].
#'
#' @format A data frame with 28 rows and 11 columns:
#' \describe{
#'   \item{sample_id}{sample identifier}
#'   \item{As.left, As.right}{arsenic lower/upper bound (ug/L)}
#'   \item{Cd.left, Cd.right}{cadmium lower/upper bound (ug/L)}
#'   \item{Pb.left, Pb.right}{lead lower/upper bound (ug/L)}
#'   \item{Cu.left, Cu.right}{copper lower/upper bound (ug/L)}
#'   \item{Zn.left, Zn.right}{zinc lower/upper bound (ug/L)}
#' }
#' @source Simulated for illustration; see `data-raw/sample_data.R`.
#' @seealso [surfacewater] for the same measurements in a qualifier-flag layout.
#' @examples
#' cl <- as_cens_list(groundwater)
#' as.data.frame(desc_np(cl))
"groundwater"

#' Surface-water trace-metal concentrations (flagged layout)
#'
#' A small simulated example in the **two-tier value + qualifier-flag** layout,
#' the most common laboratory export: one row per sample, and for each metal a
#' value column plus a `.flag` column where `"<"` marks a non-detect (the value
#' being the detection/reporting limit) and a blank marks a quantified result.
#' Concentrations are in micrograms per litre. Read it with
#' `as_cens_list(surfacewater, format = "flagged")`.
#'
#' @format A data frame with 30 rows and 11 columns:
#' \describe{
#'   \item{sample_id}{sample identifier}
#'   \item{As, As.flag}{arsenic value (ug/L) and qualifier flag}
#'   \item{Cd, Cd.flag}{cadmium value (ug/L) and qualifier flag}
#'   \item{Pb, Pb.flag}{lead value (ug/L) and qualifier flag}
#'   \item{Cu, Cu.flag}{copper value (ug/L) and qualifier flag}
#'   \item{Zn, Zn.flag}{zinc value (ug/L) and qualifier flag}
#' }
#' @source Simulated for illustration; see `data-raw/sample_data.R`.
#' @seealso [groundwater] for the interval (`.left`/`.right`) layout.
#' @examples
#' cl <- as_cens_list(surfacewater, format = "flagged")
#' as.data.frame(desc_np(cl))
"surfacewater"
