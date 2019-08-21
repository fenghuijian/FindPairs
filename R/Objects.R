#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Create Objects

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom Matrix colSums rowSums colMeans rowMeans
#' @importFrom methods setClass setOldClass setClassUnion slot
#' slot<- setMethod new signature slotNames is
#' @importClassesFrom Matrix dgCMatrix
#'
NULL


setClassUnion(name = "AnyMatrix",
              members = c("matrix", "dgCMatrix")
)


#' The FindPairs Class
#'
#' The FindPairs Class is a repesentation of single-cell expression data for R.
#' This Clsss is set to better calculate the intercellular intercation.
#'
#' @slot Assays A list of assays of this project
#' @slot MData Contains meta-information about each cells
#' @slot IA_Type User-defined intercellular interaction need to be calculated
#' @slot P User-defined parameters
#' @slot PPI Curated Protein-Protein Interaciton Database
#' @slot DEGs A list of output of this differentially expressed genes
#' @slot STP A list of output of the "single time point" data
#' @slot TPpCC A list of output of the target: "time-point per cell-cell interaction" of "multi-time points" data
#' @slot CCpTP A list of output of the target: "cell-cell interaction per time-point" of "multi-time points" data
#' @slot Enrich A list of output of the enrcihment
#'
#' @name FindPairs-class
#' @rdname FindPairs-class
#' @exportClass FindPairs
#'
FindPairs <- setClass(Class = "FindPairs",
                      slots = c(
                        # user input data
                        Assays = "list",
                        MData = "data.frame",
                        IA_Type = "data.frame",
                        # user input parameter
                        P = "list",
                        PPI = "data.frame",
                        # deg for the ident
                        DEGs = "list",
                        #output result
                        # single time point
                        STP = "list",
                        # select time-point per cell-cell
                        TPpCC = "list",
                        # select cell-cell per time-point
                        CCpTP = "list",
                        # enrich result
                        Enrich = "list"
                      ))


#' The Assays Class
#'
#' The Assays Class is an intermediate data storage class that stores the raw normalized data
#' and the normalized data filtered by Intercellular Interaction Database.
#'
#' @slot RNData An original normalized data
#' @slot FNData An filtered normalized data
#'
#' @name Assays-class
#' @rdname Assays-class
#' @exportClass Assays
#'
Assays <- setClass(Class = "Assays",
                   slots = c(
                     RNData = "AnyMatrix",
                     FNData = "AnyMatrix"
                   ))


#' The Output object
#' The Output class store the caculated result
#'
#' @slot output.pvalue output the pvalue result
#' @slot output.strength output the interactive strength result
#'
#' @name Output-class
#' @rdname Output-class
#' @exportClass Output
#'
Output <- setClass(Class = "Output",
                   slots = c(
                     output.pvalue = "AnyMatrix",
                     output.strength = "AnyMatrix"
                   ))


#' Setup the FindPairs object
#'
#' @importFrom utils read.csv
#' @importFrom methods new
#'
#' @param data the normalized data of gene expression
#' @param meta.data annotation of cells. "Cell types" must be annotated. If it is multi-times-point, the times
#' information should be added to the dataframe.
#' @param interactive.type The \strong{feature} expressed by \strong{cell types}. The dataframe must contain two variables.
#' This parameter needs to be set manually, and it must follow the  format below.
#' \itemize{
#'   \item the \strong{cell types} must be in the first column
#'   \item the \strong{gene feature} must be in the second column, and \strong{gene feature} have two varibales:
#'   \strong{"SYMBOLA"} means this gene feature is "SYMBOLA" feature and \strong{"SYMBOLB"} means this gene feature is
#'   "SYMBOLB" feature.
#' }
#' @param ident Determine whether the input data contains multi-times information. There are two arguments, mtp:
#' multi-times-point or stp: single-time-point.
#' @export
#'
BuildFP <- function(data,
                    meta.data,
                    interactive.type,
                    ident = c("mtp", "stp")){
  if(length(ident) > 1) {
    stop("Please determine a attribute")
  }
  else if(ident %in% c("mtp", "stp") == FALSE) {
    stop("Please enter a valid attribute")
  }
  else{
    # filter data
    filter.data <- data[rownames(data) %in% c(ppidb$symbol_a, ppidb$symbol_b), ]
    assays <- new(Class = "Assays",
                  RNData = data,
                  FNData = filter.data)
    # parmeter
    ob <- new(Class = "FindPairs",
              Assays = list(NData = assays),
              MData = meta.data,
              P = list(Ident = ident),
              IA_Type = interactive.type,
              PPI = ppidb)
    return(ob)
  }
}





