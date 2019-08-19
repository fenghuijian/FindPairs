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
#' @slot PPi Curated Protein-Protein Interaciton Database
#' @slot Output A list of Output of this project
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
                        # load the data
                        PPI = "data.frame",
                        # outresult
                        Output = "list"
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
#' @param meta.data annotation of
#' @param interactive.type The \strong{feature} expressed by \strong{cell types}. The dataframe must contain two variables.
#' This parameter needs to be set manually, and it must follow the  format below.
#' \itemize{
#'   \item the \strong{cell types} must be in the first column
#'   \item the \strong{gene feature} must be in the second column, and \strong{gene feature} have two varibales:
#'   \strong{"SYMBOLA"} means this gene feature is "SYMBOLA" feature and \strong{"SYMBOLB"} means this gene feature is
#'   "SYMBOLB" feature.
#' }
#' @param ident input the cells day
#' @export
#'
BuildFP <- function(data,
                    meta.data,
                    interactive.type,
                    ident){
  options(stringsAsFactors = F)
  # filter data
  filter.data <- data[rownames(data) %in% c(ppidb$symbol_a, ppidb$symbol_b), ]
  assays <- new(Class = "Assays",
                RNData = data,
                FNData = filter.data)
  # metadata
  meta.data$index <- paste0(meta.data[, 1], "/", meta.data[, 2])
  # parmeter
  ob <- new(Class = "FindPairs",
            Assays = list(NData = assays),
            MData = meta.data,
            P = list(Ident = ident),
            IA_Type = interactive.type,
            PPI = ppidb
  )
}


