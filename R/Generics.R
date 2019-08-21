#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# S4 methods

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


setMethod(f = "show",
          signature = "FindPairs",
          definition = function(object){
            cat("An object of class", class(x = object), "\n")
            cat("\n",
                "Original Normalized Data:\n",
                "nGenes x nCells\n",
                paste0(nrow(object@Assays$NData@RNData), " x ", ncol(object@Assays$NData@RNData)), "\n",
                sep = "")
            cat("\n",
                "Decting the intercellular interaction:\n",
                sep = "")
            print(object@IA_Type)
          })


setMethod(f = "[[",
          signature = c("x" = "FindPairs"),
          definition = function(x, i, ..., value){
            if(!is.character(x = i)){
              stop("'i' must be a character", call. = FALSE)
            }
            if(length(x = i) > 1){
              stop("Please input a unique character", call. = FALSE)
            }
            if(i == "FNData") {
              return(x@Assays$NDdata@FNData)
            }
            if(i == "MData") {
              return(x@MData)
            }
            if(i == "IA_Type") {
              return(x@IA_Type)
            }
            if(i == "PPI") {
              return(x@PPI)
            }
            if(i %in% names(x@Output)){
              return(x@Output[[i]])
            }
          })


setMethod(f = "[[<-",
          signature = c("x" = "FindPairs"),
          definition = function(x, i, ..., value){
            if(!is.character(x = i)) {
              stop("'i' must be a character", call. = FALSE)
            }
            if(length(x = i) > 1) {
              stop("Please input a unique character", call. = FALSE)
            }
            if(i == "MData") {
              x@MData <- value
            }
            else if(i == "IA_Type") {
              x@IA_Type <- value
            }
            else{
              stop("There is no slot called ", i, call. = FALSE)
            }
            return(x)
          })


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# S3 methods

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Calculate the interacellular interacion score
#'
#' Calculate the cell-cell interaction score.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns FindPairs object with the score calculation stored in the Output slot
#'
#' @export
#'
#' @rdname Calculate_score
#' @export Calculate_score
#'
Calculate_score <- function(object, ...) {
  print("Calculate the interaction cellular score")
  UseMethod(generic = "Calculate_score", object = object)
}


#' Time-Point per Cell-Cell interaction Permutate Test
#'
#' Permutate Test for different time point of one Cell-Cell interaction
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname tppcc_ptest
#' @export tppcc_ptest
#'
tppcc_ptest <- function(object, ...) {
  UseMethod(generic = "tppcc_ptest", object = object)
}



#' Cell-Cell interaction per Time-Point Permutate Test
#'
#' Permutate Test for different cell-cell interaction of one time point
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname ccptp_ptest
#' @export ccptp_ptest
#'
ccptp_ptest <- function(object, ...) {
  UseMethod(generic = "ccptp_ptest", object = object)
}


#' Cell-Cell interaction Permutate Test
#'
#' Permutate Test for different cell-cell interaction without time-point
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname stp_ptest
#' @export stp_ptest
#'
stp_ptest <- function(object, ...) {
  UseMethod(generic = "stp_ptest", object = object)
}


#' FindPairs object permutate test
#' Permutate Test on all intercellular interaction
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns FindPairs object with the pvalue calculation stored in the Output slot
#'
#' @export
#'
#' @rdname FP_test
#' @export FP_test
#'
FP_test <- function(object, ...) {
  UseMethod(generic = "FP_test", object = object)
}



















