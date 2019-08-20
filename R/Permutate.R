#' @include Generics.R
#' @importFrom Matrix rowSums rowMeans
#' @importFrom stats na.omit sd pnorm
#'
NULL


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Methods for the FindPairs-defined generics

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# calculate the cell-cell interaciton
calculate_score <- function(object, ident, symA_no, symB_no, times, ...){
  if(is.null(ident) | !ident %in% c("stp", "mtp")){
    stop("Peasle input the ident of data")
  }
  else if(ident == "stp"){
    symA <- object@Assays$NData@FNData[, which(object@MData[, 1] == symA_no)]
    symB <- object@Assays$NData@FNData[, which(object@MData[, 1] == symB_no)]
  }
  else if(ident == "mtp"){
    symA <- object@Assays$NData@FNData[, which(object@MData[, 3] == paste0(symA_no, "/", times))]
    symB <- object@Assays$NData@FNData[, which(object@MData[, 3] == paste0(symB_no, "/", times))]
  }
  if(metric == "mean"){
    symA_mean <- rowMeans(as.matrix(symA))
    symB_mean <- rowMeans(as.matrix(symB))
    score <- exp(symA_mean[object@PPI$symbol_a] + symB_mean[object@PPI$symbol_b])
    names(score) <- object@PPI$symbol_ab
    score <- na.omit(score)
  }
  if(metric == "ratio"){
    symA_ratio <- rowSums(as.matrix(symA) > 0) / ncol(symA)
    symB_ratio <- rowSums(as.matrix(symB) > 0) / ncol(symB)
    score <- exp(symA_ratio[object@PPI$symbol_a] + symB_ratio[object@PPI$symbol_b])
    names(score) <- object@PPI$symbol_ab
    score <- na.omit(score)
  }
  return(score)
}


#' @rdname Calculate_score
#' @export
#'
Calculate_score.default <- function(object, ...) {
  print("You screwed up. I do not know how to handle the this object")
}


#' @param metric A metric of interaction. Set to two methods: ratio and mean.
#' \itemize{
#'   \item ratio: The percentage of cells with gene expression > 0
#'   \item mean: The mean of gene expression
#' }
#' @param target If this is the "stp" data, the "target" agrument is not considered. Calling this function will
#' directly calculate the intercellular interaciton. If this is the "mtp" data, the "target" has two potential
#' parameters. 1. "CCpTP", \strong{cell-cell per time point}: this calculate the interaciton of different cell types at
#' time poit. 2. "TPpCC", \strong{time point per cell-cell}: this calculate the cell-cell interaction of different time points
#' for each cell-cell interaciton.
#'
#' @rdname Calculate_score
#' @export
#' @method Calculate_score FindPairs
#'
Calculate_score.FindPairs <- function(object, metric = c("mean", "ratio"), target = NULL, ...) {
  symA_sum <- object@IA_Type[object@IA_Type[, 2] == "SYMBOLA", ][, 1]
  symB_sum <- object@IA_Type[object@IA_Type[, 2] == "SYMBOLB", ][, 1]
  Statistics <- list()
  object@P[["Target"]] <- c()
  # single time point intercellular interaction
  if(object@P[["Ident"]] == "stp") {
    data_stp <- data.frame()
    for(symB_no in symA_sum) {
      for(symA_no in symA_sum) {
        symA <- object@Assays$NData@FNData[, which(object@MData[, 1] == symA_no)]
        symB <- object@Assays$NData@FNData[, which(object@MData[, 1] == symB_no)]
        if(length(metric) > 1){
          stop("Please determine a measurement method")
        }
        if(metric == "mean"){
          symA_mean <- rowMeans(as.matrix(symA))
          symB_mean <- rowMeans(as.matrix(symB))
          score <- exp(symA_mean[object@PPI$symbol_a] + symB_mean[object@PPI$symbol_b])
          names(score) <- object@PPI$symbol_ab
          score <- na.omit(score)
        }
        if(metric == "ratio"){
          symA_ratio <- rowSums(as.matrix(symA) > 0) / ncol(symA)
          symB_ratio <- rowSums(as.matrix(symB) > 0) / ncol(symB)
          score <- exp(symA_ratio[object@PPI$symbol_a] + symB_ratio[object@PPI$symbol_b])
          names(score) <- object@PPI$symbol_ab
          score <- na.omit(score)
        }
        col_stp <- paste0(symA_no, "/", symB_no)
        data_stp <- add.col.self(dataframe = data_stp, new.vector = score)
        colnames(data_stp)[ncol(data_stp)] <- col_stp
      }
    }
    data_stp <- na.omit(data_stp)
    output <- new(Class = "Output",
                  output.strength = as.matrix(data.frame(data_stp)))
    Statistics["stp"] <- list(output)
    object@P[["Target"]] <- append(bject@P[["Target"]], "single time point")
    object@STP <- Statistics
  }
  # multi-times point intercellular interaction
  else if(object@P[["Ident"]] == "mtp") {
    object@MData$index <- paste0(object@MData[, 1], "/", object@MData[, 2])
    time_point <- unique(object@MData[, 2])
    # add input parament
    object@P[["Time_points"]] <- time_point
    if(is.null(target)) {
      stop("Please enter the target you want to calculate")
    }
    # time point per cell-cell interaction
    else if(target == "TPpCC") {
      for(symB_no in symB_sum) {
        for(symA_no in symA_sum) {
          data_tppcc <- data.frame()
          for(times in time_point) {
            symA <- object@Assays$NData@FNData[, which(object@MData[, 3] == paste0(symA_no, "/", times))]
            symB <- object@Assays$NData@FNData[, which(object@MData[, 3] == paste0(symB_no, "/", times))]
            if(length(metric) > 1){
              stop("Please determine a measurement method")
            }
            if(metric == "mean"){
              symA_mean <- rowMeans(as.matrix(symA))
              symB_mean <- rowMeans(as.matrix(symB))
              score <- exp(symA_mean[object@PPI$symbol_a] + symB_mean[object@PPI$symbol_b])
              names(score) <- object@PPI$symbol_ab
              score <- na.omit(score)
            }
            if(metric == "ratio"){
              symA_ratio <- rowSums(as.matrix(symA) > 0) / ncol(symA)
              symB_ratio <- rowSums(as.matrix(symB) > 0) / ncol(symB)
              score <- exp(symA_ratio[object@PPI$symbol_a] + symB_ratio[object@PPI$symbol_b])
              names(score) <- object@PPI$symbol_ab
              score <- na.omit(score)
            }
            data_tppcc <- add.col.self(dataframe = data_tppcc, new.vector = score)
          }
          colnames(data_tppcc) <- object@P[["Time_points"]]
          data_tppcc <- na.omit(data_tppcc)
          output <- new(Class = "Output",
                        output.strength = as.matrix(data.frame(data_tppcc)))
          nom <- paste0(symA_no, "/", symB_no)
          Statistics[nom] <- list(output)
        }
      }
      object@TPpCC <- Statistics
      object@P[["Target"]] <- append(object@P[["Target"]], "Time-Point per Cell-Cell interaction")
    }
    # cell-cell interaction per time point
    else if(target == "CCpTP") {
      for(tp in time_point){
        data_ccptp <- data.frame()
        for(symB_no in symB_sum){
          for(symA_no in symA_sum){
            symA <- object@Assays$NData@FNData[, which(object@MData[, 3] == paste0(symA_no, "/", tp))]
            symB <- object@Assays$NData@FNData[, which(object@MData[, 3] == paste0(symB_no, "/", tp))]
            if(length(metric) > 1){
              stop("Please determine a measurement method")
            }
            if(metric == "mean"){
              symA_mean <- rowMeans(as.matrix(symA))
              symB_mean <- rowMeans(as.matrix(symB))
              score <- exp(symA_mean[object@PPI$symbol_a] + symB_mean[object@PPI$symbol_b])
              names(score) <- object@PPI$symbol_ab
              score <- na.omit(score)
            }
            if(metric == "ratio"){
              symA_ratio <- rowSums(as.matrix(symA) > 0) / ncol(symA)
              symB_ratio <- rowSums(as.matrix(symB) > 0) / ncol(symB)
              score <- exp(symA_ratio[object@PPI$symbol_a] + symB_ratio[object@PPI$symbol_b])
              names(score) <- object@PPI$symbol_ab
              score <- na.omit(score)
            }
            data_ccptp <- add.col.self(dataframe = data_ccptp, new.vector = score)
            colnames(data_ccptp)[ncol(data_ccptp)] <- paste0(symA_no, "/", symB_no)
          }
        }
        data_ccptp <- na.omit(data_ccptp)
        output <- new(Class = "Output",
                      output.strength = as.matrix(data.frame(data_ccptp)))
        Statistics[tp] <- list(output)
      }
      object@CCpTP <- Statistics
      object@P[["Target"]] <- append(object@P[["Target"]], "Cell-Cell interaction per Time-Point")
    }
  }
  # add parameter
  object@P[["metric"]] <- metric
  return(object)
}



#' @rdname permutate_test
#' @export
#'
permutate_test.default <- function(object, ...) {
  print("You screwed up. I do not know how to handle the this object")
}



#' @param number Permutation numbers
#' @param i The names of intercellular interaction, such as "cellA/CellB"
#' @param Asample_size The number of cell expressing symbol A virtually. We generalized the concept of permutation test.
#' In the traditional permutation test, the number of random substitutions should be the total number of
#' corresponding targets. But, in our permutation test, we set this number as an adjustable parameter. This paramter
#' have two funcitons: one is to provide a fixed background distribution, and other is to reduce the impact of low number
#' cells in the permutation.
#' @param Bsample_size The number of cell expressing symbol B virtually.
#' @param seed random seed
#'
#'
#' @rdname permutate_test
#' @export
#' @method permutate_test FindPairs
#'
permutate_test.FindPairs <- function(object, number, i, Asample_size, Bsample_size, seed, ...) {
  symAcell <- unlist(strsplit(i, "/"))[1]
  symBcell <- unlist(strsplit(i, "/"))[2]
  exp_dis <- data.frame()
  for(a in 1 : number){
    set.seed(seed + a * 3)
    symA_s <- object@Assays$NData@FNData[, sample(which(object@MData[, 1] == symAcell), Asample_size)]
    symB_s <- object@Assays$NData@FNData[, sample(which(object@MData[, 1] == symBcell), Bsample_size)]
    if(object@P$metric == "ratio") {
      symA_sratio <- rowSums(as.matrix(symA_s) > 0) / ncol(symA_s)
      symB_sratio <- rowSums(as.matrix(symB_s) > 0) / ncol(symB_s)
      score <- exp(symA_sratio[object@PPI$symbol_a] + symB_sratio[object@PPI$symbol_b])
    }
    if(object@P$metric == "mean") {
      symA_smean <- rowMeans(as.matrix(symA_s))
      symB_smean <- rowMeans(as.matrix(symB_s))
      score <- exp(symA_smean[object@PPI$symbol_a] + symB_smean[object@PPI$symbol_b])
    }
    exp_dis <- add.col.self(dataframe = exp_dis, new.vector = score)
  }
  rownames(exp_dis) <- object@PPI$symbol_ab
  exp_dis <- na.omit(exp_dis)
  # pvalue
  pvl_one <- apply(object@Output[[i]]@output.strength, 2, function(x){
    u <- (x - mean(exp_dis))/sd(exp_dis)
    pvl <- 1 - pnorm(q = u, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  })
  colnames(pvl_one) <- object@P[["Ident"]]
  return(pvl_one)
}



#' @rdname FP_test
#' @export
#'
FP_test.default <- function(object, ...){
  print("You screwed up. I do not know how to handle the this object")
}



#' @param Asample_size Permutating the number of cell expressing symbol A. Default 200.We generalized the concept of
#' permutation test. In the traditional permutation test, the number of random substitutions should be the total number
#' of corresponding targets(the target cell type). But, in our permutation test, we set this number as an adjustable
#' parameter. This paramterhave two funcitons: one is to provide a fixed background distribution, and other is to
#' reduce the impact of low number cells in the permutation.
#' @param Bsample_size Permutating the number of cell expressing symbol B. Default 200'
#' @param fdr Filter for significant pvalue. Default 0.01
#' @param permutation_number The number of permutation. Default 5000
#' @param fdr_method Corrected method in multiple hypothesis testing. Default "BH".
#' @param dopar Whether Parallel computation is allowed. Default FLASE.
#' @param dopar.numbers If the \code{do.par = TRUE}, multiple cores will be used for parallel operations. Default 2.
#' @param random.seed Set the random seed. Default 42
#'
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#'
#' @return
#'
#' @rdname FP_test
#' @export
#' @method FP_test FindPairs
#'
FP_test.FindPairs <- function(
  object,
  Asample_size = 3000,
  Bsample_size = 3000,
  fdr = 0.01,
  fdr_method = "BH",
  permutation_number = 5000,
  dopar = FALSE,
  dopar.numbers = 2,
  random.seed = 42) {
  if(length(object@Output) != 0){
    if(dopar == TRUE) {
      cl <- makeCluster(dopar.numbers)
      registerDoParallel(cl)
      pvl_list <- foreach(name = names(object@Output), .export = c("permutate_test", "permutate_test.default", "permutate_test.FindPairs", "add.col.self"),
                          .packages = c("Matrix")) %dopar% {
                            pvl_one <- permutate_test(object, number = permutation_number, i = name,
                                                      seed = random.seed, Asample_size = Asample_size, Bsample_size = Bsample_size)

                          }
      stopCluster(cl)
      names(pvl_list) <- names(object@Output)
    }
    else {
      pvl_list <- list()
      for(name in names(object@Output)) {
        pvl_one <- permutate_test(object, number = permutation_number, i = name,
                                  seed = random.seed, Asample_size = Asample_size, Bsample_size = Bsample_size)
        pvl_list[name] <- list(pvl_one)
      }
    }
    for(p in names(pvl_list)) {
      object@Output[[p]]@output.pvalue <- pvl_list[[p]]
    }
    return(object)
  }
  else{
    stop("Please calculate the interaciton score by calling 'Calculate_score'")
  }
}



#' Add column without restriction
#' Add new column to the new dataframe without restriction
#'
#' @param dataframe The data.frame with unrestricted shape, including the empty dataframe
#' @param new.vector The vector whith unrestricted length
#'
#' @return A new matrix
#'
add.col.self <- function(dataframe, new.vector) {
  nl <- list(dataframe, new.vector)
  nl <- lapply(nl, as.matrix)
  no <- max(sapply(nl, nrow))
  cbinddata <- do.call(cbind, lapply(nl, function(x) rbind(x, matrix(, no - nrow(x), ncol(x)))))
  return(cbinddata)
}

#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

# plot the heatmap


