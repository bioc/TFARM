#' Generates an heatmap visualization of the mean Importance Index
#' of pairs of transcription factors.
#'
#' A square matrix (TFs x TFs) is built. The element (i,j) of such matrix
#' contains the mean Importance Index (II) of the couple of transcription
#' factors (TFi, TFj). This is rapresented in a heatmap visualization in
#' which the scale color from red (dark) to white (light) indicates low mean
#' importance of the couple (red), or high mean importance of the
#' couple (white).
#'
#' @param TFs a string vector with names of transcription factors.
#' @param I a vector of mean II of pairs of transcription factors in TFs.
#'
#' @return The Importance index of pairs of Transcription Factors as a heatmap
#' @export
#' @import fields
#' @importFrom grDevices heat.colors
#' @importFrom graphics par axis title
#' @examples
#' # Load p_TFs and I_c_2 from the data_man collection of datasets:
#' data('data_man')
#'
#' # Heatmap visualization of the mean importances of transcription factors in p
#' # and their combinations in two elements:
#' heatI(p_TFs, I_c_2)

heatI <- function(TFs, I) {
    estract <- function(x) {
        vec <- vapply(TFs, function(TFs) gregexpr(TFs, x)[[1]][1] != -1,
                      numeric(1))
        return(vec)
    }
    l <- length(TFs)
    aa <- t(sapply(I[, 1], estract))
    matrix_imp <- matrix(0, l, l)
    pos <- apply(aa == 1, 1, function(x) which(x)[1:2])
    matrix_imp[cbind(pos[1, which(is.na(pos[2, ]))], pos[1, which(is.na(pos[2,
                                                                            ]))])] <- I[which(is.na(pos[2, ])), "imp"]
    matrix_imp[cbind(pos[1, which(!is.na(pos[2, ]))], pos[2, which(!is.na(pos[2,
                                                                              ]))])] <- I[which(!is.na(pos[2, ])), "imp"]
    matrix_imp[cbind(pos[2, which(!is.na(pos[2, ]))], pos[1, which(!is.na(pos[2,
                                                                              ]))])] <- I[which(!is.na(pos[2, ])), "imp"]

    rownames(matrix_imp) = TFs
    colnames(matrix_imp) = TFs
    par(mar = c(7, 7, 4, 10))
    image(matrix_imp, axes = FALSE, xlab = "", ylab = "", col = heat.colors(20),
          ylim = c(1.04, -0.04))
    axis(1, at = seq(0, 1, length = l), labels = rownames(matrix_imp),
         las = 2, cex.axis = 0.8)
    axis(2, at = seq(0, 1, length = l), labels = rownames(matrix_imp),
         las = 2, cex.axis = 0.8)
    im <- image.plot(matrix_imp, col = heat.colors(20), legend.only = TRUE,
                     ylim = c(1.04, -0.04))
    title(main = "Importance index of pairs of Transcription Factors")
    return(im)
}
