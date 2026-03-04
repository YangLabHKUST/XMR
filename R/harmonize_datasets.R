#' @export
comple <- function(allele){

  ifelse(allele == "A","T", ifelse(allele == "T","A", ifelse(allele == "G","C", ifelse(allele == "C","G", allele)) ))

}

#' @export
# facilitate harmonization of two GWAS datasets
align_data <- function(dat1,
                      dat2,
                      trait1.name = "exposure",
                      trait2.name = "outcome",
                      LDSC = T,
                      h2.fix.intercept = F,
                      ldscore.dir = NULL,
                      ld=NULL,
                      M=NULL){
  
    dat1 <- dat1 %>% dplyr::mutate_if(is.integer, as.numeric)
    dat2 <- dat2 %>% dplyr::mutate_if(is.integer, as.numeric)
    dat1 <- dat1 %>% dplyr::mutate_if(is.factor, as.character)
    dat2 <- dat2 %>% dplyr::mutate_if(is.factor, as.character)

    message("Merge dat1 and dat2 by SNP ...")
    dat = merge(dat1, dat2, by="SNP")

    message("Harmonize the direction of SNP effects of exposure and outcome")
    flip.index = which((dat$A1.x == dat$A2.y & dat$A1.y == dat$A2.x) |
                       (dat$A1.x ==comple(dat$A2.y) & dat$A1.y == comple(dat$A2.x)))

    dat[,"A1.y"] = dat[,"A1.x"]
    dat[,"A2.y"] = dat[,"A2.x"]
    dat[flip.index ,"Z.y"] = -dat[flip.index ,"Z.y"]

    #message("Add LD scores to the harmonized data set...")
    #merged  = merge(dat, ld, by="SNP")
    merged = dat

    message("The Harmonized data set will also be used for MR analysis \n")
    dat = data.frame(SNP = merged$SNP,
                   A1 = merged$A1.x,
                   A2 = merged$A2.x,
                   b.exp = merged$Z.x/sqrt(merged$N.x),
                   b.out = merged$Z.y/sqrt(merged$N.y),
                   se.exp = 1/sqrt(merged$N.x),
                   se.out = 1/sqrt(merged$N.y),
                   pval.exp = merged$P.x,
                   pval.out = merged$P.y
                  )

    return(list(dat=dat))
}
