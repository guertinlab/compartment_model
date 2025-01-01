transform_rates <- function(combinedResultsFile, baselineName, exprName, yLabel, fileidentifier=NULL, kelong_vals = c(30,45,60), kpre_vals = c(1e-3,1e-2,1e-1,1), ALPHA=0.03) {
  combinedResults = read.table(combinedResultsFile,sep="\t", header=T)
  if (is.null(fileidentifier)) {
     fileidentifier = sprintf("%s_vs_%s", exprName, baselineName)
     print(sprintf("setting filePrefix to %s", fileidentifier))
 }
   p_b_df = combinedResults[,c("gene", sprintf("P_normalized_%s", baselineName), sprintf("B_normalized_%s", baselineName),
                                           sprintf("P_normalized_%s", exprName), sprintf("B_normalized_%s", exprName))]
   colnames(p_b_df) = c("gene", "P_base", "B_base", "P_expr", "B_expr")
  rate.df = as.data.frame(c())
   for (kelong_val in kelong_vals) {
       p_b_df$kelong = kelong_val
         estimatedRate = p_b_df
         for (kpre_val in kpre_vals) {
           estimatedRate$kpre = kpre_val
           estimatedRate$kinit_base = estimatedRate$kpre * estimatedRate$P_base + estimatedRate$kelong * estimatedRate$B_base
           estimatedRate$kinit_expr = estimatedRate$kpre * estimatedRate$P_expr + estimatedRate$kelong * estimatedRate$B_expr
           estimatedRate$krel_base = estimatedRate$kelong * ( estimatedRate$B_base / estimatedRate$P_base )
           estimatedRate$krel_expr = estimatedRate$kelong * ( estimatedRate$B_expr / estimatedRate$P_expr )
           rate.df = rbind(rate.df, estimatedRate)
         }
  }
 saveRDS(rate.df, sprintf("calculatedRateData_%s.rds", fileidentifier))

}


