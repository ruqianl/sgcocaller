## diagnostic txt to png

## phase/\
.libPaths("/mnt/beegfs/mccarthy/scratch/general/rlyu/Software/R/4.0/library/")

args <- (commandArgs(trailingOnly = TRUE))

for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

print(pngfile)
print(inputfile)
print(snpAnnotFile)


tryCatch(
 {
   library(ggplot2) 
   library(dplyr)
   library(tidyr)
   message("Loading pacakges")
  
      # for the condition handlers for warnings and error below)
  },
  error=function(cond) {
  
    message("required R pacakges for plotting diagnostic plots are not available")

    return(NA)
  })

plot_df <- read.table(file = inputfile,
                      header = T)
snpAnnot_df <- read.table(file = snpAnnotFile, header =TRUE)
snpAnnot_df <- snpAnnot_df[snpAnnot_df$Phase!="-1",]
plot_df_match <- apply(plot_df,2,function(x){
  tmp <-  (x == plot_df$templateGeno)
  tmp[x == 0] <-  NA
  tmp
})


data.frame(cbind(plot_df_match,snpAnnot_df)) %>%
  tidyr::pivot_longer(cols = colnames(plot_df)) %>%
  dplyr::filter(!is.na(value),POS > 1.3e8, POS < 1.35e8) %>%
  ggplot()+
  geom_point(mapping = aes(x = POS , y = name, color = value))+
  scale_color_manual(values = c("FALSE"='red',
                                "TRUE" = "blue"))+
    theme_bw(base_size = 18)+geom_vline(mapping = aes(xintercept = c(switch_pos[2])),
                                        color = "red",size = 1.2)

ggsave(pngfile, dpi = 100, width = 14, height = 6)
