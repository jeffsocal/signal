# Jeff Jones
# 2019.02.07
#
# progress timer wrapper

library(progress)

progress_bar <- function(size, text = "running ..."){
  cat(text, "\n")
  
  pb <- progress_bar$new(
    format = paste(size, ":percent :elapsed [:bar] :eta"),
    total = size, clear = FALSE, width= 60)
  
  return(pb)
}
