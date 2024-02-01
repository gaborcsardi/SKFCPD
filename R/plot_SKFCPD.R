##########################################################################
## plot_SKFCPD function
## 
## SKFCPD Package
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, April 2013.
##
## Copyright (C) 2023-present Hanmo Li, Yuedong Wang, Mengyang Gu
##  						  
##    
##########################################################################

plot_SKFCPD <- function(x) {
  
  model = x
  n_obs = nrow(model@design)
  output_dim_plot = 1
  
  x_test = as.matrix(as.matrix(model@design)[model@test_start:n_obs,])
  y_test = as.matrix(as.matrix(model@response)[model@test_start:n_obs,output_dim_plot])
  
  temp_df = data.frame(input = x_test,
                       output = y_test)

  plot_1 = ggplot(temp_df, aes(x = .data$input, y = .data$output))+
    geom_line() +
    theme_bw() +
    geom_vline(xintercept = model@design[model@cp], linetype="dotted", 
               color = "red", linewidth=1)
  
  run_length_posterior_mat = model@run_length_posterior_mat
  
  runlength_index = apply(run_length_posterior_mat, 2, which.max)
  runlength_df = data.frame(date = temp_df$input,
                            runlength = runlength_index)
  
  run_length_posterior_mat[run_length_posterior_mat<10^(-5)] = 0
  
  d2.df <- melt(t(run_length_posterior_mat), c("x", "Run length"), value.name = "p(Run length)")
  d2.df$input = rep(temp_df$input, length(temp_df$input))
  d2.df$runlength = rep(runlength_index, length(temp_df$input))
  d2.df$`p(Run length)` = d2.df$`p(Run length)` + 10^(-5)
  plot_2 = ggplot(data=d2.df,aes(x=.data$input,y=.data$`Run length`,fill=.data$`p(Run length)`))+
    geom_tile()+
    scale_fill_gradient(
      low = "#FFFFFF",
      high = "black",
      na.value = "#FFFFFF",
      aesthetics = "fill",
      trans = "log10",
      guide = guide_colourbar(direction = "horizontal", title.position = "top",
                                   label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                   label.theme = element_text(angle = 90))
    ) + 
    xlab("Input")+
    ylab("Run length")+
    scale_y_continuous(trans = "reverse")+
    theme_bw()+
    theme(legend.position = c(0.2, 0.3))+
    geom_line(aes(x = .data$input, y = .data$runlength), color = "red", linetype = "dashed", linewidth = 1.5)

  
  ggarrange(plot_1, plot_2, nrow = 2, ncol = 1,  heights = c(1, 1), labels = c("a", "b"))
  
}
