# general configurations for plots

#set ggplot output
ggplot2::theme_set(ggplot2::theme_minimal(base_size = 10))

ggplot2::theme_update(text = element_text(colour = "grey20", size = 10), 
                      legend.text = element_text(colour = "grey20", size = 10),
                      legend.title = element_text(colour = "grey20", size = 10),
                      axis.text = element_text(colour = "grey40", size = 10),
                      axis.ticks = element_line(colour = "grey50"),
                      strip.text = element_text(colour = "grey20", size = 10),
                      panel.grid.minor = element_blank(),  
                      panel.grid.major = element_blank(),
                      plot.title = element_text(colour = "grey20", size = 10, 
                                                face = "bold"), 
                      plot.subtitle = element_text(colour = "grey20", size = 10,
                                                   face =  "italic"))

# define output sizes
image_width <- 183
image_height <- 100
image_units <- "mm"

# define pallets

# define common color
colour_yellow = "#EEb462"
colour_purple = "#534666"
colour_mint = "#138086"
colour_grey = "grey55"
colour_coral = "#CD7672"
  


# functions ---------------------------------------------------------------

# set up model function
brm_logistic <- function(model_formula) {
  
  brm(bf(model_formula), 
      family = "bernoulli", 
      data = dat_merged, 
      seed = 1708, 
      cores = parallel::detectCores(), 
      chains = 4, 
      iter = 10000, 
      silent = 2,
      refresh = 0)
  
}

