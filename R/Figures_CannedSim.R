# 
# 
# 
# library(ggplot2)
# library(cowplot)
# library(reshape2)
# 
# # Create example data
# set.seed(123)
# n <- 100
# L1 <- matrix(runif(n*n), nrow = n)
# 
# 
# # Convert matrices to data frames for ggplot
# L1_df <- melt(L1)
# colnames(L1_df) <- c("x", "y", "value")
# # Plot L1
# ggplot(L1_df, aes(x = x, y = y, fill = value)) + 
#   geom_tile() + 
#   scale_fill_gradient(low = "blue", high = "white") + 
#   ggtitle("L1") +
#   theme_minimal() +
#   theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
# 
# 
# 
# 
# 
# 
# Receptor_Lingand_Sim <- function(save_path = NULL, width = 12, height = 5){
#   # Load required libraries
#   library(ggplot2)
#   library(cowplot)
#   library(reshape2)
#   
#   # Set dimensions and parameters
#   n <- 100
#   circle_radius <- 10
#   donut_inner_radius <- 15
#   donut_outer_radius <- 20
#   circle1_center <- c(75, 75)
#   circle2_center <- c(25, 25)
#   
#   # Generate random background
#   background <- matrix(runif(n*n), nrow = n)
#   
#   # Function to create circles
#   create_circle <- function(center, radius, n) {
#     circle <- matrix(0, nrow = n, ncol = n)
#     for (i in 1:n) {
#       for (j in 1:n) {
#         if ((i - center[1])^2 + (j - center[2])^2 <= radius^2) {
#           circle[i, j] <- 1
#         }
#       }
#     }
#     return(circle)
#   }
#   
#   # Create circles for L1
#   circle1 <- create_circle(circle1_center, circle_radius, n)
#   circle2 <- create_circle(circle2_center, circle_radius, n)
#   L1 <- background
#   L1[circle1 == 1 | circle2 == 1] <- 1
#   
#   # Function to create donuts
#   create_donut <- function(center, inner_radius, outer_radius, n) {
#     donut <- matrix(0, nrow = n, ncol = n)
#     for (i in 1:n) {
#       for (j in 1:n) {
#         dist_squared <- (i - center[1])^2 + (j - center[2])^2
#         if (dist_squared <= outer_radius^2 && dist_squared >= inner_radius^2) {
#           donut[i, j] <- 1
#         }
#       }
#     }
#     return(donut)
#   }
#   
#   # Create donuts for R
#   donut1 <- create_donut(circle1_center, donut_inner_radius, donut_outer_radius, n)
#   donut2 <- create_donut(circle2_center, donut_inner_radius, donut_outer_radius, n)
#   R <- background
#   R[donut1 == 1 | donut2 == 1] <- 1
#   
#   # Compute difference L1-R
#   diff <- L1 - R
#   
#   # Convert matrices to data frames for ggplot
#   L1_df <- melt(L1)
#   colnames(L1_df) <- c("x", "y", "value")
#   R_df <- melt(R)
#   colnames(R_df) <- c("x", "y", "value")
#   diff_df <- melt(diff)
#   colnames(diff_df) <- c("x", "y", "value")
#   
#   # Plot L1
#   p1 <- ggplot(L1_df, aes(x = x, y = y, fill = value)) + 
#     geom_tile() + 
#     scale_fill_gradient(low = "blue", high = "green") + 
#     theme_nothing()
# 
#   
#   p1
#   
#   # Plot R
#   p2 <- ggplot(R_df, aes(x = x, y = y, fill = value)) + 
#     geom_tile() + 
#     scale_fill_gradient(low = "blue", high = "green") + 
#     theme_nothing()
# 
#   
#   # p2
#   
#   # Plot L1-R
#   p3 <- ggplot(diff_df, aes(x = x, y = y, fill = value)) + 
#     geom_tile() + 
#     scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + 
#     ggtitle("L1-R") +
#     theme_minimal() +
#     theme_nothing()
#   
#   # p3
#   
#   # Arrange plots
#   final_plot <- plot_grid(p1, p2, p3, nrow = 1)
#   
#   print(final_plot)
#   
#   # Save plot
#   if(!is.null(save_path)) ggsave(save_path, final_plot, width = width, height = height) 
# }
# 
# # Receptor_Lingand_Sim(save_path = "./Receptor_Ligand_sim_pdf.pdf")
