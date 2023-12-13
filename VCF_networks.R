library(igraph)
library(ggplot2)

args <- commandArgs()
wd <- args[6]
setwd(wd)

ls <- list.files(pattern="edgelist.txt")

num_random_networks = 1000

for(i in ls) {
edge <- read.table(i, header = FALSE, col.names = c("TF_ID", "from", "Gene_ID", "to", "Coord"))
edges <- edge[c(2, 4)]
l1 <- as.data.frame(edges)
graph <- graph_from_data_frame(l1)
l2 <- vcount(graph)
l3 <- ecount(graph)

#Obtaining degree frequency distribution for ecotype graphs
degrees <- degree(graph)
degree_dist <- data.frame(degree = degrees)
degree_dist <- as.data.frame(table(degree_dist))
colnames(degree_dist) <- c("Degree", "Frequency")
degree_dist$Degree <- as.numeric(degree_dist$Degree)
degree_dist$log_degree <- log(degree_dist$Degree)
degree_dist$log_frequency <- log(degree_dist$Frequency)

#Obtaining hubs from ecotype graphs
sorted_nodes <- sort(degrees, decreasing = TRUE)
hubs <- sorted_nodes[1:100]
write.table(hubs, file = paste(i,"hubs.txt", sep = "_"), col.names = FALSE)

#Obtaining weighted edgelists from ecotype graphs
edge_counts <- table(paste(edge$from, edge$to))
weighted_edgelist <- cbind(edge, weight = edge_counts[match(paste(edge$from, edge$to), names(edge_counts))])
write.csv(weighted_edgelist, file = paste(i,"weighted_edgelist.csv", sep = "_"), row.names = FALSE)

#Fitting a linear model to the data for ecotype graphs
model <- lm(log_frequency ~ log_degree, data = degree_dist)
coefficients <- coef(model)
power_law_exponent <- coefficients["log_degree"]
power_law_curve <- data.frame(log_degree = degree_dist$log_degree,
                              predicted_log_frequency = predict(model, degree_dist))

# Plotting the log-log graph and the power-law curve for ecotype graphs
l4 <- ggplot(data = degree_dist, aes(x = log_degree, y = log_frequency)) +
  geom_point() +
  geom_line(data = power_law_curve, aes(x = log_degree, y = predicted_log_frequency), color = "red") +
  labs(title = paste("Power Law Fit (alpha =", round(power_law_exponent, 2), ")"),
       x = "Log(Degree)",
       y = "Log(Frequency)") +
  theme_minimal()

svg(paste(i,"gr_svg", sep="_"))
plot(l4)
dev.off()

#Plotting degree distribution graph (L-graph) for ecotype graphs
hist <- as.data.frame(table(degrees))
hist[,1] <- as.numeric(paste(hist[,1]))
L_g <- ggplot(hist, aes(x = degrees, y = Freq)) + geom_point() +
scale_x_continuous("Degree") +
scale_y_continuous("Frequency") + ggtitle("Degree Distribution") + theme_bw()

svg(paste(i, "lg.svg", sep="_"))
plot(L_g)
dev.off()

networks <- list()

for (x in 1:num_random_networks) {
graphr <- erdos.renyi.game(l2, l3, type = "gnm")

#Obtaining degree frequency distribution for random graphs
degrees_r <- degree(graphr)
degree_dist_r <- data.frame(degree = degrees_r)
degree_dist_r <- as.data.frame(table(degree_dist_r))
colnames(degree_dist_r) <- c("Degree", "Frequency")
degree_dist_r$Degree <- as.numeric(degree_dist_r$Degree)
degree_dist_r$log_degree <- log(degree_dist_r$Degree)
degree_dist_r$log_frequency <- log(degree_dist_r$Frequency)

#Fitting a linear model to the data for random graphs
model_r <- lm(log_frequency ~ log_degree, data = degree_dist_r)
coefficients_r <- coef(model_r)
power_law_exponent_r <- coefficients_r["log_degree"]
power_law_curve_r <- data.frame(log_degree = degree_dist_r$log_degree,
                              predicted_log_frequency_r = predict(model_r, degree_dist_r))

# Plotting the log-log graph and the power-law curve for random graphs
l5 <- ggplot(data = degree_dist_r, aes(x = log_degree, y = log_frequency)) +
  geom_point() +
  geom_line(data = power_law_curve_r, aes(x = log_degree, y = predicted_log_frequency_r), color = "red") +
  labs(title = paste("Power Law Fit (alpha =", round(power_law_exponent_r, 2), ")"),
       x = "Log(Degree)",
       y = "Log(Frequency)") +
  theme_minimal()

svg(paste(i, x, '_random.svg', sep ="_"))
plot(l5)
dev.off()

#Plotting degree distribution graph (gauss-graph) for random graphs
degrees_gr <- degree(graphr)
hist_gr <- as.data.frame(table(degrees_gr))
hist_gr[,1] <- as.numeric(paste(hist_gr[,1]))
networks[[x]] <- hist_gr 
}

# Create an empty ggplot object
gauss_random <- ggplot() + labs(title = paste("Superimposed Random Networks (n=",num_random_networks,")", sep = ""), x = "Degree", y = "Frequency") + theme_minimal()

# Superimpose the graphs using ggplot2
for (x in 1:num_random_networks) {
gauss_random <- gauss_random + geom_line(data = networks[[x]], aes(x = degrees_gr, y = Freq), color = rainbow(num_random_networks)[x])
}

svg(paste(i, x, '_si.svg', sep ="_"))
plot(gauss_random)
dev.off()


matrix <- as.matrix(get.adjacency(graph.data.frame(edges)))
write.table(matrix,file=paste0(i,"_matrix.csv"),sep=",")

}



