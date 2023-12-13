library(igraph)
library(ggplot2)

args <- commandArgs()
wd <- args[6]
setwd(wd)

ls <- list.files(pattern = "\\.txt$")

merged_network <- NULL

num_random_networks = 1000

# Looping through the edgelists and merging the networks
for (file in ls) {
  edge <- read.table(file, header = FALSE)
  edges <- edge[c(2, 4)]
  r1 <- as.data.frame(edges)
  graph <- graph_from_data_frame(d = r1, directed = FALSE)
  
  if (is.null(merged_network)) {
    merged_network <- graph
  } else {
    merged_network <- intersection(merged_network, graph)
  }

# Generating weighted edgelist for the merged network
weighted_edgelist <- get.data.frame(merged_network, what = "edges") 
colnames(weighted_edgelist) = c("from", "to")
edge_counts <- table(paste(weighted_edgelist$from, weighted_edgelist$to))
merged_weighted_edgelist <- cbind(weighted_edgelist, weight = edge_counts[match(paste(weighted_edgelist$from, weighted_edgelist$to), names(edge_counts))])
write.csv(merged_weighted_edgelist, file = "Merged_Network_Weighted_Edgelist_1.csv" , row.names = FALSE)

count1 <- vcount(merged_network)
count2 <- ecount(merged_network)

#Obtaining degree frequency distribution for merged graph
degrees_merged <- degree(merged_network)
degree_dist_merged <- data.frame(degree = degrees_merged)
degree_dist_merged <- as.data.frame(table(degree_dist_merged))
colnames(degree_dist_merged) <- c("Degree", "Frequency")
degree_dist_merged$Degree <- as.numeric(degree_dist_merged$Degree)
degree_dist_merged$log_degree <- log(degree_dist_merged$Degree)
degree_dist_merged$log_frequency <- log(degree_dist_merged$Frequency)

#Obtaining hubs from merged graph
sorted_nodes <- sort(degrees_merged, decreasing = TRUE)
hubs_merged <- sorted_nodes[1:100]
write.table(hubs_merged, file="merged_hubs.txt", col.names = FALSE)

#Fitting a linear model to the data for ecotype graphs
model_merged <- lm(log_frequency ~ log_degree, data = degree_dist_merged)
coefficients_merged <- coef(model_merged)
power_law_exponent_merged <- coefficients_merged["log_degree"]
power_law_curve_merged <- data.frame(log_degree = degree_dist_merged$log_degree,
                              predicted_log_frequency_merged = predict(model_merged, degree_dist_merged))

# Plotting the log-log graph and the power-law curve for ecotype graphs
merged_graph <- ggplot(data = degree_dist_merged, aes(x = log_degree, y = log_frequency)) +
  geom_point() +
  geom_line(data = power_law_curve_merged, aes(x = log_degree, y = predicted_log_frequency_merged), color = "red") +
  labs(title = paste("Power Law Fit (alpha =", round(power_law_exponent_merged, 2), ")"),
       x = "Log(Degree)",
       y = "Log(Frequency)") +
  theme_minimal()

svg("merged_graph.svg")
plot(merged_graph)
dev.off()

#Plotting degree distribution graph (L-graph) for ecotype graphs
hist <- as.data.frame(table(degrees_merged))
hist[,1] <- as.numeric(paste(hist[,1]))
L_g <- ggplot(hist, aes(x = degrees_merged, y = Freq)) + geom_point() +
scale_x_continuous("Degree") +
scale_y_continuous("Frequency") + ggtitle("Degree Distribution") + theme_bw()

svg("Merged_lg.svg")
plot(L_g)
dev.off()

networks <- list()

for (x in 1:num_random_networks) {
graphr_merged <- erdos.renyi.game(count1, count2, type = "gnm")

#Obtaining degree frequency distribution for random graphs
degrees_r_merged <- degree(graphr_merged)
degree_dist_r_merged <- data.frame(degree = degrees_r_merged)
degree_dist_r_merged <- as.data.frame(table(degree_dist_r_merged))
colnames(degree_dist_r_merged) <- c("Degree", "Frequency")
degree_dist_r_merged$Degree <- as.numeric(degree_dist_r_merged$Degree)
degree_dist_r_merged$log_degree <- log(degree_dist_r_merged$Degree)
degree_dist_r_merged$log_frequency <- log(degree_dist_r_merged$Frequency)

#Fitting a linear model to the data for random graphs
model_r_merged <- lm(log_frequency ~ log_degree, data = degree_dist_r_merged)
coefficients_r_merged <- coef(model_r_merged)
power_law_exponent_r_merged <- coefficients_r_merged["log_degree"]
power_law_curve_r_merged <- data.frame(log_degree = degree_dist_r_merged$log_degree,
                              predicted_log_frequency_r_merged = predict(model_r_merged, degree_dist_r_merged))

# Plotting the log-log graph and the power-law curve for random graphs
random_merged <- ggplot(data = degree_dist_r_merged, aes(x = log_degree, y = log_frequency)) +
  geom_point() +
  geom_line(data = power_law_curve_r_merged, aes(x = log_degree, y = predicted_log_frequency_r_merged), color = "red") +
  labs(title = paste("Power Law Fit (alpha =", round(power_law_exponent_r_merged, 2), ")"),
       x = "Log(Degree)",
       y = "Log(Frequency)") +
  theme_minimal()

svg(paste(x, 'random_merged.svg', sep ="_"))
plot(random_merged)
dev.off()

#Plotting degree distribution graph (gauss-graph) for random graphs
degrees_gr <- degree(graphr_merged)
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

svg(paste(x, 'merged_si.svg', sep ="_"))
plot(gauss_random)
dev.off()

}

