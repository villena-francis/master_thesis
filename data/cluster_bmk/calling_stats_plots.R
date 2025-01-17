################################################################################
# COMPUTATIONAL ANALYSIS OF STRUCTURAL VARIANT CALLERS
# This script analyzes performance metrics and computational resource usage
# of different SV callers at various sequencing depths
################################################################################

# Load required libraries
library(dplyr)        # For data manipulation
library(tidyr)        # For data reshaping
library(ggplot2)      # For creating plots
library(lubridate)    # For time calculations
library(colorspace)   # For color manipulation
library(gridExtra)    # For combining plots
library(grid)         # For advanced graphics
library(patchwork)    # For plot composition

################################################################################
# SECTION 1: VISOR ANALYSIS
################################################################################

# Read VISOR performance data
visor_data <- read.csv("hpc_data/hpc_visor_laser.csv")
visor_data$Sequencing_deph <- factor(visor_data$Sequencing_deph)

# Time conversion function: converts HH:MM:SS to days
convert_time_to_days <- function(time_str) {
  # Separate days if they exist
  parts <- strsplit(time_str, "-")[[1]]
  
  if (length(parts) == 2) {
    days <- as.numeric(parts[1])
    time <- parts[2]
  } else {
    days <- 0
    time <- parts[1]
  }
  
  # Convert HH:MM:SS format
  time_components <- as.numeric(strsplit(time, ":")[[1]])
  hours <- time_components[1]
  minutes <- time_components[2]
  seconds <- time_components[3]
  
  # Convert everything to days
  total_days <- days + 
    hours/24 + 
    minutes/(24*60) + 
    seconds/(24*60*60)
  
  return(round(total_days, 3))
}

# Process time data and calculate CPU cores
visor_data$Time <- sapply(visor_data$Time, convert_time_to_days)
visor_data <- visor_data %>%
  mutate(CPU_cores = CPU_asigned * (CPU_usage / 100))

# Define color schemes for visualization
fill_colors <- c(
  Time = "#D0D8DB",      # Light shade derivatives
  CPU_cores = "#B0C7C9", 
  RAM_usage = "#A4C5CC"  
)

point_colors <- c(
  Time = "#90A3A8",      # Base colors for points
  CPU_cores = "#496F73", 
  RAM_usage = "#5E9CA3"  
)

# Calculate statistics for plotting
avg_CPU_asigned <- mean(visor_data$CPU_asigned)
global_y_max <- max(visor_data$CPU_cores)

# Create CPU usage plot
visor_cpu <- ggplot(visor_data, aes(x = Sequencing_deph, y = CPU_cores)) +
  geom_boxplot(
    fill = fill_colors["CPU_cores"],
    color = point_colors["CPU_cores"],
    alpha = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    color = point_colors["CPU_cores"],
    width = 0.2,
    height = 0,
    alpha = 0.5
  ) +
  geom_hline(
    yintercept = avg_CPU_asigned,
    linetype = "dashed",
    color = "darkgray"
  ) +
  geom_text(
    aes(
      x = Inf,
      y = avg_CPU_asigned,
      label = paste("asigned =", round(avg_CPU_asigned, 2))
    ),
    hjust = 1.1,
    vjust = ifelse(avg_CPU_asigned > (global_y_max * 0.9), 1.5, -0.5),
    color = "darkgray",
    size = 3
  ) +
  scale_x_discrete(breaks = c(30, 50, 100, 200)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, NA)) +
  labs(
    x = "Sequencing Depth",
    y = "CPU Cores"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9)
  )

# Create memory usage plot
visor_mem <- ggplot(visor_data, aes(x = Sequencing_deph, y = RAM_usage)) +
  geom_boxplot(
    fill = fill_colors["RAM_usage"],
    color = point_colors["RAM_usage"],
    alpha = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    color = point_colors["RAM_usage"],
    width = 0.2,
    height = 0,
    alpha = 0.5
  ) +
  scale_x_discrete(breaks = c(30, 50, 100, 200)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, NA)) +
  labs(
    x = "Sequencing Depth",
    y = "RAM Usage (GB)"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9)
  )

# Create execution time plot
global_y_max <- max(visor_data$Time)

visor_times <- ggplot(visor_data, aes(x = Sequencing_deph, y = Time)) +
  geom_boxplot(
    fill = fill_colors["Time"],
    color = point_colors["Time"],
    alpha = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    color = point_colors["Time"],
    width = 0.2,
    height = 0,
    alpha = 0.5
  ) +
  scale_x_discrete(breaks = c(30, 50, 100, 200)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, NA)) +
  labs(
    x = "Sequencing Depth",
    y = "Time (days)"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9)
  )

# Combine VISOR plots
visor_combined_plot <- (visor_cpu + plot_spacer() + visor_mem + plot_spacer() + visor_times) +
  plot_layout(
    widths = c(1, 0.1, 1, 0.1, 1),
    guides = "collect"
  ) &
  theme(aspect.ratio = 1)

print(visor_combined_plot)

################################################################################
# SECTION 2: SV CALLERS ANALYSIS
################################################################################

# Read SV callers performance data
hpc_data <- read.csv("hpc_data/hpc_callers.csv")
hpc_data$Sequencing_deph <- factor(hpc_data$Sequencing_deph)

# Convert time format and calculate CPU cores
hpc_data$Time <- as.numeric(hms(hpc_data$Time)) / 60
hpc_data <- hpc_data %>%
  mutate(CPU_cores = CPU_asigned * (CPU_usage / 100))

# Transform data to long format
hpc_data_long <- hpc_data %>%
  pivot_longer(
    cols = c(CPU_cores, RAM_usage, Time),
    names_to = "Metric", 
    values_to = "Value"
  )

# Calculate CPU statistics
cpu_assigned_avg <- hpc_data %>%
  group_by(Caller) %>%
  summarise(avg_CPU_asigned = mean(CPU_asigned, na.rm = TRUE))

# Join average CPU data
hpc_data_long <- hpc_data_long %>%
  left_join(cpu_assigned_avg, by = "Caller")

# Calculate global maximum for CPU cores
global_y_max <- hpc_data_long %>%
  filter(Metric == "CPU_cores") %>%
  summarise(global_y_max = max(Value, na.rm = TRUE)) %>%
  pull(global_y_max)

# Create SV callers performance plot
hpc_plot <- ggplot(hpc_data_long, aes(x = Sequencing_deph, y = Value)) +
  geom_boxplot(
    aes(fill = Metric, color = Metric),
    alpha = 0.5
  ) +
  geom_jitter(
    aes(color = Metric),
    width = 0.2,
    height = 0,
    alpha = 0.5
  ) +
  geom_hline(
    data = subset(hpc_data_long, Metric == "CPU_cores"), 
    aes(yintercept = avg_CPU_asigned),
    linetype = "dashed",
    color = "darkgray"
  ) +
  geom_text(
    data = subset(hpc_data_long, Metric == "CPU_cores"),
    aes(
      x = Inf,
      y = avg_CPU_asigned,
      label = paste("asigned =", round(avg_CPU_asigned, 2)),
      vjust = ifelse(avg_CPU_asigned > (global_y_max * 0.9), 1.5, -0.5)
    ),
    hjust = 1.1,
    color = "darkgray",
    size = 3
  ) +
  facet_grid(rows = vars(Metric), cols = vars(Caller), scales = "free_y") +
  scale_x_discrete(breaks = c(30, 50, 100, 200)) +
  labs(y = "") +
  theme(strip.text = element_text(face = "bold")) +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(
    values = c(point_colors, sapply(fill_colors, function(color) darken(color, 0.2)))
  )

print(hpc_plot)

################################################################################
# SECTION 3: CALLER ACCURACY ANALYSIS
################################################################################

# Read accuracy metrics data
calls_data <- read.csv("calls_data/callings.csv")
calls_data$Sequencing_deph <- factor(calls_data$Sequencing_deph)

# Transform to long format
calls_data_long <- calls_data %>%
  pivot_longer(
    cols = c(Precision, Recall, F1_score),
    names_to = "Metric",
    values_to = "Value"
  )

# Order metric levels
calls_data_long$Metric <- factor(
  calls_data_long$Metric, 
  levels = c("Precision", "Recall", "F1_score")
)

# Define colors for accuracy metrics
fill_colors <- c(
  Recall = "#E6F3FF",
  F1_score = "#B3D7DC",
  Precision = "#E6F0E6"
)

point_colors <- c(
  Recall = "#193C40",
  F1_score = "#539DA6",
  Precision = "#A8B5BF"
)

# Create accuracy metrics plot
calls_plot <- ggplot(calls_data_long, aes(x = Sequencing_deph, y = Value)) +
  geom_boxplot(
    aes(fill = Metric, color = Metric),
    alpha = 0.5
  ) +
  geom_jitter(
    aes(color = Metric),
    width = 0.2,
    height = 0,
    alpha = 0.5
  ) +
  facet_grid(rows = vars(Metric), cols = vars(Caller), scales = "free_y") +
  scale_x_discrete(breaks = c(30, 50, 100, 200)) +
  labs(y = "") +
  theme(strip.text = element_text(face = "bold")) +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = point_colors)

print(calls_plot)
