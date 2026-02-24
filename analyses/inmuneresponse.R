response_data_raw = read.csv("data/flies_data.csv",sep = ";")
key = read.csv("data/flies_info.csv",sep = ";")

my_color_palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")

response_data_pre <- response_data_pre |>
  dplyr::mutate(treatment = factor(treatment, levels = c("T1", "T2", "T3", "T4")))

response_data_pre <- response_data_raw |>
  tidyr::pivot_longer(
    cols = -c(Date, Time, Light1_Dark0),
    names_to = "FlyID",
    values_to = "Activity"
  ) |>
  dplyr::mutate(FlyID = stringr::str_remove(FlyID, "^.")) |>
  dplyr::group_by(FlyID) |>
  dplyr::mutate(Minute = dplyr::row_number()) |>
  dplyr::ungroup() |>
  dplyr::left_join(
    key |> dplyr::select(Monitor_TubeLoc, Fly_Sex, Treatment, Status_dead0_alive1),
    by = dplyr::join_by(FlyID == Monitor_TubeLoc)
  ) |> dplyr::filter(Treatment != "BLANK") |> remove_trailing_zeros() |> 
  dplyr::group_by(FlyID) |>
  dplyr::mutate(Last_mov = length(Activity)) |>
  dplyr::ungroup() |>
  dplyr::mutate(Treatment = factor(Treatment, levels = c("naive", "heatkill", "PBS", "Pentomophila")))

response_data_pre |>
  dplyr::group_by(FlyID) |>
  dplyr::summarise(
    last_mov  = dplyr::first(Last_mov) / 60,
    treatment = dplyr::first(Treatment),
    Fly_Sex   = dplyr::first(Fly_Sex),
    .groups   = "drop"
  ) |>
  ggplot(aes(x = treatment, y = last_mov, fill = treatment)) +
  geom_boxplot(
    outlier.shape = NA, alpha = 0.8, na.rm = TRUE,
    position = position_dodge()
  ) +
  geom_point(
    na.rm = TRUE,
    position = position_jitterdodge(0.27),
    pch = 21
  ) +
  facet_grid(rows = vars(Fly_Sex)) +
  labs(
    x = "Treatment",
    y = "Last Movement Time (hours)"
  ) +
  scale_fill_manual(values = my_color_palette) +
  theme_light() +
  theme(
    legend.position = "none",
    strip.text      = element_text(face = "bold"),
    axis.title      = element_text(size = 12),  # axis labels
    axis.text       = element_text(size = 10)   # tick labels
  )






























fly_data_hourly = response_data |>
  mutate(hour = floor(time / 60), Status = paste0(age,"/",treatment)) |> filter(hour<2) |>
  group_by(FlyID, Status, hour) |>
  summarise(activity = mean(activity, na.rm = TRUE), .groups = "drop") 

mean_data_hourly = fly_data_hourly |>
  group_by(Status, hour) |>
  summarise(mean_activity = mean(activity), .groups = "drop")

p1 = ggplot() +
  geom_line(data = fly_data_hourly, aes(x = hour, y = activity, group = FlyID, color = Status), alpha = 0.1) +
  geom_line(data = mean_data_hourly, aes(x = hour, y = mean_activity, color = "black"), size = 1.0) +
  facet_wrap(~Status, ncol = 2) +  # Separate panels for each Status
  labs(x = "hour",
       y = "mean activity",
       color = "Status") +
  theme_minimal() +
  theme(legend.position = "none")#+
  #scale_color_manual(values = c("infected" = "steelblue", "dead" = "darkred")) 

p1

ggplot() +
  geom_line(data = response_data, aes(x = time, y = activity, group = FlyID), alpha = 0.1) +
  #facet_wrap(~Status, ncol = 2) +  # Separate panels for each Status
  labs(x = "hour",
       y = "mean activity",
       color = "Status") +
  theme_minimal() +
  theme(legend.position = "none")

fly_data_response = split(response_data,response_data$FlyID)

############ bout analyses ############

calculate_prop_active <- function(activity_data, threshold = 0) {
  (sum(activity_data > threshold) / length(activity_data))
}

# Apply to each fly in the dataset
fly_data_response <- lapply(fly_data_response, function(fly_df) {
  fly_df$prop_active <- calculate_prop_active(fly_df$activity)
  return(fly_df)
})

fly_data_features <- do.call(rbind,fly_data_response) |> 
  group_by(FlyID, age, sex, treatment) |>
  mutate(
    
    # Bout detection (active = activity > threshold)
    state = ifelse(activity > 0, "active", "inactive"),
    bout_id = data.table::rleid(state),
    
    # Bout features
    bout_length = sequence(rle(state)$lengths)
  ) |>
  ungroup()

transition_stats <- fly_data_features |>
  group_by(FlyID, age, sex, treatment) |>
  summarise(
    # Bout duration distributions (using median for non-normal data)
    mean_active_bout = mean(bout_length[state == "active"], na.rm = TRUE),
    mean_inactive_bout = mean(bout_length[state == "inactive"], na.rm = TRUE),
    
    # activity rate
    activity_rate = sum(activity[state == "active"], na.rm = TRUE),#/sum(state == "active"),
    
    # Proportion of time active (should be constant per fly)
    prop_active = first(prop_active),
    
    # Transition matrix
    #transition_prob = log(as.numeric(markovchainFit(data = state)$estimate[1])[2]/as.numeric(markovchainFit(data = state)$estimate[2])[1]),
    .groups = "drop" # Important for preventing downstream issues
  )

# --- Gamma GLM for Active Bout Duration ---
# The formula tests all main effects and their interactions (*)
glm_active_bout <- glm(
  mean_active_bout ~ treatment * sex * age,
  data = transition_stats,
  family = Gamma(link = "log") # Use Gamma for skewed, positive data
)

Anova(glm_active_bout, type = "II")

# Use emmeans to compare treatment effects, separately for each sex
posthoc_results <- emmeans(glm_prop_active, pairwise ~ age | treatment)

# View the results
print(posthoc_results)

# --- Binomial GLM for Proportion of Time Active ---
# First, you need to calculate the number of active and total minutes per fly.
# Let's assume your experiment ran for 14 days (20160 minutes).
total_minutes <- 200

transition_stats <- transition_stats |>
  mutate(
    active_minutes = round(prop_active * total_minutes),
    inactive_minutes = total_minutes - active_minutes
  )

# Build the model using the cbind() syntax for successes/failures
glm_prop_active <- glm(
  cbind(active_minutes, inactive_minutes) ~ treatment * sex * age,
  data = transition_stats,
  family = quasibinomial(link = "logit")
)

Anova(glm_prop_active, type = "II")

# Use emmeans to compare treatment effects, separately for each sex
posthoc_results <- emmeans(glm_prop_active, pairwise ~ age | treatment)

# View the results
print(posthoc_results)

# --- Gamma GLM for activity rate ---

glm_activity_rate <- glm(
  activity_rate ~ treatment * sex * age,
  data = transition_stats,
  family = Gamma(link = "log")
)

Anova(glm_activity_rate, type = "II")

# Use emmeans to compare treatment effects, separately for each sex
posthoc_results <- emmeans(glm_activity_rate, pairwise ~ sex | treatment)

# View the results
print(posthoc_results)

########## transition probs #######

glm_transition <- glm(
  transition_prob ~ treatment * sex * age,
  data = transition_stats,
  family = quasibinomial(link = "logit")
)

# Get the significance tests for each factor and their interactions
Anova(glm_transition, type = "II")

# Use emmeans to compare treatment effects, separately for each sex
#posthoc_results <- emmeans(glm_transition, pairwise ~ sex | age)

# View the results
#print(posthoc_results)

# Use emmeans to compare treatment effects, separately for each sex
#posthoc_results <- emmeans(glm_transition, pairwise ~ treatment | age)

# View the results
#print(posthoc_results)

############# boxplots ###########

ggplot(transition_stats, aes(x = treatment, y = mean_active_bout, fill = treatment)) + geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, height = 0) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 4, fill = "white") +
  facet_grid(sex~age) +
  labs(
    x = "age",
    y = "mean active bout duration (min)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


ggplot(transition_stats, aes(x = treatment, y = prop_active, fill = treatment)) + geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, height = 0) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 4, fill = "white") +
  facet_grid(sex~age) +
  labs(
    x = "age",
    y = "mean active bout duration (min)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(transition_stats, aes(x = treatment, y = activity_rate, fill = treatment)) + geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, height = 0) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 4, fill = "white") +
  facet_grid(sex~age) +
  labs(
    x = "age",
    y = "activity rate"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

######### barplot #########

transition_stats |>
  group_by(age, sex, treatment) |>
  summarise(
    total_activity = mean(activity_rate, na.rm = TRUE),
    se_activity = sd(activity_rate, na.rm = TRUE) / sqrt(n())
  ) |>
  ungroup() |>
  ggplot(aes(x = treatment, y = total_activity, fill = treatment)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = total_activity - se_activity, ymax = total_activity + se_activity),
    position = position_dodge(width = 0.8), width = 0.2
  ) +
  facet_grid(sex ~ age) +
  labs(
    x = "treatment",
    y = "total Activity",
    fill = "Treatment"
  ) +
  scale_fill_manual(values = c("UC" = "#F8766D", "PBS" = "#7CAE00", "HK" = "#00BFC4", "P.ento" = "#C77CFF")) +
  theme_minimal() +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

transition_stats |>
  group_by(age, sex, treatment) |>
  summarise(
    total_activity = mean(prop_active, na.rm = TRUE),
    se_activity = sd(prop_active, na.rm = TRUE) / sqrt(n())
  ) |>
  ungroup() |>
  ggplot(aes(x = treatment, y = total_activity, fill = treatment)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = total_activity - se_activity, ymax = total_activity + se_activity),
    position = position_dodge(width = 0.8), width = 0.2
  ) +
  facet_grid(sex ~ age) +
  labs(
    x = "treatment",
    y = "proportion of time active",
    fill = "Treatment"
  ) +
  scale_fill_manual(values = c("UC" = "#F8766D", "PBS" = "#7CAE00", "HK" = "#00BFC4", "P.ento" = "#C77CFF")) +
  theme_minimal() +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

transition_stats |>
  group_by(age, sex, treatment) |>
  summarise(
    total_activity = mean(mean_active_bout, na.rm = TRUE),
    se_activity = sd(mean_active_bout, na.rm = TRUE) / sqrt(n())
  ) |>
  ungroup() |>
  ggplot(aes(x = treatment, y = total_activity, fill = treatment)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = total_activity - se_activity, ymax = total_activity + se_activity),
    position = position_dodge(width = 0.8), width = 0.2
  ) +
  facet_grid(sex ~ age) +
  labs(
    x = "treatment",
    y = "mean active bout",
    fill = "Treatment"
  ) +
  scale_fill_manual(values = c("UC" = "#F8766D", "PBS" = "#7CAE00", "HK" = "#00BFC4", "P.ento" = "#C77CFF")) +
  theme_minimal() +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

#ggplot(transition_stats, aes(x = age, y = transition_prob, fill = age)) + geom_boxplot(alpha = 0.7, outlier.shape = NA) +
 # geom_jitter(width = 0.15, alpha = 0.6, height = 0) +
#  stat_summary(fun = mean, geom = "point", shape = 21, size = 4, fill = "white") +
 # labs(
#    x = "age",
 #   y = "inactivity/activity log ratio"
  #) +
  #theme_minimal() +
  #theme(legend.position = "none")

########## correlations ########

# --- 3. Calculate Kendall's Tau For Each Group ---

# Helper function to calculate the mean of the off-diagonal elements of a matrix.
# This gives us a single "average synchronization" score for the group.
mean_off_diagonal <- function(mat) {
  # Set the diagonal (self-correlation) to NA so it's ignored
  diag(mat) <- NA
  # Calculate the mean of all other pairwise correlations
  mean(mat, na.rm = TRUE)
}

# The main workflow
group_correlation_summary <- do.call(rbind,fly_data_response) |>
  # Group by your experimental conditions
  group_by(sex, age, treatment) |>
  
  # Nest the data: creates a list-column 'data' containing the time series for each group
  nest() |>
  
  # Create a new column 'corr_matrix' by applying the correlation function to each nested dataset
  mutate(
    corr_matrix = map(data, function(group_df) {
      
      # For correlation, data must be in "wide" format (flies as columns)
      wide_df <- group_df |>
        pivot_wider(names_from = FlyID, values_from = activity, id_cols = time) |>
        dplyr::select(-time) # Remove the time column
      
      # Calculate the Kendall's tau correlation matrix for all flies in the group
      cor(wide_df, method = "kendall", use = "pairwise.complete.obs")
    })
  ) |>
  
  # Create the final summary metric by applying our helper function to each matrix
  mutate(
    mean_kendall_tau = map_dbl(corr_matrix, mean_off_diagonal)
  )

# --- 4. Display the Final Summary Matrix ---
# Select only the relevant columns to create your final results table
final_summary_matrix <- group_correlation_summary |>
  dplyr::select(sex, age, treatment, mean_kendall_tau)

print("Summary Matrix of Average Group Synchronization (Kendall's Tau):")
print(final_summary_matrix)

ggplot(final_summary_matrix, aes(x = treatment, y = mean_kendall_tau, fill = age)) +
  geom_col(position = "dodge") +
  facet_wrap(~sex) + 
  labs(
    x = "Treatment",
    y = "Mean Kendall's Tau"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


######## Phase locking value ########

# Align all flies to the shortest length in the group
align_time_series <- function(fly_list) {
  min_len <- min(sapply(fly_list, nrow))
  lapply(fly_list, function(df) df[1:min_len, ])
}

compute_plv <- function(fly1_df, fly2_df) {
  n <- min(nrow(fly1_df), nrow(fly2_df))
  if (n == 0) return(NA_real_)
  t1 <- as.matrix(fly1_df[1:n, c("time", "activity")])
  t2 <- as.matrix(fly2_df[1:n, c("time", "activity")])
  plv_result <- phase.sync(t1, t2, nrands = 0, mod = 1, nbreaks = 100, method = "fft", quiet = TRUE)
  abs(mean(exp(1i * plv_result$deltaphase$phasediff), na.rm = TRUE))
}

is_valid_fly <- function(df, min_length = 100) {
  nrow(df) >= min_length &&
    !all(is.na(df$activity)) &&
    length(unique(df$activity)) > 1 &&
    !all(df$activity == 0)
}

compute_group_plv_matrix <- function(group_flies, min_length = 100) {
  # Filter out flies with too few data points
  group_flies <- Filter(function(df) nrow(df) >= min_length, group_flies)
  # Align all to the shortest length
  group_flies <- align_time_series(group_flies)
  n <- length(group_flies)
  if (n < 2) return(matrix(NA, nrow = n, ncol = n))
  plv_mat <- matrix(NA, nrow = n, ncol = n)
  rownames(plv_mat) <- names(group_flies)
  colnames(plv_mat) <- names(group_flies)
  for (i in 1:n) {
    for (j in i:n) {
      if (i == j) {
        plv_mat[i, j] <- 1
      } else {
        plv <- compute_plv(group_flies[[i]], group_flies[[j]])
        plv_mat[i, j] <- plv
        plv_mat[j, i] <- plv
      }
    }
  }
  plv_mat
}

# Extract metadata for all flies
fly_metadata <- map_dfr(fly_data_response, ~ .x[1, c("FlyID", "sex", "age", "treatment")])

grouped_plv <- fly_metadata |>
  group_by(sex, age, treatment) |>
  group_split() |>
  set_names(map_chr(., ~ paste(unique(.$sex), unique(.$age), unique(.$treatment), sep = "_"))) |>
  map(function(meta_group) {
    flies_in_group <- fly_data_response[meta_group$FlyID]
    compute_group_plv_matrix(flies_in_group)
  })

mean_plv_per_group <- map_dbl(grouped_plv, function(mat) {
  n <- nrow(mat)
  if (n < 2) return(NA_real_)
  mean(mat[lower.tri(mat)], na.rm = TRUE)
})

# Create a summary dataframe and split group names
group_sync_summary <- tibble(
  group = names(mean_plv_per_group),
  mean_PLV = mean_plv_per_group
) |>
  separate(group, into = c("sex", "age", "treatment"), sep = "_")

ggplot(group_sync_summary, aes(x = age, y = mean_PLV, fill = treatment)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~sex) +
  labs(
    x = "age",
    y = "mean PLV"
  ) +
  scale_fill_brewer(palette = "Set2") +
  ylim(0, 0.25) +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
