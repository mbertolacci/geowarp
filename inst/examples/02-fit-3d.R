library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(patchwork)
library(geowarp)

theme_set(theme_bw())

# Simulate a 3-D grid with dependencies in depth and space, and
# which have mean and variance that vary with depth.
simulated_depth_mean_beta <- c(0, 0.1, 0.02)
simulated_depth_variance_beta <- c(0, 0.1)
horizontal_length_scale <- 40
vertical_length_scale <- 1

simulated_grid <- expand.grid(
  x = seq(0, 100, by = 10),
  y = seq(0, 100, by = 10),
  z = seq(0, 20, by = 1)
)
coords_grid <- cbind(
  simulated_grid$x / horizontal_length_scale,
  simulated_grid$y / horizontal_length_scale,
  simulated_grid$z / vertical_length_scale
)
d_grid <- fields::rdist(coords_grid)
Sigma_grid <- exp(-d_grid)
epsilon_grid <- as.vector(crossprod(chol(Sigma_grid), rnorm(nrow(Sigma_grid))))

simulated_grid$target <- (
  simulated_depth_mean_beta[1]
  + simulated_depth_mean_beta[2] * simulated_grid$z
  + simulated_depth_mean_beta[3] * simulated_grid$z ^ 2
  + sqrt(exp(
    simulated_depth_variance_beta[1]
    + simulated_depth_variance_beta[2] * simulated_grid$z
  )) * epsilon_grid
)

# Generate some sounding locations and depths
n_simulated_soundings <- 8
simulated_soundings <- tibble(
  name = sprintf('Sounding_%02d', seq_len(n_simulated_soundings)),
  x = runif(n_simulated_soundings, 0, 100),
  y = runif(n_simulated_soundings, 0, 100)
) %>%
  group_by(name, x, y) %>%
  group_modify(~ {
    tibble(
      z = seq(0.1, 19.9, by = 0.1)
    )
  }) %>%
  ungroup()

# Generate sounding data consistent with the simulated grid
coords_data <- cbind(
  simulated_soundings$x / horizontal_length_scale,
  simulated_soundings$y / horizontal_length_scale,
  simulated_soundings$z / vertical_length_scale
)
Sigma_data <- exp(-fields::rdist(coords_data))
Sigma_grid_to_data <- exp(-fields::rdist(coords_grid, coords_data))
epsilon_data_mean <- crossprod(Sigma_grid_to_data, solve(
  Sigma_grid,
  epsilon_grid
))
Sigma_data_conditional <- Sigma_data - crossprod(
  Sigma_grid_to_data,
  solve(Sigma_grid, Sigma_grid_to_data)
)
epsilon_data <- as.vector(
  epsilon_data_mean
  + crossprod(
    chol(Sigma_data_conditional),
    rnorm(nrow(Sigma_data_conditional))
  )
)
simulated_soundings <- simulated_soundings %>%
  mutate(
    target = (
      simulated_depth_mean_beta[1]
      + simulated_depth_mean_beta[2] * z
      + simulated_depth_mean_beta[3] * z ^ 2
      + sqrt(exp(
        simulated_depth_variance_beta[1]
        + simulated_depth_variance_beta[2] * z
      )) * epsilon_data
      # Add measurement error
      + rnorm(n(), 0, sd = 0.1)
    )
  )

# Plot the true grid and the sounding locations
ggplot() +
  geom_raster(
    data = simulated_grid,
    mapping = aes(x, y, fill = target)
  ) +
  geom_point(
    data = simulated_soundings %>%
      distinct(x, y),
    mapping = aes(x, y),
    colour = 'white'
  ) +
  scale_fill_viridis_c() +
  coord_fixed() +
  facet_wrap(~ z)

# Plot the soundings
ggplot() +
  geom_line(
    data = simulated_soundings,
    mapping = aes(z, target, colour = name)
  ) +
  scale_x_reverse() +
  coord_flip() +
  labs(x = 'Depth', y = 'Target')

# This code block specifies the model to fit.
model <- geowarp_model(
  # These quantities correspond to the target variable, and the spatial
  # coordinates
  variable = 'target',
  horizontal_coordinates = c('x', 'y'),
  horizontal_domains = list(c(0, 100), c(0, 100)),
  vertical_coordinate = 'z',
  vertical_domain = c(0, 20),
  # The following block configures the model for the mean-process, which
  # is shared between all the soundings.
  mean_model = geowarp_mean_model(
    # This configures the deterministic part of the mean model to be a linear
    # function of depth.
    fixed_formula = ~ z,
    # This adds a spline component to the mean model, allowing for
    # non-stationarity in the vertical correlation length.
    vertical_basis_functions = TRUE,
    # This is the spacing of the knots (in depth units) for the vertical basis
    # functions that determine the depth-varying mean.
    vertical_basis_function_delta = 2
  ),
  # The following block configures the model for the deviation process, which
  # handles deviation from the mean. This special version of the deviation
  # model is used to fit a vertical-only model in which horizontal dependencies
  # are ignored.
  deviation_model = geowarp_deviation_model(
    axial_warping_units = list(
      # These two allow for a single length scale for x and y
      geowarp_linear_awu(),
      geowarp_linear_awu(),
      geowarp_bernstein_awu(
        # This allows for non-stationarity in the vertical correlation length.
        # The order roughly corresponds to resolution, so order 10 means we can
        # resolve features up to 20/10 = 2 depth units. This can be changed but
        # if it's increased too much there is a trade off between this and the
        # resolution of the variance process
        order = 10
      )
    ),
    # This prevents the correlation structure from being rotated, which can
    # be hard to estimate
    geometric_warping_unit = NULL,
    variance_model = geowarp_variance_model(
      vertical_basis_functions = TRUE,
      # This is the spacing of the knots (in depth units) for the vertical basis
      # functions that determine the depth-varying variance.
      vertical_basis_function_delta = 2
    )
  )
)

fit <- geowarp_optimise(
  simulated_soundings,
  model,
  # This determines the number of parents used in the Vecchia approximation;
  # it needs to be high enough that the approximation is good, but low enough
  # that the optimisation completes in a reasonable time. It should be at least
  # twice as large as the number of soundings.
  n_parents = 32,
  # This determines the number of times the optimisation routine is run from
  # different starting points; the best result is returned
  best_of = 1,
  # This turns on detailed logging - 0 gives no logging, otherwise it determines
  # how frequently messages are printed (1 = every iteration, 2 = every second
  # iteration, etc)
  trace = 5,
  # This sets the number of threads to use for parallelisation
  threads = 4,
  # This sets the grain size for parallelisation. This should be neither too
  # high nor too low, and 16 is a good default.
  grain_size = 16
)

# These are the estimated horizontal (x, y) length scales
print(1 / fit$parameters$gamma_deviation_horizontal)

# Calculate some summary statistics over the depth profile
profile_df <- data.frame(
  x = 0,
  y = 0,
  z = seq(0, 20, by = 0.1)
) %>%
  mutate(
    estimated_mean = mean_profile(fit, df = .),
    estimated_stdev = sqrt(marginal_variance_profile(fit, df = .)),
    true_mean = (
      simulated_depth_mean_beta[1]
      + simulated_depth_mean_beta[2] * z
      + simulated_depth_mean_beta[3] * z ^ 2
    ),
    true_stdev = sqrt(exp(
      simulated_depth_variance_beta[1]
      + simulated_depth_variance_beta[2] * z
    ))
  )

# Plot the true and estimated mean profile, and the true and estimated standard
# deviation profile
wrap_plots(
  ggplot(profile_df) +
    geom_line(aes(z, true_mean), col = 'red') +
    geom_line(aes(z, estimated_mean), col = 'blue') +
    scale_x_reverse() +
    coord_flip() +
    labs(x = 'Depth', y = 'Mean', title = 'Mean profile'),
  ggplot(profile_df) +
    geom_line(aes(z, true_stdev), col = 'red') +
    geom_line(aes(z, estimated_stdev), col = 'blue') +
    scale_x_reverse() +
    coord_flip() +
    labs(x = 'Depth', y = 'Standard deviation', title = 'Standard deviation profile'),
  nrow = 1
)

# Plot the marginal mean and variance bounds over the data (these are not
# predictions per se, but a measure of the spread of the data according to the
# estimated model)
ggplot() +
  geom_line(
    data = simulated_soundings,
    mapping = aes(z, target, group = name),
    alpha = 0.1
  ) +
  geom_line(
    data = profile_df,
    mapping = aes(z, estimated_mean),
    col = 'blue'
  ) +
  geom_line(
    data = profile_df,
    mapping = aes(z, estimated_mean - 2 * estimated_stdev),
    col = 'blue',
    lty = 2
  ) +
  geom_line(
    data = profile_df,
    mapping = aes(z, estimated_mean + 2 * estimated_stdev),
    col = 'blue',
    lty = 2
  ) +
  scale_x_reverse() +
  coord_flip() +
  labs(x = 'Depth', y = 'Target', title = 'Target profile')

# The model allows for non-stationarity of the length scale with depth,
# which is determined by a "warping function". In practice, this should
# be linear here
warping_df <- data.frame(
  x = 0,
  y = 0,
  z = seq(0, 20, by = 0.1)
) %>%
  mutate(
    z_warped = warped_coordinates(fit, df = .)[, 3]
  )

# This should look roughly linear, since the data are actually stationary.
# The way this works is that the data are stationary with a length scale of
# one in the warped space
ggplot(warping_df, aes(z, z_warped)) +
  geom_line() +
  scale_x_reverse() +
  coord_flip() +
  labs(x = 'Depth', y = 'Warped depth', title = 'Warping function')

# Perform predictions from the soundings back to the full grid
prediction_df <- simulated_grid %>%
  select(x, y, z)

depth_delta <- 1
easting_delta <- 10
northing_delta <- 10
scaling <- c(depth_delta / easting_delta, depth_delta / northing_delta, 1)

parent_structure <- vecchia_parent_structure(
  fit$observed_df,
  fit$model,
  n_parents = 64,
  prediction_df = prediction_df,
  prediction_type = 'regular'
)

predictions <- predict(
  fit,
  prediction_df,
  include_mean = TRUE,
  include_samples = TRUE,
  n_samples = 256L,
  nugget = FALSE,
  parent_structure = parent_structure
)
prediction_df$target_mean <- predictions$mean
prediction_df$target_sd <- matrixStats::rowSds(predictions$samples)

# Plot the true and estimated 3-D grid values
bind_rows(
  simulated_grid %>%
    mutate(type = 'true', value = target),
  prediction_df %>%
    mutate(type = 'estimated', value = target_mean)
) %>%
  ggplot() +
    geom_raster(aes(x, y, fill = value)) +
    scale_fill_viridis_c() +
    coord_fixed() +
    facet_grid(z ~ type)

# Plot the prediction stdev
ggplot() +
  geom_raster(
    data = prediction_df,
    mapping = aes(x, y, fill = target_sd)
  ) +
  geom_point(
    data = simulated_soundings %>% distinct(x, y),
    mapping = aes(x, y),
    colour = 'white'
  ) +
  scale_fill_viridis_c() +
  coord_fixed() +
  facet_wrap(~ z)
