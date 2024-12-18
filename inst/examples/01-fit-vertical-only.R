library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(patchwork)
library(geowarp)

theme_set(theme_bw())

# Generate some simulated sounding profiles. These have a mean
# given by a quadratic function of depth, and the variance also
# increases with depth. They are spatially independent from each
# (just for  the convenience of the demonstration), and they have
# stationary vertical dependencies. Here, the horizontal coordinates
# are the same within a sounding, but the code does not assume that
simulated_depths <- seq(0.1, 19.9, by = 0.1)
simulated_depth_mean_beta <- c(0, 0.1, 0.02)
simulated_depth_variance_beta <- c(0, 0.1)
n_simulated_soundings <- 8
simulated_sounding_metadata <- tibble(
  name = sprintf('Sounding_%02d', seq_len(n_simulated_soundings)),
  x = runif(n_simulated_soundings, 0, 100),
  y = runif(n_simulated_soundings, 0, 100)
)
ar_rho <- 0.6
simulated_soundings <- simulated_sounding_metadata %>%
  group_by(name, x, y) %>%
  group_modify(~ {
    tibble(
      z = simulated_depths,
      target = (
        simulated_depth_mean_beta[1]
        + simulated_depth_mean_beta[2] * z
        + simulated_depth_mean_beta[3] * z ^ 2
        + sqrt(exp(
          simulated_depth_variance_beta[1]
          + simulated_depth_variance_beta[2] * z
        )) * as.vector(arima.sim(
          model = list(order = c(1, 0, 0), ar = ar_rho),
          n = length(simulated_depths),
          sd = sqrt(1 - ar_rho ^ 2)
        ))
      )
    )
  }, .groups = 'drop')

# Plot the simulated sounding profiles
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
  deviation_model = geowarp_vertical_only_deviation_model(
    axial_warping_unit = geowarp_bernstein_awu(
      # This allows for non-stationarity in the vertical correlation length.
      # The order roughly corresponds to resolution, so order 10 means we can
      # resolve features up to 20/10 = 2 depth units. This can be changed but
      # if it's increased too much there is a trade off between this and the
      # resolution of the variance process
      order = 10
    ),
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
  # that the optimisation completes in a reasonable time
  n_parents = 16,
  # This determines the number of times the optimisation routine is run from
  # different starting points; the best result is returned
  best_of = 1,
  # This turns on detailed logging - 0 gives no logging, otherwise it determines
  # how frequently messages are printed (1 = every iteration, 2 = every second
  # iteration, etc)
  trace = 10
)

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

# Plot the prediction bounds over the data
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
