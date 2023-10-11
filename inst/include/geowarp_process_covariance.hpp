#ifndef GEOWARP_PROCESS_COVARIANCE_HPP
#define GEOWARP_PROCESS_COVARIANCE_HPP

template <typename T_ss, typename T_sd, typename T_x>
class geowarp_process_covariance_vari : public stan::math::vari {
public:
  const size_t size_;
  const size_t D_;
  const size_t size_ltri_;

  const double smoothness_;

  double *tmp_;

  stan::math::vari *sigma_squared_nugget_vari_;
  stan::math::vari **deviation_sd_vari_;
  stan::math::vari **x_vari_;

  stan::math::vari **cov_lower_;
  stan::math::vari **cov_diag_;

  geowarp_process_covariance_vari(
    const T_ss& sigma_squared_nugget,
    const T_sd& deviation_sd,
    const std::vector<Eigen::Matrix<T_x, Eigen::Dynamic, 1> >& x,
    const double& smoothness
  ) : vari(0.0),
      size_(x.size()),
      D_(x[0].rows()),
      size_ltri_(size_ * (size_ - 1) / 2),
      smoothness_(smoothness),
      tmp_(
        ChainableStack::instance_->memalloc_.alloc_array<double>(size_ltri_)
      ),
      sigma_squared_nugget_vari_(sigma_squared_nugget.vi_),
      deviation_sd_vari_(
        ChainableStack::instance_->memalloc_.alloc_array<stan::math::vari *>(size_)
      ),
      x_vari_(ChainableStack::instance_->memalloc_.alloc_array<stan::math::vari *>(
        D_ * size_
      )),
      cov_lower_(
        ChainableStack::instance_->memalloc_.alloc_array<stan::math::vari *>(size_ltri_)
      ),
      cov_diag_(
        ChainableStack::instance_->memalloc_.alloc_array<stan::math::vari *>(size_)
      ) {

    for (size_t j = 0; j < size_; ++j) {
      deviation_sd_vari_[j] = deviation_sd.coeff(j).vi_;
      for (size_t k = 0; k < D_; ++k) {
        x_vari_[j * D_ + k] = x[j].coeff(k).vi_;
      }
    }

    size_t pos = 0;
    if (smoothness_ == 0.5) {
      for (size_t j = 0; j < size_ - 1; ++j) {
        for (size_t i = j + 1; i < size_; ++i) {
          tmp_[pos] = distance(value_of(x[i]), value_of(x[j]));
          cov_lower_[pos] = new stan::math::vari(
            value_of(deviation_sd.coeff(j))
            * value_of(deviation_sd.coeff(i))
            * std::exp(-tmp_[pos]),
            false);
          ++pos;
        }
      }
    } else if (smoothness_ == 1.5) {
      for (size_t j = 0; j < size_ - 1; ++j) {
        for (size_t i = j + 1; i < size_; ++i) {
          double sqrt3_dist = std::sqrt(3) * distance(value_of(x[i]), value_of(x[j]));
          tmp_[pos] = (
            value_of(deviation_sd.coeff(j))
            * value_of(deviation_sd.coeff(i))
            * std::exp(-sqrt3_dist)
          );
          cov_lower_[pos] = new stan::math::vari(
            (1 + sqrt3_dist) * tmp_[pos],
            false
          );
          ++pos;
        }
      }
    }

    for (size_t i = 0; i < size_; ++i) {
      cov_diag_[i] = new stan::math::vari(
        value_of(deviation_sd.coeff(i)) * value_of(deviation_sd.coeff(i))
        + value_of(sigma_squared_nugget),
        false
      );
    }
  }

  virtual void chain() {
    size_t pos = 0;
    if (smoothness_ == 0.5) {
      for (size_t j = 0; j < size_ - 1; ++j) {
        for (size_t i = j + 1; i < size_; ++i) {
          vari *el_low = cov_lower_[pos];
          double base = el_low->adj_ * el_low->val_;

          deviation_sd_vari_[j]->adj_ += base / deviation_sd_vari_[j]->val_;
          deviation_sd_vari_[i]->adj_ += base / deviation_sd_vari_[i]->val_;

          for (size_t k = 0; k < D_; ++k) {
            double diff = x_vari_[j * D_ + k]->val_ - x_vari_[i * D_ + k]->val_;
            double left = -diff * base / tmp_[pos];
            x_vari_[j * D_ + k]->adj_ += left;
            x_vari_[i * D_ + k]->adj_ -= left;
          }

          ++pos;
        }
      }
    } else if (smoothness_ == 1.5) {
      for (size_t j = 0; j < size_ - 1; ++j) {
        for (size_t i = j + 1; i < size_; ++i) {
          vari *el_low = cov_lower_[pos];
          double base = el_low->adj_ * el_low->val_;

          deviation_sd_vari_[j]->adj_ += base / deviation_sd_vari_[j]->val_;
          deviation_sd_vari_[i]->adj_ += base / deviation_sd_vari_[i]->val_;

          for (size_t k = 0; k < D_; ++k) {
            double diff = x_vari_[j * D_ + k]->val_ - x_vari_[i * D_ + k]->val_;
            double left = -3 * el_low->adj_ * diff * tmp_[pos];
            x_vari_[j * D_ + k]->adj_ += left;
            x_vari_[i * D_ + k]->adj_ -= left;
          }

          ++pos;
        }
      }
    }

    for (size_t i = 0; i < size_; ++i) {
      vari *el = cov_diag_[i];
      sigma_squared_nugget_vari_->adj_ += el->adj_;
      deviation_sd_vari_[i]->adj_ += 2 * el->adj_ * deviation_sd_vari_[i]->val_;
    }
  }
};

template <typename T_ss, typename T_sd, typename T_x>
Eigen::Matrix<
  stan::promote_args_t<T_ss, stan::value_type_t<T_sd>, T_x>,
  Eigen::Dynamic,
  Eigen::Dynamic
>
geowarp_process_covariance(
  const T_ss& sigma_squared_nugget,
  const T_sd& deviation_sd,
  const std::vector<Eigen::Matrix<T_x, Eigen::Dynamic, 1> >& x,
  const double& smoothness,
  std::ostream* pstream__
) {
  int N = deviation_sd.size();
  Eigen::Matrix<
    stan::promote_args_t<T_ss, stan::value_type_t<T_sd>, T_x>,
    Eigen::Dynamic,
    Eigen::Dynamic
  > K_block(N, N);
  if (smoothness == 0.5) {
    for (int j = 0; j < N; ++j) {
      K_block.coeffRef(j, j) = (
        deviation_sd[j] * deviation_sd[j] + sigma_squared_nugget
      );
      for (int k = j + 1; k < N; ++k) {
        auto value_i_j_k = (
          deviation_sd[j]
          * deviation_sd[k]
          * exp(-distance(x[j], x[k]))
        );
        K_block.coeffRef(j, k) = value_i_j_k;
        K_block.coeffRef(k, j) = value_i_j_k;
      }
    }
  } else if (smoothness == 1.5) {
    for (int j = 0; j < N; ++j) {
      K_block.coeffRef(j, j) = (
        deviation_sd[j] * deviation_sd[j] + sigma_squared_nugget
      );
      for (int k = j + 1; k < N; ++k) {
        auto s3_distance_i_j_k = std::sqrt(3) * distance(x[j], x[k]);
        auto value_i_j_k = (
          deviation_sd[j]
          * deviation_sd[k]
          * (1 + s3_distance_i_j_k)
          * exp(-s3_distance_i_j_k)
        );
        K_block.coeffRef(j, k) = value_i_j_k;
        K_block.coeffRef(k, j) = value_i_j_k;
      }
    }
  }
  return K_block;
}

template<>
inline
Eigen::Matrix<
  stan::math::var,
  Eigen::Dynamic,
  Eigen::Dynamic
>
geowarp_process_covariance(
  const stan::math::var& sigma_squared_nugget,
  const Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>& deviation_sd,
  const std::vector<Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> >& x,
  const double& smoothness,
  std::ostream* pstream__
) {
  using stan::math::var;
  size_t x_size = x.size();

  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> output(x_size, x_size);
  if (x_size == 0) {
    return output;
  }

  geowarp_process_covariance_vari<var, Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>, var> *baseVari
    = new geowarp_process_covariance_vari<var, Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>, var>(
      sigma_squared_nugget,
      deviation_sd,
      x,
      smoothness
    );

  size_t pos = 0;
  for (size_t j = 0; j < x_size - 1; ++j) {
    for (size_t i = (j + 1); i < x_size; ++i) {
      output.coeffRef(i, j).vi_ = baseVari->cov_lower_[pos];
      output.coeffRef(j, i).vi_ = output.coeffRef(i, j).vi_;
      ++pos;
    }
    output.coeffRef(j, j).vi_ = baseVari->cov_diag_[j];
  }
  output.coeffRef(x_size - 1, x_size - 1).vi_ = baseVari->cov_diag_[x_size - 1];
  return output;
}

// 1-D version

template <typename T_ss, typename T_sd, typename T_x>
class geowarp_process_covariance_1d_vari : public stan::math::vari {
public:
  const size_t size_;
  const size_t size_ltri_;

  const double smoothness_;

  double *tmp_;

  stan::math::vari *sigma_squared_nugget_vari_;
  stan::math::vari **deviation_sd_vari_;
  stan::math::vari **x_vari_;

  stan::math::vari **cov_lower_;
  stan::math::vari **cov_diag_;

  geowarp_process_covariance_1d_vari(
    const T_ss& sigma_squared_nugget,
    const T_sd& deviation_sd,
    const std::vector<T_x>& x,
    const double& smoothness
  ) : vari(0.0),
      size_(x.size()),
      size_ltri_(size_ * (size_ - 1) / 2),
      smoothness_(stan::math::value_of(smoothness)),
      tmp_(
        ChainableStack::instance_->memalloc_.alloc_array<double>(size_ltri_)
      ),
      sigma_squared_nugget_vari_(sigma_squared_nugget.vi_),
      deviation_sd_vari_(
        ChainableStack::instance_->memalloc_.alloc_array<stan::math::vari *>(size_)
      ),
      x_vari_(ChainableStack::instance_->memalloc_.alloc_array<stan::math::vari *>(
        size_
      )),
      cov_lower_(
        ChainableStack::instance_->memalloc_.alloc_array<stan::math::vari *>(size_ltri_)
      ),
      cov_diag_(
        ChainableStack::instance_->memalloc_.alloc_array<stan::math::vari *>(size_)
      ) {

    for (size_t j = 0; j < size_; ++j) {
      deviation_sd_vari_[j] = deviation_sd.coeff(j).vi_;
      x_vari_[j] = x[j].vi_;
    }

    size_t pos = 0;
    if (smoothness_ == 0.5) {
      for (size_t j = 0; j < size_ - 1; ++j) {
        for (size_t i = j + 1; i < size_; ++i) {
          tmp_[pos] = abs(value_of(x[i]) - value_of(x[j]));
          cov_lower_[pos] = new stan::math::vari(
            value_of(deviation_sd.coeff(j))
            * value_of(deviation_sd.coeff(i))
            * std::exp(-tmp_[pos]),
            false);
          ++pos;
        }
      }
    } else if (smoothness_ == 1.5) {
      for (size_t j = 0; j < size_ - 1; ++j) {
        for (size_t i = j + 1; i < size_; ++i) {
          double sqrt3_dist = std::sqrt(3) * abs(value_of(x[i]) - value_of(x[j]));
          tmp_[pos] = (
            value_of(deviation_sd.coeff(j))
            * value_of(deviation_sd.coeff(i))
            * std::exp(-sqrt3_dist)
          );
          cov_lower_[pos] = new stan::math::vari(
            (1 + sqrt3_dist) * tmp_[pos],
            false
          );
          ++pos;
        }
      }
    }

    for (size_t i = 0; i < size_; ++i) {
      cov_diag_[i] = new stan::math::vari(
        value_of(deviation_sd.coeff(i)) * value_of(deviation_sd.coeff(i))
        + value_of(sigma_squared_nugget),
        false
      );
    }
  }

  virtual void chain() {
    size_t pos = 0;
    if (smoothness_ == 0.5) {
      for (size_t j = 0; j < size_ - 1; ++j) {
        for (size_t i = j + 1; i < size_; ++i) {
          vari *el_low = cov_lower_[pos];
          double base = el_low->adj_ * el_low->val_;

          deviation_sd_vari_[j]->adj_ += base / deviation_sd_vari_[j]->val_;
          deviation_sd_vari_[i]->adj_ += base / deviation_sd_vari_[i]->val_;

          double diff = x_vari_[j]->val_ - x_vari_[i]->val_;
          double left = -diff * base / tmp_[pos];
          x_vari_[j]->adj_ += left;
          x_vari_[i]->adj_ -= left;

          ++pos;
        }
      }
    } else if (smoothness_ == 1.5) {
      for (size_t j = 0; j < size_ - 1; ++j) {
        for (size_t i = j + 1; i < size_; ++i) {
          vari *el_low = cov_lower_[pos];
          double base = el_low->adj_ * el_low->val_;

          deviation_sd_vari_[j]->adj_ += base / deviation_sd_vari_[j]->val_;
          deviation_sd_vari_[i]->adj_ += base / deviation_sd_vari_[i]->val_;

          double diff = x_vari_[j]->val_ - x_vari_[i]->val_;
          double left = -3 * el_low->adj_ * diff * tmp_[pos];
          x_vari_[j]->adj_ += left;
          x_vari_[i]->adj_ -= left;

          ++pos;
        }
      }
    }

    for (size_t i = 0; i < size_; ++i) {
      vari *el = cov_diag_[i];
      sigma_squared_nugget_vari_->adj_ += el->adj_;
      deviation_sd_vari_[i]->adj_ += 2 * el->adj_ * deviation_sd_vari_[i]->val_;
    }
  }
};

template <typename T_ss, typename T_sd, typename T_x>
Eigen::Matrix<
  stan::promote_args_t<T_ss, stan::value_type_t<T_sd>, T_x>,
  Eigen::Dynamic,
  Eigen::Dynamic
>
inline geowarp_process_covariance_1d(
  const T_ss& sigma_squared_nugget,
  const T_sd& deviation_sd,
  const std::vector<T_x>& x,
  const double& smoothness,
  std::ostream* pstream__
) {
  int N = deviation_sd.size();
  Eigen::Matrix<
  stan::promote_args_t<T_ss, stan::value_type_t<T_sd>, T_x>,
    Eigen::Dynamic,
    Eigen::Dynamic
  > K_block(N, N);
  if (smoothness == 0.5) {
    for (int j = 0; j < N; ++j) {
      K_block.coeffRef(j, j) = (
        pow(deviation_sd[j], 2) + sigma_squared_nugget
      );
      for (int k = j + 1; k < N; ++k) {
        auto value_i_j_k = (
          deviation_sd[j]
          * deviation_sd[k]
          * exp(-abs(x[j] - x[k]))
        );
        K_block.coeffRef(j, k) = value_i_j_k;
        K_block.coeffRef(k, j) = value_i_j_k;
      }
    }
  } else if (smoothness == 1.5) {
    for (int j = 0; j < N; ++j) {
      K_block.coeffRef(j, j) = (
        pow(deviation_sd[j], 2) + sigma_squared_nugget
      );
      for (int k = j + 1; k < N; ++k) {
        auto s3_distance_i_j_k = std::sqrt(3) * abs(x[j] - x[k]);
        auto value_i_j_k = (
          deviation_sd[j]
          * deviation_sd[k]
          * (1 + s3_distance_i_j_k)
          * exp(-s3_distance_i_j_k)
        );
        K_block.coeffRef(j, k) = value_i_j_k;
        K_block.coeffRef(k, j) = value_i_j_k;
      }
    }
  }
  return K_block;
}

template<>
Eigen::Matrix<
  stan::math::var,
  Eigen::Dynamic,
  Eigen::Dynamic
>
inline geowarp_process_covariance_1d(
  const stan::math::var& sigma_squared_nugget,
  const Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>& deviation_sd,
  const std::vector<stan::math::var>& x,
  const double& smoothness,
  std::ostream* pstream__
) {
  using stan::math::var;
  size_t x_size = x.size();

  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> output(x_size, x_size);
  if (x_size == 0) {
    return output;
  }

  geowarp_process_covariance_1d_vari<var, Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>, var> *baseVari
    = new geowarp_process_covariance_1d_vari<var, Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>, var>(
      sigma_squared_nugget,
      deviation_sd,
      x,
      smoothness
    );

  size_t pos = 0;
  for (size_t j = 0; j < x_size - 1; ++j) {
    for (size_t i = (j + 1); i < x_size; ++i) {
      output.coeffRef(i, j).vi_ = baseVari->cov_lower_[pos];
      output.coeffRef(j, i).vi_ = output.coeffRef(i, j).vi_;
      ++pos;
    }
    output.coeffRef(j, j).vi_ = baseVari->cov_diag_[j];
  }
  output.coeffRef(x_size - 1, x_size - 1).vi_ = baseVari->cov_diag_[x_size - 1];
  return output;
}

#endif  // GEOWARP_PROCESS_COVARIANCE_HPP
