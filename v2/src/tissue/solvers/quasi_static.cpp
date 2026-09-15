//
// QuasiStatic: growth stepping with mechanical equilibrium solved, not
// integrated.
//
// Why this exists. The framework's mechanics are overdamped, dx/dt = F(x),
// so an explicit integrator's step size is bounded by the *stiffest* elastic
// mode: h < ~2/K. On the hook shell that is h ~ 2e-5 h, while the process
// being modelled (differential growth unrolling the organ) takes hours -
// a stiffness ratio near 1e5. Nearly all of the work resolves elastic
// transients that the trajectory does not depend on, because the tissue is
// always essentially at force balance.
//
// It is also why the model needs Y and P scaled up together: the *softest*
// mode (whole-arm rotation) has a tiny restoring force but full per-vertex
// drag, so at biological stiffness the arm barely moves over the growth
// timescale. Scaling Y and P leaves the equilibrium shape untouched and
// speeds that mode up - at the cost of making the stiff modes stiffer still.
//
// Solving equilibrium instead of integrating toward it removes both problems
// at once. The equilibrium of F(x) = 0 does not depend on the drag or on the
// overall scale of Y and P (only on their ratios), so the mobility fudge is
// unnecessary, and the growth step is limited by growth accuracy alone.
//
// Relaxation uses FIRE (Bitzek et al., PRL 97:170201, 2006): damped inertial
// dynamics with an adaptive timestep. Inertial dynamics is the point - its
// stable step scales as 1/sqrt(K) rather than 1/K, and it converges in
// O(sqrt(kappa)) iterations rather than O(kappa). It needs nothing but the
// force already computed by derivs(), so every reaction works unchanged.
//
// Solver file:
//   QuasiStatic
//   <startTime> <endTime>
//   <printFlag> <numPrint>
//   <h_growth> <force_tol> <max_relax_iterations>
//
// Comparing against another solver: QuasiStatic relaxes to force balance
// *before* its first print, so its frame 0 is the relaxed configuration while
// the explicit solvers print the initial file as given. Absolute values at
// t = 0 therefore differ by whatever the initial state was out of balance by;
// compare growth factors against each run's own frame 0, not absolute values
// across solvers.
//
// force_tol is relative: relaxation stops when the largest vertex force falls
// below force_tol times the largest force seen at the start of that step's
// relaxation (plus an absolute floor), so it adapts to the force scale of the
// model rather than needing units.
//
// IMPORTANT - this does not reproduce RK5Adaptive on an existing model, and
// is not meant to. Integrating dx/dt = F with per-vertex drag leaves the
// tissue lagging behind force balance whenever the drag timescale approaches
// the growth timescale, which on the hook shell it does: the same model that
// RK5Adaptive opens from 159 to 122 degrees over 2 h opens to 6 degrees here,
// and halving h_growth and tightening force_tol tenfold moves that by under
// 3 degrees, so the difference is the lag, not discretisation error. A model
// tuned against the lagged dynamics therefore has its growth rate absorbing
// a mobility artefact and must be recalibrated (roughly, k_growth divided by
// the ratio of opening rates) before its fit means the same thing here.
//
// Note also that uniformly scaling Y, Y_fiber and P buys nothing here: FIRE
// converges in O(sqrt(condition number)) iterations and uniform scaling
// leaves the condition number unchanged (measured: 62 s scaled vs 74 s at
// biological stiffness for the same 2 h). What it does buy is that the
// scaling is no longer *needed*, since equilibrium does not depend on it.
//
#include <algorithm>
#include <cmath>
#include <iostream>

#include <vector>
#include "tissue/solvers/base_solver.h"
#include "tissue/solvers/solvers.h"

namespace tissue {

QuasiStatic::QuasiStatic(Tissue *T, std::istream &in) : BaseSolver(T) {
  in >> startTime_ >> endTime_;
  in >> printFlag_ >> numPrint_;
  in >> hGrowth_ >> forceTol_ >> maxRelax_;
  t_ = startTime_;
  if (hGrowth_ <= 0.0) {
    std::cerr << "QuasiStatic: growth step must be positive." << std::endl;
    std::exit(-1);
  }
  if (forceTol_ <= 0.0)
    forceTol_ = 1e-3;
  if (maxRelax_ <= 0)
    maxRelax_ = 2000;
}

// Largest |force| over every mechanical degree of freedom. NaN propagates
// deliberately: std::max would swallow it and make a diverged state look
// converged.
double QuasiStatic::maxForce() const {
  double m = 0.0;
  auto f = vertexDerivs_.flat();
  for (size_t k = 0; k < f.size(); ++k) {
    const double a = std::fabs(f[k]);
    if (std::isnan(a))
      return a;
    m = std::max(m, a);
  }
  auto cf = cellDerivs_.flat();
  for (size_t idx : posIndex_) {
    const double a = std::fabs(cf[idx]);
    if (std::isnan(a))
      return a;
    m = std::max(m, a);
  }
  return m;
}

// FIRE relaxation of the vertex positions to force balance. Returns the
// number of force evaluations used; velocity_ is kept between calls so a
// growth step that barely perturbs the equilibrium restarts almost converged.
size_t QuasiStatic::relax() {
  constexpr size_t kNMin = 5;
  constexpr double kFInc = 1.1, kFDec = 0.5;
  constexpr double kAlpha0 = 0.1, kFAlpha = 0.99;

  if (!velocity_.sameShape(vertexData_))
    velocity_.reshapeLike(vertexData_);
  auto v = velocity_.flat();
  std::fill(v.begin(), v.end(), 0.0);
  posVel_.assign(posIndex_.size(), 0.0);

  forceOnly();
  size_t evals = 1;
  const double f0 = maxForce();
  if (f0 <= 0.0)
    return evals;
  // Absolute floor keeps a nearly force-free state from chasing round-off.
  const double target = std::max(forceTol_ * f0, 1e-10 * scaleHint_);

  double dt = dt0_;
  double alpha = kAlpha0;
  size_t nPos = 0;

  for (size_t it = 0; it < static_cast<size_t>(maxRelax_); ++it) {
    auto f = vertexDerivs_.flat();
    auto x = vertexData_.flat();
    v = velocity_.flat();

    auto cf = cellDerivs_.flat();
    auto cx = cellData_.flat();

    double power = 0.0, vNorm = 0.0, fNorm = 0.0;
    for (size_t k = 0; k < f.size(); ++k) {
      power += f[k] * v[k];
      vNorm += v[k] * v[k];
      fNorm += f[k] * f[k];
    }
    for (size_t j = 0; j < posIndex_.size(); ++j) {
      const double fj = cf[posIndex_[j]];
      power += fj * posVel_[j];
      vNorm += posVel_[j] * posVel_[j];
      fNorm += fj * fj;
    }
    vNorm = std::sqrt(vNorm);
    fNorm = std::sqrt(fNorm);

    if (power > 0.0) {
      if (++nPos > kNMin) {
        dt = std::min(dt * kFInc, dtMax_);
        alpha *= kFAlpha;
      }
    } else {
      // Moving uphill: discard the velocity and shorten the step.
      nPos = 0;
      dt *= kFDec;
      alpha = kAlpha0;
      std::fill(v.begin(), v.end(), 0.0);
      std::fill(posVel_.begin(), posVel_.end(), 0.0);
      vNorm = 0.0;
    }

    // Semi-implicit Euler MD step with the FIRE velocity mixing, over the
    // vertices and the center-triangulation vertices together.
    for (size_t k = 0; k < f.size(); ++k)
      v[k] += dt * f[k];
    for (size_t j = 0; j < posIndex_.size(); ++j)
      posVel_[j] += dt * cf[posIndex_[j]];
    if (fNorm > 0.0) {
      double vn = 0.0;
      for (size_t k = 0; k < v.size(); ++k)
        vn += v[k] * v[k];
      for (size_t j = 0; j < posVel_.size(); ++j)
        vn += posVel_[j] * posVel_[j];
      vn = std::sqrt(vn);
      const double mix = alpha * vn / fNorm;
      for (size_t k = 0; k < v.size(); ++k)
        v[k] = (1.0 - alpha) * v[k] + mix * f[k];
      for (size_t j = 0; j < posVel_.size(); ++j)
        posVel_[j] = (1.0 - alpha) * posVel_[j] + mix * cf[posIndex_[j]];
    }
    for (size_t k = 0; k < x.size(); ++k)
      x[k] += dt * v[k];
    for (size_t j = 0; j < posIndex_.size(); ++j)
      cx[posIndex_[j]] += dt * posVel_[j];

    forceOnly();
    ++evals;
    if (maxForce() <= target)
      return evals;
  }
  ++relaxNotConverged_;
  return evals;
}

// Forces only: prescribed-velocity reactions are held out, because a term
// that does not vanish at equilibrium cannot be relaxed - FIRE would drive it
// without bound instead of converging.
void QuasiStatic::forceOnly() {
  if (hasPrescribed_)
    T_->derivsSplit(cellData_, wallData_, vertexData_, cellDerivs_, wallDerivs_,
                    vertexDerivs_, vertexVel_);
  else
    T_->derivs(cellData_, wallData_, vertexData_, cellDerivs_, wallDerivs_,
               vertexDerivs_);
}

// Imposed motion over one growth step, applied before relaxation.
void QuasiStatic::applyPrescribed(double h) {
  if (!hasPrescribed_)
    return;
  T_->derivsSplit(cellData_, wallData_, vertexData_, cellDerivs_, wallDerivs_,
                  vertexDerivs_, vertexVel_);
  auto x = vertexData_.flat();
  auto v = vertexVel_.flat();
  for (size_t k = 0; k < x.size(); ++k)
    x[k] += h * v[k];
}

void QuasiStatic::simulate() {
  if (endTime_ - startTime_ <= 0.0) {
    std::cerr << "QuasiStatic::simulate() Wrong time borders. No simulation "
                 "performed." << std::endl;
    std::exit(-1);
  }
  T_->initiateReactions(cellData_, wallData_, vertexData_, cellDerivs_,
                        wallDerivs_, vertexDerivs_);
  cellDerivs_.reshapeLike(cellData_);
  wallDerivs_.reshapeLike(wallData_);
  vertexDerivs_.reshapeLike(vertexData_);
  T_->initiateDirection(cellData_, wallData_, vertexData_, cellDerivs_,
                        wallDerivs_, vertexDerivs_);
  initPrintSchedule(1e-12);
  hasPrescribed_ = T_->hasPrescribedVelocity();
  vertexVel_.reshapeLike(vertexData_);
  if (hasPrescribed_)
    std::cerr << "QuasiStatic: model prescribes vertex velocity; integrating "
                 "that separately from the relaxation." << std::endl;

  // Set the FIRE timestep from the fastest elastic mode, estimated once from
  // the initial force response: dt_stable ~ 2/sqrt(K) for inertial dynamics.
  // A single probe displacement gives K without needing the Jacobian.
  calibrateStep();

  t_ = startTime_;
  numOk_ = numBad_ = 0;
  size_t totalEvals = 0;

  // Start from equilibrium so the first growth step sees balanced forces.
  totalEvals += relax();

  for (;;) {
    if (doPrint_ && t_ >= printTime_) {
      printTime_ += printDeltaTime_;
      print();
    }
    double h = hGrowth_;
    double tMin = endTime_ < printTime_ ? endTime_ : printTime_;
    if (t_ + h > tMin)
      h = tMin - t_;
    if (h <= 0.0)
      h = hGrowth_;

    // --- growth / chemistry step, vertices held fixed ------------------
    // These variables are non-stiff: their rates change over hours, so a
    // Heun (predictor-corrector) step at the growth step size is ample.
    applyPrescribed(h);
    forceOnly();
    ++totalEvals;
    cellStart_.copyFrom(cellData_);
    wallStart_.copyFrom(wallData_);
    cellK1_.copyFrom(cellDerivs_);
    wallK1_.copyFrom(wallDerivs_);
    addScaled(cellData_, cellK1_, h);
    addScaled(wallData_, wallK1_, h);
    restorePositional(cellData_, cellStart_);

    forceOnly();
    ++totalEvals;
    // y = y0 + h/2 (k1 + k2)
    combineHeun(cellData_, cellStart_, cellK1_, cellDerivs_, h);
    combineHeun(wallData_, wallStart_, wallK1_, wallDerivs_, h);
    // Positions are relaxed below, never integrated as growth.
    restorePositional(cellData_, cellStart_);

    // --- re-establish mechanical equilibrium ---------------------------
    totalEvals += relax();

    t_ += h;
    ++numOk_;
    postStep(h);
    if (!velocity_.sameShape(vertexData_))
      velocity_.reshapeLike(vertexData_);

    if (t_ >= endTime_ - 1e-12) {
      if (doPrint_) {
        forceOnly();
        print();
      }
      std::cerr << "Simulation done. growth steps: " << numOk_
                << ", force evaluations: " << totalEvals << " ("
                << (totalEvals / std::max(1u, numOk_)) << " per step)";
      if (relaxNotConverged_)
        std::cerr << "\n  WARNING: " << relaxNotConverged_ << " of " << numOk_
                  << " growth steps hit the relaxation cap (" << maxRelax_
                  << ") without reaching force balance - those steps are NOT"
                     " converged. Raise the cap or loosen force_tol.";
      std::cerr << std::endl;
      return;
    }
  }
}

// One probe displacement gives the stiffest response without a Jacobian:
// push every vertex along its current force direction by a small amount and
// measure how much the force changes. K ~ |dF|/|dx|, and FIRE is stable for
// dt < 2/sqrt(K).
void QuasiStatic::calibrateStep() {
  // Flat indices of the cell columns that are positions rather than
  // concentrations (center-triangulation vertices).
  posIndex_.clear();
  const std::vector<size_t> cols = T_->positionalCellVariables();
  for (size_t c = 0; c < cellData_.rows(); ++c) {
    const size_t base = cellData_.rowBegin(c);
    for (size_t col : cols)
      if (col < cellData_.rowSize(c))
        posIndex_.push_back(base + col);
  }
  std::cerr << "QuasiStatic: " << cols.size()
            << " positional cell variables per cell (" << posIndex_.size()
            << " mechanical DOFs in cell rows)" << std::endl;

  cellStart_.reshapeLike(cellData_);
  wallStart_.reshapeLike(wallData_);
  cellK1_.reshapeLike(cellData_);
  wallK1_.reshapeLike(wallData_);
  velocity_.reshapeLike(vertexData_);

  T_->derivs(cellData_, wallData_, vertexData_, cellDerivs_, wallDerivs_,
             vertexDerivs_);
  auto x = vertexData_.flat();

  double xScale = 0.0;
  for (size_t k = 0; k < x.size(); ++k)
    xScale = std::max(xScale, std::fabs(x[k]));
  scaleHint_ = std::max(xScale, 1.0);

  // Largest Hessian eigenvalue by matrix-free power iteration. FIRE is stable
  // for dt < 2/sqrt(lambda_max), so what matters is the *stiffest* mode - the
  // shortest edge, typically - and nothing else.
  //
  // A single probe along the force direction (what this used to do) measures
  // the stiffness the force happens to sample, which is dominated by whichever
  // modes carry large forces. Those are the driven, soft ones. On a fine mesh
  // that underestimates lambda_max badly, dt0 comes out too large, and FIRE
  // spends its budget thrashing: every overshoot trips the P < 0 branch, which
  // zeroes the velocity and halves dt, so it never builds up speed and hits
  // the iteration cap without reaching force balance.
  //
  // -H v is obtained from a finite difference of the force, since F = -grad U.
  const double eps = 1e-6 * scaleHint_;
  std::vector<double> v(x.size()), f0(x.size()), saved(x.begin(), x.end());
  {
    auto f = vertexDerivs_.flat();
    std::copy(f.begin(), f.end(), f0.begin());
  }
  // Seed off the force direction, falling back to a fixed pattern if the
  // state is already balanced (a deterministic seed keeps runs reproducible).
  double vNorm = 0.0;
  for (size_t k = 0; k < v.size(); ++k) {
    v[k] = f0[k];
    vNorm += v[k] * v[k];
  }
  if (std::sqrt(vNorm) <= 0.0) {
    for (size_t k = 0; k < v.size(); ++k)
      v[k] = (k % 2 == 0) ? 1.0 : -1.0;
    vNorm = static_cast<double>(v.size());
  }
  vNorm = std::sqrt(vNorm);
  for (size_t k = 0; k < v.size(); ++k)
    v[k] /= vNorm;

  double lambda = 1.0, prev = 0.0;
  size_t iters = 0;
  for (; iters < 24; ++iters) {
    for (size_t k = 0; k < x.size(); ++k)
      x[k] = saved[k] + eps * v[k];
    T_->derivs(cellData_, wallData_, vertexData_, cellDerivs_, wallDerivs_,
               vertexDerivs_);
    auto f1 = vertexDerivs_.flat();
    double n = 0.0;
    for (size_t k = 0; k < v.size(); ++k) {
      v[k] = (f0[k] - f1[k]) / eps; // -H v
      n += v[k] * v[k];
    }
    n = std::sqrt(n);
    if (!(n > 0.0))
      break;
    lambda = n; // |H v| with |v| = 1 -> Rayleigh-like estimate
    if (iters == 0)
      firstProbe_ = lambda; // what a single force-direction probe would give
    for (size_t k = 0; k < v.size(); ++k)
      v[k] /= n;
    if (iters > 2 && std::fabs(lambda - prev) <= 0.02 * lambda) {
      ++iters;
      break;
    }
    prev = lambda;
  }
  for (size_t k = 0; k < x.size(); ++k)
    x[k] = saved[k];
  const double K = std::max(lambda, 1e-12);

  // Stability limit for the inertial step, then a margin: power iteration
  // converges from below, so leave room for what it has not yet resolved.
  const double dtStable = 2.0 / std::sqrt(K);
  dt0_ = 0.10 * dtStable;
  dtMax_ = 0.50 * dtStable;
  std::cerr << "QuasiStatic: lambda_max ~ " << K << " (" << iters
            << " power iterations), FIRE dt0 = " << dt0_
            << ", dt_max = " << dtMax_ << std::endl;
  if (firstProbe_ > 0.0)
    std::cerr << "  (a single force-direction probe would have said "
              << firstProbe_ << ", i.e. " << (100.0 * (K / firstProbe_ - 1.0))
              << "% low - the gap grows with mesh fineness)" << std::endl;
}

// Put the positional cell columns back to their pre-growth-step values.
void QuasiStatic::restorePositional(Matrix &y, const Matrix &y0) const {
  auto yf = y.flat();
  auto y0f = y0.flat();
  for (size_t idx : posIndex_)
    yf[idx] = y0f[idx];
}

void QuasiStatic::addScaled(Matrix &y, const Matrix &k, double h) {
  auto yf = y.flat();
  auto kf = k.flat();
  for (size_t i = 0; i < yf.size(); ++i)
    yf[i] += h * kf[i];
}

void QuasiStatic::combineHeun(Matrix &y, const Matrix &y0, const Matrix &k1,
                              const Matrix &k2, double h) {
  auto yf = y.flat();
  auto y0f = y0.flat();
  auto k1f = k1.flat();
  auto k2f = k2.flat();
  for (size_t i = 0; i < yf.size(); ++i)
    yf[i] = y0f[i] + 0.5 * h * (k1f[i] + k2f[i]);
}

} // namespace tissue
