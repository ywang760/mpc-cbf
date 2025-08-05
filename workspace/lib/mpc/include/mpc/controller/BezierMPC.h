//
// Created by lishuo on 8/22/24.
//

#ifndef MPC_BEZIERMPC_H
#define MPC_BEZIERMPC_H

#include <cmath>
#include <math/Helpers.h>
#include <math/collision_shapes/CollisionShape.h>
#include <mpc/optimization/PiecewiseBezierMPCQPGenerator.h>
#include <qpcpp/Problem.h>
#include <qpcpp/solvers/CPLEX.h>
#include <separating_hyperplanes/Voronoi.h>

namespace mpc {

template <typename T, unsigned int DIM>
class BezierMPC {
  public:
    using Params = typename mpc::PiecewiseBezierMPCQPOperations<T, DIM>::Params;
    using State = model::State<T, DIM>;
    using DoubleIntegrator = model::DoubleIntegrator<T, DIM>;
    using StatePropagator = model::StatePropagator<T>;
    using SingleParameterPiecewiseCurve = splines::SingleParameterPiecewiseCurve<T, DIM>;
    using PiecewiseBezierMPCQPOperations = mpc::PiecewiseBezierMPCQPOperations<T, DIM>;
    using PiecewiseBezierMPCQPGenerator = mpc::PiecewiseBezierMPCQPGenerator<T, DIM>;
    using Problem = qpcpp::Problem<T>;
    using CPLEXSolver = qpcpp::CPLEXSolver<T>;
    using SolveStatus = qpcpp::SolveStatus;
    using Vector = math::Vector<T>;
    using Matrix = math::Matrix<T>;
    using Hyperplane = math::Hyperplane<T, DIM>;
    using VectorDIM = math::VectorDIM<T, DIM>;
    using AlignedBox = math::AlignedBox<T, DIM>;
    using CollisionShape = math::CollisionShape<T, DIM>;

    BezierMPC(Params& p, std::shared_ptr<DoubleIntegrator> model_ptr,
              uint64_t bezier_continuity_upto_degree,
              std::shared_ptr<const CollisionShape> collision_shape_ptr);
    ~BezierMPC() = default;

    bool optimize(SingleParameterPiecewiseCurve& result_curve, const State& current_state,
                  const std::vector<VectorDIM>& other_robot_positions, const Vector& ref_positions);

    Vector generatorDerivativeControlInputs(uint64_t derivative_degree);

  private:
    TuningParams<T> mpc_tuning_;
    std::shared_ptr<DoubleIntegrator> model_ptr_;
    StatePropagator A0_;
    StatePropagator Lambda_;
    Vector ts_samples_;

    PiecewiseBezierMPCQPGenerator qp_generator_;

    // the derivative degree that the resulting trajectory must be
    // continuous upto.
    uint64_t bezier_continuity_upto_degree_;
    std::shared_ptr<const CollisionShape> collision_shape_ptr_;
};

} // namespace mpc

#endif //MPC_BEZIERMPC_H
