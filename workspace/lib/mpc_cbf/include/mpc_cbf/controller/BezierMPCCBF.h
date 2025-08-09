//
// Created by lishuo on 9/21/24.
//

#ifndef MPC_CBF_BEZIERMPCCBF_H
#define MPC_CBF_BEZIERMPCCBF_H

#include <math/collision_shapes/CollisionShape.h>
#include <mpc_cbf/optimization/PiecewiseBezierMPCCBFQPGenerator.h>
#include <qpcpp/Problem.h>
#include <qpcpp/solvers/CPLEX.h>

namespace mpc_cbf {
template <typename T, unsigned int DIM>
class BezierMPCCBF {
  public:
    using Params = typename mpc_cbf::PiecewiseBezierMPCCBFQPOperations<T, DIM>::Params;
    using DoubleIntegrator =
        typename mpc_cbf::PiecewiseBezierMPCCBFQPOperations<T, DIM>::DoubleIntegrator;
    using PiecewiseBezierMPCCBFQPOperations = mpc_cbf::PiecewiseBezierMPCCBFQPOperations<T, DIM>;
    using PiecewiseBezierMPCCBFQPGenerator = mpc_cbf::PiecewiseBezierMPCCBFQPGenerator<T, DIM>;
    using FovCBF = typename mpc_cbf::PiecewiseBezierMPCCBFQPOperations<T, DIM>::FovCBF;
    using State = typename mpc_cbf::PiecewiseBezierMPCCBFQPOperations<T, DIM>::State;
    using TuningParams = mpc::TuningParams<T>;

    using SingleParameterPiecewiseCurve = splines::SingleParameterPiecewiseCurve<T, DIM>;
    using VectorDIM = math::VectorDIM<T, DIM>;
    using Vector = math::Vector<T>;
    using Matrix = math::Matrix<T>;
    using CollisionShape = math::CollisionShape<T, DIM>;
    using Problem = qpcpp::Problem<T>;
    using CPLEXSolver = qpcpp::CPLEXSolver<T>;
    using SolveStatus = qpcpp::SolveStatus;

    BezierMPCCBF(Params& p, std::shared_ptr<DoubleIntegrator> model_ptr,
                 std::shared_ptr<FovCBF> fov_cbf_ptr,
                 std::shared_ptr<const CollisionShape> collision_shape_ptr);
    ~BezierMPCCBF() = default;

    bool optimize(SingleParameterPiecewiseCurve& result_curve, const State& current_state,
                  const std::vector<VectorDIM>& other_robot_positions, const Vector& ref_positions);

    Vector generatorDerivativeControlInputs(uint64_t derivative_degree);

  private:
    TuningParams mpc_tuning_;
    Vector ts_samples_;

    PiecewiseBezierMPCCBFQPGenerator qp_generator_;
    // the derivative degree that the resulting trajectory must be
    // continuous upto.
    uint64_t bezier_continuity_upto_degree_;
    std::shared_ptr<const CollisionShape> collision_shape_ptr_;
};

} // namespace mpc_cbf

#endif // MPC_CBF_BEZIERMPCCBF_H
