//
// Created by lishuo on 3/31/24.
//

#include <math/Types.h>
#include <math/Helpers.h>
#include <splines/curves/SingleParameterCurve.h>
#include <splines/curves/Bezier.h>
#include <splines/curves/SingleParameterPiecewiseCurve.h>
#include <splines/curves/Helpers.h>
#include <splines/optimization/SingleParameterPiecewiseCurveQPGenerator.h>
#include <splines/optimization/BezierQPOperations.h>
#include <qpcpp/solvers/CPLEX.h>
#include <iostream>

int main() {
    constexpr unsigned int DIM = 2U;
    using VectorDIM = math::VectorDIM<double, DIM>;
    using AlignedBox = math::AlignedBox<double, DIM>;
    using SingleParameterPiecewiseCurve = splines::SingleParameterPiecewiseCurve<double, DIM>;
    using Bezier = splines::Bezier<double, DIM>;
    using SingleParameterCurve = splines::SingleParameterCurve<double, DIM>;
    using SingleParameterPiecewiseCurveQPGenerator = splines::SingleParameterPiecewiseCurveQPGenerator<double, DIM>;

    // define the maximum parameter
    double max_parameter = 1.;
    // create single parameter curves
    std::vector<VectorDIM> control_pts_0 = {{0,0}, {2,2}, {5,2}, {7,0}};
    std::unique_ptr<Bezier> bezier_ptr_0 = std::make_unique<Bezier>(max_parameter, control_pts_0);
    std::vector<VectorDIM> control_pts_1 = {{7,0}, {9,-2}, {12,-2}, {14,0}};
//    std::vector<VectorDIM> control_pts_1 = {{7,0}, {8,-2}, {12,-2}, {14,0}};
    std::unique_ptr<Bezier> bezier_ptr_1 = std::make_unique<Bezier>(max_parameter, control_pts_1);
    //
    std::cout << "bezier 1 eval at t = 0.5: " << bezier_ptr_1->eval(0.5, 0) << "\n";
    // append them to be the piecewise curve
    SingleParameterPiecewiseCurve piecewise_curve;
    piecewise_curve.addPiece(std::move(bezier_ptr_0));
    piecewise_curve.addPiece(std::move(bezier_ptr_1));

    // log the control points to json file
    std::string json_filename = "../tools/python/bezier_piecewise.json";
    splines::writePiecewiseCurveToJson(piecewise_curve, json_filename);

    // test piecewise curve
    std::cout << "number of pieces: " << piecewise_curve.numPieces() << "\n";
    std::cout << "piecewise eval at t = 1: " << piecewise_curve.eval(1, 0) << "\n";
    std::cout << "piecewise eval at t = 1.5: " << piecewise_curve.eval(1.5, 0) << "\n";
    VectorDIM derivative_t = piecewise_curve.eval(1.0, 1);
    VectorDIM derivative_t_ = piecewise_curve.eval(1 + 1e-10, 1);
    std::cout << "piecewise eval 1deg at t = 1.: " << derivative_t << "\n";
    std::cout << "piecewise eval 1deg at t = 1+e-10: " << derivative_t_ << "\n";
    std::cout << "derivative at joint point of piecewise curve is similar: " << math::isApproximatelyEqual<double, DIM>(derivative_t, derivative_t_, 1e-3);

    // compose the QP problem
    using BezierQPOperations = splines::BezierQPOperations<double, DIM>;
    using Problem = qpcpp::Problem<double>;
    using CPLEXSolver = qpcpp::CPLEXSolver<double>;
    using SolveStatus = qpcpp::SolveStatus;

    size_t bezier_degree = 3;
    size_t num_control_points = bezier_degree + 1;
    BezierQPOperations::Params p({num_control_points, max_parameter});
    std::unique_ptr<BezierQPOperations> bezier_qp_ptr_0 = std::make_unique<BezierQPOperations>(p);
    std::unique_ptr<BezierQPOperations> bezier_qp_ptr_1 = std::make_unique<BezierQPOperations>(p);

    // ---- optimize for single curve ----
    SingleParameterPiecewiseCurveQPGenerator single_curve_qp_generator;
    single_curve_qp_generator.addPiece(std::move(bezier_qp_ptr_0));
    single_curve_qp_generator.addPiece(std::move(bezier_qp_ptr_1));

    // constraint0: bounding box constraints for control point positions
    VectorDIM workspace_min = {-1, -1};
    VectorDIM workspace_max = {10, 10};
    AlignedBox workspace(workspace_min, workspace_max);
    single_curve_qp_generator.addBoundingBoxConstraint(workspace);

    // constraint1: on the start and goal target
    VectorDIM start_0 = {0,0};
    VectorDIM goal_0 = {2,0};
    single_curve_qp_generator.addEvalConstraint(0, 0, start_0);
    single_curve_qp_generator.addEvalConstraint(1, 0, goal_0);
    single_curve_qp_generator.addEvalConstraint(1+1e-2, 0, goal_0);
    single_curve_qp_generator.addEvalConstraint(1.5, 0, {3, 2});
    single_curve_qp_generator.addEvalConstraint(2, 0, {4,0});

    // constraint2: meet at a middle point
    VectorDIM middle_0 = {1, 1};
    single_curve_qp_generator.addEvalConstraint(0.5, 0, middle_0);
    VectorDIM middle_1 = {1.5, 0.2};
    single_curve_qp_generator.addEvalConstraint(0.75, 0, middle_1);

    // optimize for the energy cost
    single_curve_qp_generator.addIntegratedSquaredDerivativeCost(1, 1);

    // access the optimization problem and solve it
    CPLEXSolver cplex_solver;
    Problem& opt_problem = single_curve_qp_generator.problem();
    SolveStatus solve_status = cplex_solver.solve(opt_problem);
    std::cout << "optimization status: " << SolveStatusToStr(solve_status) << "\n";

    // recover the curve and visualize it
    SingleParameterPiecewiseCurve result_curve;
    result_curve = single_curve_qp_generator.generateCurveFromSolution();

    // log the control points to json file
    splines::writePiecewiseCurveToJson(result_curve, json_filename);

    // eval resulting curve
    std::cout << "optimized curve eval(must pass middle point): " << result_curve.eval(0.5, 0) << "\n";
    std::cout << "optimized curve eval(must pass start point): " << result_curve.eval(0, 0) << "\n";
    std::cout << "optimized curve eval(must pass goal point): " << result_curve.eval(1, 0) << "\n";

    return 0;
}
