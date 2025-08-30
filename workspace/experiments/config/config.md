# Understanding configuration parameters

Configuration parameters for the MPC-CBF Formation Control system. Each parameter group serves a specific purpose in the control system.

```json
{
    "mpc_params": {
        "h": 0.1,        // MPC discretization timestep (sec) - prediction temporal resolution
        "Ts": 0.01,      // Control timestep (sec) - must be ≤ h and h must be integer multiple of Ts  
        "k_hor": 16,     // Prediction horizon steps - total prediction time = k_hor * h
                         // CONSTRAINT: k_hor ≤ (num_pieces * piece_max_parameter) / h + 1
        "mpc_tuning": {
            "w_pos_err": 10,  // Position error weight - increase for tighter reference tracking
            "w_u_eff": 1,     // Control effort weight - increase for smoother trajectories (less jerky motion)
            "spd_f": 8        // Terminal cost horizon - position penalty on last spd_f steps, must be ≤ k_hor
        }
    },
    "physical_limits": {
        "p_min": [-5, -5],           // Workspace position bounds [x_min, y_min] 
        "p_max": [5, 5],             // Workspace position bounds [x_max, y_max]
        "v_min": [-2, -2, -2.618],   // Velocity limits [vx_min, vy_min, omega_min] TODO: change to vz_min
        "v_max": [2, 2, 2.618],      // Velocity limits [vx_max, vy_max, omega_max] TODO: change to vz_max
        "a_min": [-5, -5, -3.142],   // Acceleration limits [ax_min, ay_min, alpha_min] TODO: change to az_min
        "a_max": [5, 5, 3.142],      // Acceleration limits [ax_max, ay_max, alpha_max] TODO: change to az_max
        "pos_std": 0.001,            // Position measurement noise std (not used in MPC-CBF)
        "vel_std": 0.01              // Velocity measurement noise std (not used in MPC-CBF)
    },
    "pid_params": {                  // ⚠️  NOT USED in MPC-CBF example, but in CBF example (as a nominal control alternative)
        "kp": 3,                     // Proportional gain (would be used in PID controller)
        "ki": 0.1,                   // Integral gain (would be used in PID controller) 
        "kd": 0.3                    // Derivative gain (would be used in PID controller)
    },
    "robot_params": {
        "collision_shape": {               // ⚠️ Currently not used -> safety is enforced via CBF (d_min), Voronoi-based geometric collision avoidance is removed
            "aligned_box": [0.2, 0.2, 0],  // Robot bounding box half-extents [x_half, y_half, z_half]
            "radius": 0.5                  // Robot collision radius (alternative to box model)
        }
    },
    "cbf_params": {                       // Control Barrier Function parameters for safety/connectivity
        "d_min": 2.0,                     // Min safety distance - robots must maintain this distance
        "d_max": 4.0,                     // Max connectivity distance - robots must connect when farther  
        "cbf_horizon": 2,                 // CBF prediction horizon - must be ≤ k_hor
        "impc_iter": 2,                   // Iterative MPC-CBF iterations - more = better constraints but slower
        "slack_config": {                 // Separate slack controls for different constraint types
            "safety_slack": true,         // Allow safety constraint violations (collision avoidance)
            "clf_slack": false,           // Allow CLF constraint violations (formation performance)  
            "connectivity_slack": true,   // Allow connectivity constraint violations (d_min/d_max)
            "safety_slack_cost": 100000,  // Slack cost for safety violations - high cost to avoid collisions
            "clf_slack_cost": 50000,      // Slack cost for CLF violations
            "connectivity_slack_cost": 25000, // Slack cost for connectivity violations - formation flexibility
            "slack_decay_rate": 0.1       // Slack decay rate ∈ (0,1] - lower = longer violation tolerance
        }
    },
    "bezier_params": {                    // Piecewise Bezier trajectory representation
        "num_pieces": 3,                  // Number of Bezier pieces - increase if k_hor parameter error occurs
        "num_control_points": 4,          // Control points per piece (4 = cubic Bezier)
        "piece_max_parameter": 0.5,       // Parameter span per piece - total range = num_pieces * piece_max_parameter
        "bezier_continuity_upto_degree": 3 // Continuity degree between pieces (3 = C³ smooth)
    }
}
```

## Common Tuning Scenarios

**🔧 Trajectory too jerky/unsmooth:** Increase `w_u_eff`, reduce acceleration limits, or increase `bezier_continuity_upto_degree`

**🔧 Poor reference tracking:** Increase `w_pos_err` or `spd_f` for stronger position penalties  

**🔧 System too sluggish:** Decrease `w_u_eff`, increase velocity/acceleration limits, or reduce `k_hor`

**🔧 Parameter range error:** Either reduce `k_hor ≤ (num_pieces * piece_max_parameter) / h + 1` or increase `num_pieces`/`piece_max_parameter`

**🔧 Formation too tight/loose:** Adjust `d_min`/`d_max` connectivity distances

**🔧 Constraint violations:** Increase specific `slack_cost` values, disable relevant slack types, or increase `cbf_horizon` for stricter safety
