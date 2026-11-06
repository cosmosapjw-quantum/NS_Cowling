! Reusable methods: callers supply the differential equations and their units.
module ns_numerical_methods
 use ns_bisection, only: bisect
 use ns_rk4, only: rk4_step
 use ns_adaptive_rk, only: adaptive_integrate
 use ns_numerics_types, only: root_options_t, root_result_t, ode_options_t, ode_result_t, step_result_t
 implicit none
 public
end module ns_numerical_methods
