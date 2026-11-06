! Public entry points. Every calculation calls the explicit physics modules.
module ns_stellar_methods
 use ns_cell_shooting, only: integrate_cells, match_cells, cell_run_t
 use ns_cell_eos, only: cell_state, cell_pressure_coordinate, table_pressure_coordinate
 use ns_table_eos, only: table_eos_t, table_audit_t, init_table_eos
 use ns_controlled_table_variant, only: build_controlled_variant, controlled_variant_audit_t
 implicit none
 public
end module ns_stellar_methods
