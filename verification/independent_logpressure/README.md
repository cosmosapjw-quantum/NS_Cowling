# Independent log-pressure stellar check

This standalone Fortran module takes a cold barotrope in strictly increasing pressure order within each smooth phase. Pressure and energy density are in MeV per cubic femtometre. At a Maxwell interface include the same pressure twice, in low-to-high energy-density order. Supply positive pressures, a central pressure inside the grid, and a physical positive-slope barotrope. The first row defines the surface cutoff. The module uses RK4, resolves every EOS interval, and interpolates linearly in log pressure and log energy density. The independent variable is the square root of the logarithm of central pressure divided by local pressure.

Compile and run the supplied smooth-star check with GNU Fortran.

```sh
gfortran -O3 -ffree-line-length-none stellar.f90 smoke.f90 -o smoke
./smoke
```

The C-interoperable entry point is star(n,p,e,pc,frequency,offset,limit,res). The first calculation uses zero frequency and zero offset. The result contains gravitational mass in solar masses, radius in km, quadrupolar Love number, dimensionless tidal deformability, the central lapse-potential correction, the Cowling surface residual, and the outward W-node count. A nonzero positive trial frequency in Hz enables slow-conversion Cowling propagation. Pass the lapse correction from the equilibrium call as the offset. Determine mode identity with a separate converged two-sided solve and phasewise W and V counts. A production mode classification needs both displacement components. Recheck the pressure grid, cutoff and step cap for each EOS.
