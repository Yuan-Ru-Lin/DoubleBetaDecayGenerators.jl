# Experiment for juliac

Compilation succeeded with the following command.

`julialauncher +1.12 --startup-file=no $JULIA_DEPOT_PATH/juliaup/julia-1.12.0-rc1+0.aarch64.apple.darwin14/share/julia/juliac/juliac.jl --output-exe a.out juliac_trial.jl`

It took me 1.5 minutes to compile and the output executable is of size 256MB. It ran.

Note that you may not set `$JULIA_DEPOT_PATH` and the Julia installation you have may not be `julia-1.12.0-rc1+0.aarch64.apple.darwin14`.

## Plan

For now this feature is too experimental to use in production. In order for external packages e.g. remage to utilize the generators defined in this repo, one should write a script that outputs an LH5 file where the table inside the file follows a TBD spec.
