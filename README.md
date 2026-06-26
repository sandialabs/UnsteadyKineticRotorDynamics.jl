# UnsteadyKineticRotorDynamics.jl

**Summary**: Unsteady kinetic rotor dynamics for wind and hydrokinetic rotors,
implemented as low-pass time-domain dynamics around steady CCBlade solves.

**Author**: Kevin R. Moore and contributors

**Features**:

- Time-domain low-pass filtering of CCBlade outputs
- Optional OpenFAST inflow (OWENSOpenFASTWrappers)
- Compatible with AD tools when `ifw=false`

**Installation**:

```julia
using Pkg
Pkg.add(url = "https://github.com/kevmoor/UnsteadyKineticRotorDynamics.jl")
```

For local development from a checkout:

```julia
using Pkg
Pkg.develop(path = "/path/to/UnsteadyKineticRotorDynamics.jl")
```

**Quick Start**:

```julia
import CCBlade
import UnsteadyKineticRotorDynamics

# rotor/sections/ops setup using CCBlade
params = UnsteadyKineticRotorDynamics.UnsteadyParams(0.5, 2.0)
state = UnsteadyKineticRotorDynamics.UnsteadyState(sections, ops; V_wake_old=Vinf)

for k = 1:nsteps
    out = UnsteadyKineticRotorDynamics.unsteady_step!(state, rotor, sections, ops, params; dt=dt)
end
```

**Run Unit Tests**:

```julia
pkg> activate .
pkg> test
```

**Documentation**:

```julia
julia --project=docs docs/make.jl
```

Sources live in `docs/src/`.

**Notes**:

- This package depends on `CCBlade.jl` for steady solves.
- `OWENSOpenFASTWrappers` is only needed when `ifw=true`.
