struct Generator
  delta::Float64
  prop_0::HEOMPropagator
  prop_C::Vector{HEOMPropagator}
end

struct Trajectory
  ado::SparseVector
  ado_tgt::SparseVector
  weight::Float64
  gen::Generator

  function Trajectory(rho::Matrix{ComplexF64}, rho_tgt::Matrix{ComplexF64}, weight, gen::Generator)
    hdim = size(rho, 1)
    dim = size(gen.prop_0.mat, 1)
    ado = spzeros(ComplexF64, dim)
    ado[1:hdim^2] = reshape(rho, (hdim^2, 1))
    ado_tgt = spzeros(ComplexF64, dim)
    ado_tgt[1:hdim^2] = reshape(rho_tgt, (hdim^2, 1))

    new(ado, ado_tgt, weight, gen)
  end
  function Trajectory(ado::SparseVector, ado_tgt::SparseVector, weight, gen::Generator)
    new(ado, ado_tgt, weight, gen)
  end
end


struct ControlProblem
  tgrid::Vector{Float64}
  lambda_a::Float64
  iterations::Int
  shape::Matrix{Float64}
  trajectories::Vector{Trajectory}
  isparametrized::Bool
  function ControlProblem(tgrid, lambda_a, iterations, shape, trajectories; isparametrized=false)
    new(tgrid, lambda_a, iterations, shape, trajectories, isparametrized)
  end
end

abstract type Config end

struct Policy <: Config
  eps::Matrix{Float64}
  function Policy(eps::Matrix{Float64})
    new(eps)
  end
end

struct Parametrization <: Config
  A::Matrix{Float64}
  tau::Matrix{Float64}
  sigma::Matrix{Float64}
end

