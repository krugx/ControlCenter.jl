function grape_obj(eps_vec, prob)
  eps = reshape(eps_vec, (length(prob.shape[:, 1]), :))
  config = Policy(eps)
  S, F = action(config, prob; pwc=true)
  return S
end
function grape_grad!(storage_vec, eps_vec, prob)
  N = length(prob.tgrid) - 1
  dt = prob.tgrid[2] - prob.tgrid[1]
  # Cast to matrix
  eps = reshape(eps_vec, (length(prob.shape[:, 1]), :))
  config = Policy(eps)
  u = pulse(config, prob)
  storage = zeros(size(eps))

  ## Running cost
  for i in eachindex(prob.shape[:, 1])
    for j in 1:N
      storage[i, j] += prob.lambda_a * dt * eps[i, j]
    end
  end
  ## Final Cost
  for traj in prob.trajectories
    phi = copy(traj.ado)
    @assert phi isa SparseVector
    rho = Vector{SparseVector{ComplexF64}}(undef, N)
    for j in 1:N
      phi = rk4_step(phi, u[:, j], u[:, j], dt, traj.gen)
      rho[j] = copy(phi)
    end

    lambda = copy(traj.ado_tgt)
    for j in N:-1:1
      for i in eachindex(prob.shape[:, 1])
        storage[i, j] -= traj.weight * prob.shape[i, j] * dt * real(lambda' * traj.gen.prop_C[i].mat * rho[j])
      end
      lambda = rk4_step(lambda, u[:, j], u[:, j], dt, traj.gen; backwards=true)
    end

  end
  storage_vec[:] .= reshape(storage, :)
end

function grape_control(config::Policy, prob::ControlProblem)
  discretize_on_midpoints!(config::Policy, prob::ControlProblem)

  f = e -> grape_obj(e, prob)
  g! = (g, e) -> grape_grad!(g, e, prob)


  eps_vec_init = reshape(copy(config.eps), :)

  opt_st = Optim.Options(f_abstol=0e-5, g_tol=1e-8, iterations=prob.iterations, allow_f_increases=true, show_trace=true)

  res_st = Optim.optimize(f, g!, eps_vec_init, LBFGS(), opt_st)
  display(res_st)
  eps_vec = Optim.minimizer(res_st)

  eps = reshape(eps_vec, (length(prob.shape[:, 1]), :))
  config_opt = Policy(eps)


  return config_opt
end
