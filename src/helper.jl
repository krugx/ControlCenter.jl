function action(conf::Config, prob::ControlProblem; pwc=false)
  N = length(prob.tgrid) - 1
  dt = prob.tgrid[2] - prob.tgrid[1]
  F_T = 0.0
  J_R = 0.0
  u = pulse(conf, prob)

  for i in 1:N
    if pwc
      J_R += dt / 2 * u[:, i]' * prob.lambda_a * u[:, i]
    else
      J_R += dt / 2 * (u[:, i]' + u[:, i+1]') / 2 * prob.lambda_a * (u[:, i] + u[:, i+1]) / 2
    end
  end

  for traj in prob.trajectories
    phi = copy(traj.ado)
    for i in 1:N
      if pwc
        phi = rk4_step(phi, u[:, i], u[:, i], dt, traj.gen)
      else
        phi = rk4_step(phi, u[:, i], u[:, i+1], dt, traj.gen)
      end
    end
    traj_prop = Trajectory(phi, traj.ado_tgt, traj.weight, traj.gen)
    F_T += fidelity(traj_prop)
  end

  return 1 - F_T + J_R, 1 - F_T
end

function fidelity(traj::Trajectory)
  return traj.weight * real(traj.ado_tgt' * traj.ado)
end

function pulse(conf::Policy, prob::ControlProblem)
  @assert size(prob.shape) == size(conf.eps)
  return conf.eps .* prob.shape
end

function pulse(conf::Parametrization, prob::ControlProblem)
  C = length(conf.A[:, 1])
  m = length(conf.A[1, :])
  function u(t)
    res = zeros(C)
    for k in 1:C
      for j in 1:m
        res[k] += conf.A[k, j] * exp(-(t - conf.tau[k, j])^2 / conf.sigma[k, j]^2)
      end
    end
    return res
  end
  ugrid = u.(prob.tgrid)
  umat = reduce(vcat, transpose.(ugrid))'
  return umat .* prob.shape
end

function H_total(u, gen)
  H_tot = gen.prop_0.mat
  for i in 1:length(gen.prop_C)
    H_tot = H_tot + u[i] * gen.prop_C[i].mat
  end
  return H_tot
end

function rk4_step(phi, u, w, dt, gen::Generator; backwards=false)
  Hfwd(v) = H_total(v, gen)
  Hrev(v) = H_total(v, gen)'

  H = backwards ? Hrev : Hfwd

  k1 = H(u) * phi
  k2 = H((u + w) / 2) * (phi + k1 * dt / 2)
  k3 = H((u + w) / 2) * (phi + k2 * dt / 2)
  k4 = H(w) * (phi + k3 * dt)

  return phi + dt / 6 * (k1 + 2k2 + 2k3 + k4)
end


function discretize_on_midpoints!(config::Policy, prob::ControlProblem)
  if length(config.eps[1, :]) != length(prob.tgrid) - 1
    eps_new = copy(config.eps)
    shape_new = copy(prob.shape)
    eps_new[:, 1] = config.eps[:, 1]
    eps_new[:, end] = config.eps[:, end]
    shape_new[:, 1] = prob.shape[:, 1]
    shape_new[:, end] = prob.shape[:, end]
    for i = 2:length(prob.tgrid)-1
      eps_new[:, i] = 2 * config.eps[:, i] - eps_new[:, i-1]
      shape_new[:, i] = 2 * prob.shape[:, i] - shape_new[:, i-1]
    end
    config.eps .= eps_new
    prob.shape .= shape_new
  end
end
