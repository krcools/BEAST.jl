struct DPVectorialTerm{A,K}
    alpha::A
    gamma::K
end

struct DPScalarTerm{A,K}
    alpha::A
    gamma::K
end

const DPTerm = Union{DPVectorialTerm, DPScalarTerm}

quaddata(pot::DPTerm, local_space, charts) = quadpoints(local_space, charts, (2,))
quadrule(pot::DPTerm, local_space, test_index, test_point, trial_index, trial_chart, quad_data) = quad_data[1,trial_index]

function kernelvals(pot::DPTerm, test_point, trial_neighborhood)

    γ = pot.gamma
    α = pot.alpha
    r = test_point - cartesian(trial_neighborhood)
    R = norm(r)

    γR = γ*R
    expn = exp(-γR)
    green = expn / (4π*R)
    gradgreen = -(γ + 1/R) * green / R * r

    return (vector=r, distance=R, green=green, gradgreen=gradgreen)
end

function integrand(pot::DPVectorialTerm, kernel_values, test_point, trial_values, trial_neighborhood)

    α = pot.alpha
    G = kernel_values.green
    f = trial_values.value

    return α*G*f
end

function integrand(pot::DPScalarTerm, kernel_values, test_point, trial_values, trial_neighborhood)

    α = pot.alpha
    ∇G = kernel_values.gradgreen
    f = trial_values.value

    return α*∇G*f
end
