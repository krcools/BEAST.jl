using .Variational



mutable struct DiscreteEquation
  equation
  trial_space_dict # dictionary mapping indices into trial space to FE spaces
  test_space_dict  # dictionary mapping indices into test space to FE spaces
end


function discretise(eq, space_mappings::Pair...)
    trial_space_dict = Dict()
    test_space_dict = Dict()
    for sm in space_mappings

        found = false
        sm.first.space == eq.lhs.trial_space && (dict = trial_space_dict; found = true)
        sm.first.space == eq.lhs.test_space  && (dict = test_space_dict;  found = true)
        @assert found "Vector $(sm.first) neither in test nor in trial space"

        @assert !haskey(dict, sm.first.idx) "multiple mappings for $(sm.first)"
        dict[sm.first.idx] = sm.second
    end

    # check that all symbols where mapped
    for p in 1:length(eq.lhs.trial_space) @assert haskey(trial_space_dict,p) end
    for p in 1:length(eq.lhs.test_space)  @assert haskey(test_space_dict, p) end

    DiscreteEquation(eq, trial_space_dict, test_space_dict)
end


"""
    discr(eq, pairs...)

This macro provides syntactical sugar for the definition of a discretisation
of a varational formulation. Given a variational equation EQ: Find j ∈ X such
that for all k ∈ Y a(k,j) = f(k) can be discretised by stating:

    eq = @discretise EQ j∈X k∈Y
"""
macro discretise(eq, pairs...)
    r = :(BEAST.discretise($eq))
    for p in pairs
        x = p.args[2]
        X = p.args[3]
        push!(r.args, :($x=>$X))
    end
    return esc(r)
end


sysmatrix(eq::DiscreteEquation) = assemble(eq.equation.lhs, eq.test_space_dict, eq.trial_space_dict)
rhs(eq::DiscreteEquation) = assemble(eq.equation.rhs, eq.test_space_dict)

function assemble(lform::LinForm, test_space_dict)

    terms = lform.terms
    T = ComplexF64

    I = Int[1]
    for p in 1:length(lform.test_space)
        X = test_space_dict[p]
        push!(I, last(I) + numfunctions(X))
    end
    B = zeros(T, last(I)-1)

    for t in terms

        α = t.coeff
        a = t.functional
        m = t.test_id
        X = test_space_dict[m]
        o = t.test_ops

        # act with the various ops on X
        for op in reverse(o)
            Y = X;
            X = op[end](op[1:end-1]..., Y)
        end

        b = assemble(a, X)
        B[I[m] : I[m+1]-1] = b
    end

    return B
end

function assemble(bilform::BilForm, test_space_dict, trial_space_dict)

  lhterms = bilform.terms
  T = ComplexF64 # TDOD: Fix this

  # determine the offsets of the different blocks in the sys matrix
  I = Int[1]
  J = Int[1]

  for p in 1:length(bilform.test_space)
    X = test_space_dict[p]
    push!(I, last(I) + numfunctions(X))
  end

  for q in 1:length(bilform.trial_space)
    Y = trial_space_dict[q]
    push!(J, last(J) + numfunctions(Y))
  end

  # allocate the memory for the matrices
  Z = zeros(T, last(I)-1, last(J)-1)

  # For each block, compute the interaction matrix
  for t in lhterms

      α = t.coeff
      a = t.kernel

      m = t.test_id
      x = test_space_dict[m]
      for op in reverse(t.test_ops)
          x = op[end](op[1:end-1]..., x)
      end

      n = t.trial_id
      y = trial_space_dict[n]
      for op in reverse(t.trial_ops)
          y = op[end](op[1:end-1]..., y)
      end

      r = I[m] : (I[m+1] - 1)
      c = J[n] : (J[n+1] - 1)

      z = assemble(a, x, y)
      Z[r,c] += α * z
  end

  return Z
end
