using BEAST

sol = 1.0
Δt, Nt = 0.12/sol,400
@show tmax = (Nt-1)*Δt
T = timebasisshiftedlagrange(Δt, Nt, 3)

timeels, timead = assemblydata(T)
