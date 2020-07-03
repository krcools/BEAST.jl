using BEAST

Δt, Nt = 1.0, 200
T0 = timebasisshiftedlagrange(Δt, Nt, 0)
T1 = timebasisshiftedlagrange(Δt, Nt, 1)
T2 = timebasisshiftedlagrange(Δt, Nt, 2)
T3 = timebasisshiftedlagrange(Δt, Nt, 3)
iT0 = integrate(T0)
iT1 = integrate(T1)
iT2 = integrate(T2)
iT3 = integrate(T3)
δ = timebasisdelta(Δt, Nt)

Id = BEAST.Identity()
Z = assemble(Id, δ, iT1)
Z = reshape(Z,1,1,length(Z))

taxis = range(0, step=Δt, length=Nt)
width, delay = 2.0, 6.0
f(t) = -2*(t-delay)/width^2*exp(-((t-delay)/width)^2)
F(t) = exp(-((t-delay)/width)^2)

b = F.(taxis)'
u = marchonintime(inv(Z[:,:,1]), Z, b, Nt)

plot(taxis, u[1,:])

taxis2 = range(-3Δt, 5Δt, length=200)
plot()
plot!(taxis2, iT3.(taxis2))

kmax = 2
cv = BEAST.ConvOp(
    reshape(Z[1,1,1:kmax],kmax,1,1),
    fill(1,(1,1)),
    fill(kmax,(1,1)),
    fill(1.0,(1,1)),
    Nt)

w, Q = BEAST.polyeig(cv)

u2 = marchonintime(inv(Z[:,:,1]), cv, b, Nt)