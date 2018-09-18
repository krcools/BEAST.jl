function show_which_quadrule(op,X,Y)

    test_charts, test_ad = assemblydata(X)
    trial_charts, trial_ad = assemblydata(Y)

    x = refspace(X)
    y = refspace(Y)

    @show @which quaddata(op, x, y, test_charts, trial_charts)
    qd = quaddata(op, x, y, test_charts, trial_charts)
    @show @which quadrule(op, x, y, 1, first(test_charts), 1, first(trial_charts), qd)


    nothing
end
