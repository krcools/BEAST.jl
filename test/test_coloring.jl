using Test
using CompScienceMeshes, BEAST
using OhMyThreads, GraphsColoring

@testset "Coloring tests" begin

    for m in [readmesh(joinpath(dirname(pathof(BEAST)), "../examples/sphere2.in")), meshrectangle(1.0, 1.0, 0.1)]

        for X in [
            raviartthomas(m),
            n × raviartthomas(m),
            BEAST.raviartthomas2(m),
            buffachristiansen(m),
            lagrangec0d1(m),
            lagrangec0d2(m),
            BEAST.gwpdiv(m; order=2),
            BEAST.gwpcurl(m; order=2),
            BEAST.gwpdiv(m; order=3),
            BEAST.gwpcurl(m; order=3),
            brezzidouglasmarini(m),
            BEAST.ncrossbdm(m),
            BEAST.nedelec2(m),
            BEAST.nedelec(m)
        ]
            elements, ad, _ = assemblydata(X; onlyactives=true)

            g = GraphsColoring.conflictmatrix(X; addata=(elements, ad, nothing))

            colors = BEAST.color(X, DynamicScheduler())

            for i in eachindex(colors)
                @test iszero(nnz(g[colors[i], colors[i]]))
            end

            for color in colors
                basisindices = Int[]
                els = Int[]

                for element in color

                    localbasisindices = Set{Int}()
                    for i in 1:length(ad[element])
                        for (m, a) in ad[element, i]
                            iszero(a) && continue
                            push!(localbasisindices, m)
                            push!(els, element)
                        end
                    end
                    append!(basisindices, collect(localbasisindices))
                end

                sort!(basisindices)
                @test unique(basisindices) == basisindices
            end
        end
    end

    m = tetmeshsphere(1.0, 0.2)
    X = nedelecd3d(m)

    for X in [nedelecd3d(m), nedelecc3d(m), raviartthomas(m), brezzidouglasmarini3d(m)]
        elements, ad, _ = assemblydata(X; onlyactives=true)

        g = GraphsColoring.conflictmatrix(X; addata=(elements, ad, nothing))

        colors = BEAST.color(X, DynamicScheduler())

        for i in eachindex(colors)
            @test iszero(nnz(g[colors[i], colors[i]]))
        end

        for color in colors
            basisindices = Int[]
            els = Int[]

            for element in color

                localbasisindices = Set{Int}()
                for i in 1:length(ad[element])
                    for (m, a) in ad[element, i]
                        iszero(a) && continue
                        push!(localbasisindices, m)
                        push!(els, element)
                    end
                end
                append!(basisindices, collect(localbasisindices))
            end

            sort!(basisindices)
            @test unique(basisindices) == basisindices
        end
    end
end


@testitem "cellcoloring vs dofsplitting" begin
    using CompScienceMeshes, OhMyThreads, LinearAlgebra

    fn = dirname(pathof(BEAST)) * "/../test/assets/sphere45.in"
    m = readmesh(fn)

    T = Maxwell3D.singlelayer(wavenumber=3.0)
    X = raviartthomas(m)
    Y = buffachristiansen(m)

    Txx1 = assemble(T, X, X;
        threading=:cellcoloring,
        scheduler=OhMyThreads.DynamicScheduler())
    Txx2 = assemble(T, X, X;
        threading=:dofsplitting,
        scheduler=OhMyThreads.DynamicScheduler())
    @test norm(Txx1 - Txx2) < 1e-12

    Tyy1 = assemble(T, Y, Y;
        threading=:cellcoloring,
        scheduler=OhMyThreads.DynamicScheduler())
    Tyy2 = assemble(T, Y, Y;
        threading=:dofsplitting,
        scheduler=OhMyThreads.DynamicScheduler())
    @test norm(Tyy1 - Tyy2) < 1e-12
end