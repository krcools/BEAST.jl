
# [Contributing](@id contribRef)

In order to contribute to this package directly create a pull request against the `master` branch. Before doing so please: 

- Follow the style of the surrounding code.
- Supplement the documentation.
- Write tests and check that no errors occur.


---
## Style

For a consistent style the [JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl) package is used which enforces the style defined in the *.JuliaFormatter.toml* file. To follow this style simply run
```julia
using JuliaFormatter
format(pkgdir(BEAST))
```

!!! note
    That all files follow the JuliaFormatter style is tested during the unit tests. Hence, do not forget to execute the two lines above. Otherwise, the tests are likely to not pass.


---
## Documentation

Add documentation for any changes or new features following the style of the existing documentation. For more information you can have a look at the [Documenter.jl](https://documenter.juliadocs.org/stable/) documentation.


---
## [Tests](@id tests)

Write tests for your code changes and verify that no errors occur, e.g., by running
```julia
using Pkg
Pkg.test("BEAST")
```

For more detailed information on which parts are tested the coverage can be evaluated on your local machine, e.g., by
```julia
using Pkg
Pkg.test("BEAST"; coverage=true, julia_args=`--threads 4`)

# determine coverage
using Coverage
src_folder = pkgdir(BEAST) * "/src"
coverage   = process_folder(src_folder)
LCOV.writefile("path-to-folder-you-like" * "BEAST.lcov.info", coverage)

clean_folder(src_folder) # delete .cov files

# extract information about coverage
covered_lines, total_lines = get_summary(coverage)
@info "Current coverage:\n$covered_lines of $total_lines lines ($(round(Int, covered_lines / total_lines * 100)) %)"
```

In Visual Studio Code the [Coverage Gutters](https://marketplace.visualstudio.com/items?itemName=ryanluker.vscode-coverage-gutters) plugin can be used to visualize the tested lines of the code by inserting the path of the *BEAST.lcov.info* file in the settings.

!!! note
    From Julia 1.11 onwards the the coverage can be displayed in Visual Studio Code directly, since the [TestItemRunner.jl](https://www.julia-vscode.org/docs/stable/userguide/testitems/#Code-coverage) package is employed.