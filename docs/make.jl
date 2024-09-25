using SchrodingerToolkit
using Documenter

DocMeta.setdocmeta!(SchrodingerToolkit, :DocTestSetup, :(using SchrodingerToolkit);
                    recursive = true)

makedocs(;
         modules = [SchrodingerToolkit],
         authors = "walexaindre <walexaindre@hotmail.com
       > and contributors",
         sitename = "SchrodingerToolkit.jl",
         format = Documenter.HTML(;
                                  canonical = "https://walexaindre.github.io/SchrodingerToolkit.jl",
                                  edit_link = "",
                                  assets = String[],),
         pages = ["Package Overview" => "index.md",
                  "Usage" => ["Backend selection" => "./Usage/backend.md",
                              "Method selection" => "./Usage/methods.md",
                              "How to define a problem" => "./Usage/construction.md"],
                  "Examples" => ["Notes" => "./Examples/start.md"
                                 "Examples without Josephson junction and trapping potential" => ["2D" => "./Examples/NoTrappingPotentialNoJunction/2D.md"]
                                 "Examples without trapping potential" => ["2D" => "./Examples/NoTrappingPotential/2D.md"]],
                  "Plotting recipes" => ["Prelude" => "./Plot/howto.md"],
                  "Extending the package" => ["Prelude" => "./Extensions/prelude.md",
                                              "Adding a new finite difference scheme" => "./Extensions/fdschemes.md"],
                  "API" => ["Index" => "./API/index.md",
                            "Types" => "./API/API.md"],
                  "Theory" => ["Problem" => "./Theory/Problem.md",
                               "DFP" => "./Theory/DFP.md"],
                  "About" => "about.md"],)

deploydocs(;
           repo = "github.com/walexaindre/SchrodingerToolkit.jl.git",
           devbranch = "main",)
