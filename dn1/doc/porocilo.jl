using Weave
ENV["TEXMFOUTPUT"] = joinpath(pwd(), "porocilo")
ENV["TEXMFOUTPUT"] = joinpath(pwd(), "porocilo")
Weave.weave("dn1/scripts/demo.jl", doctype="minted2pdf", out_path="dn1/pdf")