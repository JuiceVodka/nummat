using Weave
ENV["TEXMFOUTPUT"] = joinpath(pwd(), "porocilo")
Weave.weave("dn2/scripts/demo.jl", doctype="minted2pdf", out_path="dn2/pdf")