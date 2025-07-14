using Weave
ENV["TEXMFOUTPUT"] = joinpath(pwd(), "porocilo")
Weave.weave("dn3/scripts/demo.jl", doctype="minted2pdf", out_path="dn3/pdf")