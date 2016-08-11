ENV["JULIA_PKGDIR"] = dirname(@__FILE__) * "/Packages"

lines = ""
case = nothing
meshcase = nothing
for case in [1,2,3,4]
    lines *= "\\hline \n"
    for meshcase in [3,1,2,5,6]
        println("meshcase=$meshcase")
        include("fielddiary.jl")
        println("N = $(size(points,2)) \n")

        line = "&$a $b $c & $(size(points,2)) & $(round(mean(rpsi),4)) & $(round(sum(abs(psi - psit))/sum(abs(psit))*100,4)) & $(round(quantile(rpsi,1/4),4)) & $(round(quantile(rpsi,1/2),4)) & $(round(quantile(rpsi,3/4),4)) \\\\ \n"
        lines = lines*line
    end
end

println(lines)



lines = ""
case = nothing
meshcase = nothing
for case in [5,6,7,8,9]
    lines *= "\\hline \n"
    for meshcase in [3,1,2,5,6]
        println("meshcase=$meshcase")
        include("fielddiary.jl")
        println("N = $(size(points,2)) \n")

        line = "&$a $b $c & $(size(points,2)) & $(round(mean(rHn),4)) & $(round(sum(abs(Hn-Htn))/sum(abs(Htn))*100,4))  & $(round(quantile(rHn,1/4),4)) & $(round(quantile(rHn,1/2),4)) & $(round(quantile(rHn,3/4),4)) \\\\ \n"
        lines = lines*line
    end
end

println(lines)


Htlines = ""
Hnlines = ""

case = 10
meshcase = nothing

for meshcase in [3,1,2,5,6]
    println("meshcase=$meshcase")
    include("fielddiary.jl")
    println("N = $(size(points,2)) \n")

    Hnline = "&$a $b $c & $(size(points,2)) & $(round(mean(rHn),4)) & $(round(sum(abs(Hn-Htn))/sum(abs(Htn))*100,4))  & $(round(quantile(rHn,1/4),4)) & $(round(quantile(rHn,1/2),4)) & $(round(quantile(rHn,3/4),4)) \\\\ \n"
    Hnlines *= Hnline

    Htline = "&$a $b $c & $(size(points,2)) & $(round(mean(rHt),4)) & $(round(sum(abs(Ht-Htt))/sum(abs(Htt))*100,4))   & $(round(quantile(rHt,1/4),4)) & $(round(quantile(rHt,1/2),4)) & $(round(quantile(rHt,3/4),4)) \\\\ \n"
    Htlines *= Htline
    
end

println("\\hline \n"*Htlines*"\\hline \n"*Hnlines)
