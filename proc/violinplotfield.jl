#case = 10
case = 10
meshcase = nothing

allerrors = Any[]

for meshcase in [1,5,6]
    println("meshcase=$meshcase")
    include("fielddiary.jl")
    println("N = $(size(points,2)) \n")

    allerrors = Any[allerrors...,copy(rHt),copy(rHn)]
end


using PyPlot; const plt = PyPlot


# using PyCall
# @pyimport matplotlib.pyplot as plt


plt.rc("text",usetex=true)
plt.rc("font", family="serif", serif="Times")
plt.rc("xtick",labelsize=10)
plt.rc("ytick",labelsize=10)
plt.rc("axes",labelsize=10)


width = 3.487 
height = width / 1.618
fig = plt.figure(figsize=(width,height))

ax = fig[:add_subplot](1, 1, 1)

ax[:set_xticks]([1,2,4,5,7,8])
ax[:set_xticklabels](["","","","","",""])


# n, p = 40,8
# d = map(abs,randn((p,n)))
# #d = d .+ log(collect(1:p)) * -5 + 10
# d = d'

# d = Any[d[:,1][1:30],d[:,2][1:10],d[:,3],d[:,4],d[:,5],d[:,6],d[:,7],d[:,8]]

#pos = [1,2,4,5,7,8]
pos = [1,2,7,8,4,5]

#pos = [1,2,3,5,6,7,9,10] # [1,2,3,4,5,6,7,8]

colorcode = ["darksage","steelblue","darksage","steelblue","darksage","steelblue"]

for i in 1:6
    ri = allerrors[i]
    plt.plot(ones(ri)*pos[i],ri,"+",color=colorcode[i])
end

#plt.plot([],[],"+",color="darksage",label=L"\boldsymbol n \cdot \boldsymbol H")

#plt.legend(loc=2)

#violinparts = plt.violinplot(allerrors, pos, points=40, widths=0.5,showmeans=false, showextrema=false, showmedians=false,bw_method="silverman")
#Any[[1,2,3],[4,55,5,5,6]]

# for (c,pc) in zip(colorcode,violinparts["bodies"])
#     pc[:set_facecolor](c)
#     pc[:set_edgecolor]("k")
# end

plt.ylim(0,10)
plt.xlim(0,9)

plt.ylabel("Releative error %")

plt.tight_layout()
fig[:savefig]("figures/violinplot.pdf")
plt.close("all")
