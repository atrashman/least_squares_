include("Least_squares.jl")
using .Least_squares#导入模块
sample,target = read_txt("xixi1.txt")
B = Least_squares.least_squares(sample,target)
#查看预测结果
sample*B
