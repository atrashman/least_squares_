# 最小二乘法基于矩阵的计算：(记得在样本属性加多一列1)
# 公式推导：
# 令：AX = Y = B
# X为要求的系数  A为给定样本  B为目标值则在X不确定的情况下，Y可以看成一个平面，是由样本向量构为基构成的
# A的每一个属性为一个长的多维向量
# 若要AX最接近B，只需要点B距离Y平面的距离最小，即B与平面Y的一个连线正交
# 解法：
# 令C = B-AX  （C实际上就是那条连线）
# C垂直样本向量为基向量构成的平面Y
# 必有：C*每一个属性向量 = 0
# 这样等价于：A'(B-AX) = 0 或者 A'AX = A'B 或者 令D = A'A ,X = D^-1 * A'B
module Least_squares
    using DelimitedFiles
    export read_txt,inverse,least_squares
#读取文件
    function read_txt(filename)#则默认最后一列为目标值
        print("默认最后一列为目标值")

        matrix = readdlm(filename,'\t')
        len = size(matrix)[1]
        return ([matrix[:,1:end-1] ones(Float64,len,1)],matrix[:,end])
    end
#计算逆
    function inverse(matrix_)#自定义的矩阵求逆
        matrix = matrix_
        row,col = size(matrix)
        (row == col) || error("row is not equal to col")
        #行列式不为零
        #还要自己构造单位矩阵
        E = float([(i == j) ? 1 : 0 for i in 1:row,j in 1:row])
        for i = 1:col
            if matrix[i,i] == 0#i行i号位等于零的话把下面i号位不等于零的数加上来
                k = i + 1
                while matrix[k,i] == 0 && k < row#循环判断一种逆不存在的情况（加减后整行为零）
                    local k
                    k+=1
                end
                if matrix[k,i] == 0
                    error("no inversion")
                end
                matrix[i,:] += matrix[k,:]
                E[i,:] += E[k,:]
            end

            div_=matrix[i,i]
            matrix[i,:] = matrix[i,:]/div_
            E[i,:] = E[i,:] / div_

            for ii  = 1:row #初等行变换
                if ii != i
                    tmp = - matrix[ii,i]*matrix[i,:]
                    tmp_= - matrix[ii,i]*E[i,:]
                    matrix[ii,:] += tmp#[0,_,_,_,_...]
                    E[ii,:] += tmp_
                end
            end
        end
        return E
    end
#按公式计算最小二乘法：
    function least_squares(sample_matrix,target)#不用加1
        matrix = sample_matrix' * sample_matrix
        matrix_inversion = inverse(matrix)
        return matrix_inversion * sample_matrix' * target
        # 这样等价于：A'(B-AX) = 0 或者 A'AX = A'B 或者 令D = A'A ,X = D^-1 * A'B

    end
end
