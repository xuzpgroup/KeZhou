import numpy as np

# Total number of fes files in folder
# total_files=100
# Min and max initial guesses
# min_min=50
# min_max=90
# max_min=90
# max_max=130

file_name="gas.z.dat"
# file_name="tt.dat"
matrix=np.genfromtxt(file_name)
# print(matrix)
# sum=matrix.sum(axis=0) ## coording to row (lie)
# print(sum)
ave=matrix.mean(axis=0)
print(ave)
# for i in range(total_files):
#         file_name="fes_" + str(i) + ".dat"
#         matrix=np.genfromtxt(file_name)
#         minimum=np.amin(matrix[min_min:min_max,1])
#         maximum=np.amax(matrix[max_min:max_max,1])
#         print(str(i) + " " + str(minimum) + " " + str(maximum) + " " + str(maximum-minimum))
