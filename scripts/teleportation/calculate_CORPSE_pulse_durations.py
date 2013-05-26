f=np.linspace(7.45,7.55,11)
deg=np.array([180.,60.,300.,420.])

matrix=np.zeros([5 , 12])

matrix[0,1:]=f[:]
matrix[1:,0]=deg[:]


for i in range(5)[1:5]:
    for j in range(12)[1:12]:
        matrix[i,j]=round(deg[i-1]/180*1/(2*f[j-1])*1e3,2)
    
print matrix 