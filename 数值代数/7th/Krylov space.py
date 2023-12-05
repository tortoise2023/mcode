import math

# 计算有限维向量的2-范数
def norm2(n):
    tmp = 0
    for e in n:
        tmp += e ** 2
    return math.sqrt(tmp)

# 计算同维数向量的内积
def mul(n1, n2):
    res = 0
    dim = len(n1)
    for i in range(dim):
        res += n1[i] * n2[i]
    return res

# 验证矩阵是否为方阵
def isMartixN(A):
    n = len(A)
    for line in A:
        if len(line) != n:
            return False
    return True

def MGS(A,r,m):#Arnoldi过程（MGS）
    n=len(A);
    if n!=len(r):
        print("输入的向量与矩阵维数不匹配")
        return 0
    else:
        if not isMartixN(A):
            raise Exception('YULE ERROR 输入的矩阵非方阵 | CODE: 0x0001')
    v=[]
    # 计算v1
    for i in range(n):
        v.append(r[i]/norm2(r))
    # 把v1作为v的第1个分量
    v=[v]
    h=[]
    # j = 1, 2, ..., m - 1
    for j in range(m-1):
        #计算向量z,书上第三行
        z=[]
        for p in range(n):
            #计算z的第p个元素
            z.append(mul(A[p], v[j]))
        h.append([])
        for i in range(j+1):
            #计算h的i,j元
            h[j].append(mul(v[i], z))
            #更新z
            for p in range(n):
                z[p]-=h[j][i]*v[i][p]
        #计算h(j+1,j)
        h[j].append(norm2(z))
        # Python中计算结果为0时，会产生一个无限接近为0的数
        if h[j][j+1] < 1e-4:
            break
        #计算v_{j+1}
        v.append([])
        for p in range(n):
            v[j+1].append(z[p]/h[j][j+1])
    return v

def main1():
    A = [[1,0,0],[0,1,1],[0,0,1]]
    try:
        a=MGS(A,[1,1,1],2)
    except:
        pass
    print(a)
#如需执行Arnoldi过程（MGS）示例，则解除下行注释  
#main1() 

#生成个n行m列的维零矩阵
def zeros(n,m):
    if n==1 and isinstance(m,int):
        a=[]
        for i in range(m):
            a.append(0)
        return a
    elif n>1 and isinstance(n,int) and isinstance(m,int):
        a=[]
        for i in range(n):
            a.append([])
            for j in range(m):
                a[i].append(0)
        return a
    else:
        raise Exception('输入系数不符合要求')

#向量数乘
def mulc(L,R):
    if not isinstance(L,list):
        raise Exception("输入向量不符合要求")
    else:
        for i in L:
            if not isinstance(i,int) and not isinstance(i,float):
                raise Exception("输入向量不符合要求")
    if not isinstance(R,int) and not isinstance(R,float):
        raise Exception("输入数字不符合要求")
    n=len(L)
    E=zeros(1,n)
    for i in range(n):
        E[i]=L[i]*R
    return E

#矩阵右乘向量
def mulr(A,r):
    #判断A是否是list的双层嵌套结构
    if not isinstance(A,list):
        raise Exception("输入的并非矩阵001")
    else:
        n=len(A)
        for i in range(n):
            if not isinstance(A[i],list):
                raise Exception("非法的矩阵形式002")
    #判断A是否为长方形
    m=len(A[0])
    for i in range(n)[1:n]:
        if len(A[i])!=m:
            raise Exception("非法的矩阵形式003")
    #判断是否A为数字矩阵
    for i in range(n):
        for j in range(m):
            if not isinstance(A[i][j],int) and not isinstance(A[i][j],float):
                raise Exception("非法的矩阵形式004")
    #判断r是否合法
    if not isinstance(r,list):
        raise Exception("非法的向量形式005")
    else:
        k=len(r)
        if k!=m:
            raise Exception("向量与矩阵维数不匹配006")
        for i in range(k):
            if not isinstance(r[i],int) and not isinstance(r[i],float):
                raise Exception("非法的向量形式007")
    #计算A*r
    b=zeros(1,n)
    for i in range(n):
        b[i]=mul(A[i],r)
    return b

#向量加法
def suml(L1,L2):
    if len(L1)!=len(L2):
        raise Exception("两个向量长度不同")
    ans=[]
    for i in range(len(L1)):
        ans.append(L1[i]+L2[i])
    return ans

#Lanczos过程
def Lan(A,r,m):
    n=len(A);
    if n!=len(r):
        raise Exception("输入的向量与矩阵维数不匹配")
        return 0
    else:
        if not isMartixN(A):
            raise Exception(' 输入的矩阵非方阵 ')
    v=[]
    beta=zeros(1,m)
    #判断是否A为对称矩阵
    for i in range(n)[1:n]:
        for j in range(i):
            if A[i][j]!=A[j][i]:
                raise Exception(' 输入的矩阵非对称矩阵 ')
    #v_1,v_2,的生成
    v.append(zeros(1,n))
    v.append(mulc(r,1/norm2(r)))
    #j=1 to m-1
    for j in range(m-1):
        z=mulr(A,v[j+1])
        alpha=mul(v[j+1],z)
        
        I=mulc(v[j+1],-alpha)
        J=mulc(v[j],-beta[j])
        z=suml(suml(z,I),J)
        
        beta[j+1]=norm2(z)
        if beta[j+1]<1e-8:
            v.remove(v[0])
            return v
        v.append(mulc(z,1/beta[j+1]))
        
    v.remove(v[0])
    return v

def main2():
    A=[[1,2,2,3],[2,2,2,4],[2,2,3,5],[3,4,5,1]]
    r=[1,2,1,1]
    v=Lan(A,r,4)
    print(v)
    #print(mul(v[1],v[3]))
#如需执行Lanczos过程示例，则解除下行注释  
#main2()

#实用GMRES方法

def GMRES(A,b):
    #矩阵大小
    n=len(A)
    #选取初值
    x=zeros(1,n)
    x[0]=1
    #停机标准
    E=1e-4
    #最大迭代步数
    IterMax=1000
    r=suml(b,mulc(mulr(A,x),-1))
    beta=norm2(r)
    if beta/norm2(b)<E :
        #停止计算，输出近似解x
        return x
    v=[]
    v.append(mulc(r,1/beta))
    e1=zeros(1,n)
    e1[0]=1
    xi=mulc(e1,beta)#记录q1
    Hsize=0
    if IterMax>n:
        Hsize=n
    else:
        hsize= IterMax
    H=zeros(Hsize+1,Hsize)
    c=[]
    s=[]
    for j in range(IterMax):
        w=mulr(A,v[j])
        #Arnoldi过程
        #i=1,2,...,j
        for i in range(j+1):
            H[i][j]=mul(v[i],w)
            w=suml(w,mulc(v[i],-H[i][j]))
        H[j+1][j]=norm2(w)
        if H[j+1][j]<1e-4:
            m=j
            break #迭代中断
        v.append(mulc(w,1/H[j+1][j]))
        for i in range(j):
            [H[i][j],H[i+1][j]]=mulr([[c[i],s[i]],[-s[i],c[i]]],[H[i][j],H[i+1][j]])
        #构造Givens变换Gi
        
        if abs(H[j][j])>abs(H[j+1][j]):
            tau=H[j+1][j]/H[j][j]
            c.append(1/math.sqrt(1+tau**2))
            s.append(c[j]*tau)
        else:
            tau=H[j][j]/H[j+1][j]
            s.append(1/math.sqrt(1+tau**2))
            c.append(s[j]*tau)
        #计算Gj*Hj+1，j（1：j，j）
        H[j][j]=c[j]*H[j][j]+s[j]*H[j+1][j]
        H[j+1][j]=0
        [xi[j],xi[j+1]]=mulr([[c[j],s[j]],[-s[j],c[j]]],[xi[j],0])
        relres=abs(xi[j+1])/beta #相对残量
        if relres<E:
            m=j
            v.remove(v[j+1])
            break
    m=j
    y=zeros(1,m+1)
    #xi(1:m)左除H(1:m,1:m) 
    for i in range(m+1):
        if abs(H[m-i][m-i])<1e-8:
            raise Exception("H的第"+str(m-1-i)+"个主元的值太小了，快看看H对不对")
        y[m-i]=xi[m-i]/H[m-i][m-i]
        for k in range(i-1):
            y[m-i]-=y[m-k]*H[m-i][m-k]/H[m-i][m-i]
    for i in range(m+1):
        x=suml(x,mulc(v[i],y[i]))
    return x
A=[[1,2,0],[0,6,1],[0,0,1]]
b=[2,1,1]
r=mulr(A,[1,0,0])
r=suml(b,mulc(r,-1))
x=mulr(A,GMRES(A,b))

print(x)

