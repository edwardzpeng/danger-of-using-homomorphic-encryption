Debug=False
q=2^54
d=1024
t=83
T=100
l=floor(log(q,T).n())
P.<x>=PolynomialRing(ZZ)
f=x^d+1
R.<X>=P.quotient(f)
Zq=Zmod(q)
Zt=Zmod(t)
Pq.<x>=PolynomialRing(Zq)
Pt.<x>=PolynomialRing(Zt)
sigma = 1.0
delta=floor(q/t)
from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
D = DiscreteGaussianDistributionIntegerSampler(sigma=sigma)

def sample_e():
    return R([D() for _ in range(d)])
def sample_2():
    return R([randint(0,1) for _ in range(d)])
def sample_r():
    return R([randint(0,q-1) for _ in range(d)])

def rq(a):
    A=a.list()
    for i in range(len(A)):
        A[i]=A[i]%q
    return R(A)
def rt(a):
    A=a.list()
    for i in range(len(A)):
        A[i]=A[i]%t
    return R(A)

def Roundq(a):
    A=a.list()
    for i in range(len(A)):
        A[i]=A[i]%q
        if A[i]>(q/2):
            A[i]=A[i]-q
    return R(A)
def Roundt(a):
    A=a.list()
    for i in range(len(A)):
        A[i]=A[i]%t
        if A[i]>(t/2):
            A[i]=A[i]-t
    return R(A)



def gen():
    s=sample_2()
    a=Roundq(sample_r())
    e=Roundq(sample_e())
    pk=[Roundq(-(a*s+e)),a]
    return s,pk


def encrypt(m):
    u=sample_2()
    e1=sample_e()
    e2=sample_e()
    return (Roundq(pk[0]*u+e1+delta*m),Roundq(pk[1]*u+e2))

def decrypt(c):
    tmp=t*Roundq(c[0]+c[1]*s)
    TMP=tmp.list()
    for i in range(len(TMP)):
        TMP[i]=round(TMP[i]/q)
    tmp=R(TMP)
    return Roundt(tmp)

def add(c1,c2):
    return (Roundq(c1[0]+c2[0]),Roundq(c1[1]+c2[1]))

def mul(c1,c2):
    tmp=t*c1[0]*c2[0]
    TMP=tmp.list()
    for i in range(len(TMP)):
        TMP[i]=round(TMP[i]/q)
    tmp=R(TMP)    
    c0_p=Roundq(tmp)
    tmp=t*(c1[0]*c2[1]+c1[1]*c2[0])
    TMP=tmp.list()
    for i in range(len(TMP)):
        TMP[i]=round(TMP[i]/q)
    tmp=R(TMP)
    c1_p=Roundq(tmp)
    tmp=t*c1[1]*c2[1]
    TMP=tmp.list()
    for i in range(len(TMP)):
        TMP[i]=round(TMP[i]/q)
    tmp=R(TMP)
    c2_p=Roundq(tmp)
    return (c0_p,c1_p,c2_p)



def gen_rlk(T):
    rlk=[]
    E=[]
    for i in range(l+1):
        a=sample_r()
        e=sample_e()
        rlk.append((Roundq(-(a*s+e)+(T^i)*s*s),a))
        E.append(e)
    return rlk,E

def baseT(num,b=T):
    V=[]
    flag=1
    if(num<0):
        num=-num
        flag=-1
    while (num> 0):
        temp = num % b
        if temp > b/2: 
            temp = temp -b
        V.append(temp*flag)
        num = (num-temp)//b
    for i in range(l+1-len(V)):
        V.append(0)
    return V

def baseCT(cc):
    CC=cc.list()    
    M=[]
    for i in range(len(CC)):
        M.extend(baseT(CC[i]))
    M=matrix(ZZ,d,l+1,M)
    CC=[]
    for i in range(l+1):
        tmpM=M.column(i).list()
        CC.append(P(tmpM))
    return CC

def relin(c):
    c0=c[0]
    c1=c[1]
    CC=baseCT(c[2])
    for i in range(l+1):
        c0+=rlk[i][0]*CC[i]
        c1+=rlk[i][1]*CC[i]
    return (Roundq(c0),Roundq(c1))

def mul2(c1,c2):
    return relin(mul(c1,c2))
   

def IntegerEncoder(i,b=2):
    CC=baseT(i,b)
    return R(CC)

def IntegerDecoder(m,b=2):
    CC=m.list()
    tmp=0
    for i in range(len(CC)):
        tmp+=CC[i]*b^i
    return tmp


def baseT(num,b=T):
    V=[]
    flag=1
    if(num<0):
        num=-num
        flag=-1
    while (num> 0):
        temp = num % b
        if temp > b/2: 
            temp = temp -b
        V.append(temp*flag)
        num = (num-temp)//b
    for i in range(l+1-len(V)):
        V.append(0)
    return V
def recover_key():
    cc0=0
    cc1=delta
    print "Recover private key successfully:",(s).list()==decrypt([cc0,cc1]).list()
    print "Secret key:",Roundt(s).list()
    print "recovered key",decrypt([cc0,cc1]).list()

s,pk=gen()
rlk,E=gen_rlk(T)
recover_key()