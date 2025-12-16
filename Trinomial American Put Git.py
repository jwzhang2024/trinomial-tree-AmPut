from math import exp, sqrt
def calculate(x, isRound=False, dp=4):
    return round(x,dp) if isRound else x
def TriTree_Am_Put(S0, K, r, sigma, T, N, isRound=False, dp=4):
   dt=calculate(T/N, isRound, dp)
   u=calculate(exp(sigma*sqrt(3*dt)), isRound, dp)
   d=calculate(1/u, isRound, dp)
   pu=calculate(1/6+(r-sigma**2/2)*sqrt(dt/(12*sigma**2)), isRound, dp)
   pm=calculate(2/3, isRound, dp)
   pd=calculate(1/6-(r-sigma**2/2)*sqrt(dt/(12*sigma**2)), isRound, dp)
   f=[[0 for j in range(-i,i+1)] for i in range(0,N+1)]
   for j in range(-N,N+1):
      f[N][j+N]=calculate(max(K-S0*u**j,0), isRound, dp)
   for i in range(N-1, 0-1, -1):
      for j in range(-i, i+1):
         f[i][j+i]=calculate(max(K-S0*u**j, exp(-r*dt)*(pu*f[i+1][j+1+i+1]
                                        +pm*f[i+1][j+i+1]
                                        +pd*f[i+1][j-1+i+1])), isRound, dp)
   return f[0][0]

S0,K,r,sigma,T,N=50, 50, 0.1, 0.4, 1, 3
print(TriTree_Am_Put(S0,K,r,sigma,T,N))
print(TriTree_Am_Put(S0,K,r,sigma,T,N,isRound=True))

print(TriTree_Am_Put(S0,K,r,sigma,T,N=12)) #change N to increase simulation number
