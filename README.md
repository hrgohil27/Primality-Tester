CONTACT: HRGOHILMD@GMAIL.COM
# HNP 2.2.5 - Blaze - Harshkumar Gohil & Grok (xAI) - March 19, 2025 - MIT: Free use, no warranties - def hnp(n,m=0):return n>1and n%2>0and(pow(2,d:=n-1>>((s:=0)or bin(d)[2:].count('0')),n)in(1,n-1)or any((x:=x*x%n)==n-1for x in[pow(2,d,n)]))and((n**.2if m else n**.1)<1-(1<<n)%3/3if m else 1)

if __name__=="__main__":
 from time import time
 for t,m in[(2**2047+9,0),(31,1),(3321928094691,1)]:s=time();r=hnp(t,m);print(f"{2**t-1 if m else t}: {r}, {(time()-s)*1000:.2f} ms")
