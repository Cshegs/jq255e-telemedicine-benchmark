import csv, hashlib, math, os, random, secrets, statistics, time
from scipy.stats import shapiro, ttest_rel, wilcoxon
from cryptography.hazmat.primitives.ciphers.aead import AESGCM
from cryptography.hazmat.primitives.asymmetric import x25519, ec
Q=2**255-18651; R=2**254-131528281291764213006042413802501683931; G=(3,1); O=(1,0)
TRIALS=50; WARMUP=5; SIZES=[13,55,1000]; curves=['jq255e','X25519','P-256']
def add(P,Qp):
 if P==O:return Qp
 if Qp==O:return P
 E1,u1=P;E2,u2=Qp; a=u1*u1%Q;b=u2*u2%Q;t=8*a%Q*b%Q; D=(1-t)%Q; inv=pow(D,Q-2,Q);u3=(u1*E2+u2*E1)%Q*inv%Q; e3=((1+t)%Q*E1%Q*E2+16*u1%Q*u2%Q*((a+b)%Q))%Q*inv%Q*inv%Q;return e3,u3
def mul(k,P):
 R0=O;Qp=P
 while k:
  if k&1:R0=add(R0,Qp)
  Qp=add(Qp,Qp);k>>=1
 return R0
H=mul(2,G)
def jq_kg():
 sk=secrets.randbelow(R-1)+1;return sk,mul(sk,H)
def jq_shared(sk,pub):return mul(sk,pub)
def jq_enc(msg,pub):
 k=secrets.randbelow(R-1)+1;C=mul(k,H);S=mul(k,pub);key=hashlib.blake2s(S[0].to_bytes(32,'little')+S[1].to_bytes(32,'little')).digest();n=os.urandom(12);ct=AESGCM(key).encrypt(n,msg,None);return C,n,ct
def jq_dec(sk,p):
 C,n,ct=p;S=mul(sk,C);key=hashlib.blake2s(S[0].to_bytes(32,'little')+S[1].to_bytes(32,'little')).digest();return AESGCM(key).decrypt(n,ct,None)
def x_kg():sk=x25519.X25519PrivateKey.generate();return sk,sk.public_key()
def x_shared(sk,pub):return sk.exchange(pub)
def x_enc(msg,pub):
 k=x25519.X25519PrivateKey.generate();S=k.exchange(pub);n=os.urandom(12);ct=AESGCM(hashlib.blake2s(S).digest()).encrypt(n,msg,None);return k.public_key(),n,ct
def x_dec(sk,p):C,n,ct=p;return AESGCM(hashlib.blake2s(sk.exchange(C)).digest()).decrypt(n,ct,None)
def p_kg():sk=ec.generate_private_key(ec.SECP256R1());return sk,sk.public_key()
def p_shared(sk,pub):return sk.exchange(ec.ECDH(),pub)
def p_enc(msg,pub):
 k=ec.generate_private_key(ec.SECP256R1());S=k.exchange(ec.ECDH(),pub);n=os.urandom(12);ct=AESGCM(hashlib.blake2s(S).digest()).encrypt(n,msg,None);return k.public_key(),n,ct
def p_dec(sk,p):C,n,ct=p;return AESGCM(hashlib.blake2s(sk.exchange(ec.ECDH(),C)).digest()).decrypt(n,ct,None)
KG={'jq255e':jq_kg,'X25519':x_kg,'P-256':p_kg}; SH={'jq255e':jq_shared,'X25519':x_shared,'P-256':p_shared}; ENC={'jq255e':jq_enc,'X25519':x_enc,'P-256':p_enc}; DEC={'jq255e':jq_dec,'X25519':x_dec,'P-256':p_dec}
def tcall(fn):
 t=time.perf_counter_ns();z=fn();return (time.perf_counter_ns()-t)/1000,z
# warmups
for c in curves:
 for _ in range(WARMUP):KG[c]()
# keygen raw samples (not paired)
keygen={c:[] for c in curves}
for c in curves:
 for _ in range(TRIALS):keygen[c].append(tcall(KG[c])[0])
# paired operation benchmark
raw=[]
rng=random.Random(20261008)
for size in SIZES:
 msg=b'P'*size
 keys={c:KG[c]() for c in curves}
 payloads={c:ENC[c](msg,keys[c][1]) for c in curves}
 for op in ['ecdh','encryption','decryption']:
  vals={c:[] for c in curves}
  for i in range(TRIALS):
   order=curves.copy();rng.shuffle(order)
   for c in order:
    sk,pub=keys[c]
    if op=='ecdh': fn=lambda c=c,sk=sk,pub=pub:SH[c](sk,pub)
    elif op=='encryption': fn=lambda c=c,pub=pub,msg=msg:ENC[c](msg,pub)
    else: fn=lambda c=c,sk=sk,p=payloads[c]:DEC[c](sk,p)
    vals[c].append(tcall(fn)[0])
  for c in curves:
   for j in range(TRIALS): raw.append([size,op,c,j,vals[c][j]])
  for c in ['X25519','P-256']:
   d=[vals['jq255e'][i]-vals[c][i] for i in range(TRIALS)]
   print(size,op,'jq-',c,'mean diff',statistics.mean(d),'sd',statistics.stdev(d),'t',ttest_rel(vals['jq255e'],vals[c]),'wilcoxon',wilcoxon(d,zero_method='wilcox',alternative='two-sided',method='auto'))
   print('shapiro diff',shapiro(d))
with open('paired_benchmark_raw.csv','w',newline='') as f:
 w=csv.writer(f);w.writerow(['message_size_bytes','operation','curve','trial','time_us']);w.writerows(raw)
with open('paired_benchmark_keygen.csv','w',newline='') as f:
 w=csv.writer(f);w.writerow(['curve','trial','time_us']);
 for c in curves:
  for i,x in enumerate(keygen[c]):w.writerow([c,i,x])
