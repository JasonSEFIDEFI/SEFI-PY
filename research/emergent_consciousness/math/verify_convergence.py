from math import exp, log2

def rk4(x, k, dt):
    def f(y):
        return (k*(y[1]-y[0]), k*(y[0]-y[1]))

    def add(a, b, s):
        return tuple(u+s*v for u,v in zip(a,b))

    a = f(x)
    b = f(add(x,a,dt/2))
    c = f(add(x,b,dt/2))
    d = f(add(x,c,dt))

    return tuple(
        x[i]+dt*(a[i]+2*b[i]+2*c[i]+d[i])/6
        for i in range(2)
    )

def error(n):
    dt = 1/n
    x = (1.0,-1.0)
    for _ in range(n):
        x = rk4(x,1.0,dt)
    exact = exp(-2)
    return max(abs(x[0]-exact),abs(x[1]+exact))

errors = [error(n) for n in (10,20,40,80)]
orders = [
    log2(errors[i]/errors[i+1])
    for i in range(3)
]

print("Checkpoint 011: RK4 convergence")
for n,e in zip((10,20,40,80),errors):
    print(f"steps={n:3d}, error={e:.12e}")

for order in orders:
    print(f"observed order={order:.6f}")
    assert 3.8 < order < 4.2

print("PASS: fourth-order convergence verified")
