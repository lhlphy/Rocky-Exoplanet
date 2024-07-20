import numpy as np
import time

def decorator_timer(name):
    def run_time(func):
        def warp(*args, **kwargs):
            t1 = time.time()
            temp = func(*args, **kwargs)
            t2 = time.time()
            print(name , 'process time:', t2-t1,' s')
            return temp
        return warp
    return run_time

def f(x):
    print(x**x)
    return x**x

@decorator_timer('func1')
f(3)
