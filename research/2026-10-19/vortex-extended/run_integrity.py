import math

def eligible(result):
    if result.get('reason')!='completed':return False
    try:
        achieved=float(result['samples'][-1]['time']);requested=float(result['T_requested'])
    except (KeyError,IndexError,TypeError,ValueError):return False
    return math.isfinite(achieved) and math.isfinite(requested) and requested>0 and abs(achieved-requested)<max(1e-8,1e-8*requested)
