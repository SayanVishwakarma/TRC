import numpy as np
import matplotlib.pyplot as plt

# --- Constants & vehicle parameters ---
g = 9.81           # gravity (m/s²)
rho0 = 1.225       # sea-level air density (kg/m³)
H = 8500.0         # scale height (m)
Cd = 0.3           # drag coefficient
A = 10.0           # cross-sectional area (m²)

# --- Phase functions ---
def coast_phase(vx0, vy0, y0, m=20000, dt=0.1, t_end=10.0):
    vx, vy, x, y, t = vx0, vy0, 0.0, y0, 0.0
    hist = {'x':[], 'y':[], 'vx':[], 'vy':[]}
    while t < t_end and y > 0:
        rho = rho0 * np.exp(-y/H)
        v   = np.hypot(vx, vy)
        Fd  = 0.5 * rho * Cd * A * v**2
        axd = -Fd*(vx/v)/m if v>0 else 0
        ayd = -Fd*(vy/v)/m if v>0 else 0
        vx += axd*dt
        vy += (-g + ayd)*dt
        x  += vx*dt+0.5*axd*dt*dt
        y  += vy*dt+0.5*ayd*dt*dt
        t  += dt
        hist['x'].append(x); hist['y'].append(y)
        hist['vx'].append(vx); hist['vy'].append(vy)
    return hist


def boostback_burn(vx0, vy0, x0, y0, thrust_acc=16.0, m=160000, dt=0.1):
    vx, vy, x, y, t = vx0, vy0, x0, y0, 0.0
    target_vx = -0.45 * vx0
    dir_sign  = -np.sign(vx0)
    hist = {'x':[], 'y':[], 'vx':[], 'vy':[]}
    while (vx > target_vx if vx0>0 else vx<target_vx) and y>0:
        rho = rho0 * np.exp(-y/H)
        v   = np.hypot(vx, vy)
        Fd  = 0.5 * rho * Cd * A * v**2
        axd = -Fd*(vx/v)/m if v>0 else 0
        ayd = -Fd*(vy/v)/m if v>0 else 0
        ax = axd + thrust_acc*dir_sign
        ay = -g + ayd
        vx += ax*dt
        vy += ay*dt
        x  += vx*dt+0.5*axd*dt*dt
        y  += vy*dt+0.5*ayd*dt*dt
        t  += dt
        hist['x'].append(x); hist['y'].append(y)
        hist['vx'].append(vx); hist['vy'].append(vy)
    return hist


def final_coast(vx0, vy0, x0, y0, return_alt, m=20000, dt=0.1):
    vx, vy, x, y, t = vx0, vy0, x0, y0, 0.0
    hist = {'x':[], 'y':[], 'vx':[], 'vy':[]}
    while y > return_alt:
        rho = rho0 * np.exp(-y/H)
        v   = np.hypot(vx, vy)
        Fd  = 0.5 * rho * Cd * A * v**2
        axd = -Fd*(vx/v)/m if v>0 else 0
        ayd = -Fd*(vy/v)/m if v>0 else 0
        vx += axd*dt
        vy += (-g + ayd)*dt
        x  += vx*dt+0.5*axd*dt*dt
        y  += vy*dt+0.5*ayd*dt*dt
        t  += dt
        hist['x'].append(x); hist['y'].append(y)
        hist['vx'].append(vx); hist['vy'].append(vy)
    return hist


def entry_burn(vx0, vy0, x0, y0, thrust_acc=30.0, m=30000, dt=0.1):
    vx, vy, x, y, t = vx0, vy0, x0, y0, 0.0
    hist = {'x':[], 'y':[], 'vx':[], 'vy':[]}
    while (abs(vx)>1 or abs(vy)>1) and y>0:
        rho = rho0 * np.exp(-y/H)
        v   = np.hypot(vx, vy)
        Fd  = 0.5 * rho * Cd * A * v**2
        axd = -Fd*(vx/v)/m if v>0 else 0
        ayd = -Fd*(vy/v)/m if v>0 else 0
        ux, uy = (vx/v, vy/v) if v>0 else (0,1)
        ax = axd - thrust_acc*ux
        ay = -g + ayd - thrust_acc*uy
        vx += ax*dt
        vy += ay*dt
        x  += vx*dt+0.5*axd*dt*dt
        y  += vy*dt+0.5*ayd*dt*dt
        t  += dt
        hist['x'].append(x); hist['y'].append(y)
        hist['vx'].append(vx); hist['vy'].append(vy)
    return hist

# Unpowered free-fall from entry end down to 2 km
def unpowered_descent(vx0, vy0, x0, y0, target_alt=2500.0, m=40000, dt=0.1):
    vx, vy, x, y, t = vx0, vy0, x0, y0, 0.0
    hist = {'x':[], 'y':[], 'vx':[], 'vy':[]}
    while y > target_alt:
        rho = rho0*np.exp(-y/H)
        v   = np.hypot(vx, vy)
        Fd  = 0.5 * rho * Cd * A * v**2
        axd = -Fd*(vx/v)/m if v>0 else 0
        ayd = -Fd*(vy/v)/m if v>0 else 0
        vx += axd*dt
        vy += (-g + ayd)*dt
        x  += vx*dt+0.5*axd*dt*dt
        y  += vy*dt+0.5*ayd*dt*dt
        t += dt
        hist['x'].append(x); hist['y'].append(y)
        hist['vx'].append(vx); hist['vy'].append(vy)
    return hist

# Landing burn from 2 km until vertical velocity = 0, continuing horizontal track
def landing_burn(vx0, vy0, x0, y0, thrust_acc=25.0, m=30000, dt=0.1):
    vx, vy, x, y, t = vx0, vy0, x0, y0, 0.0
    hist = {'x':[], 'y':[], 'vx':[], 'vy':[]}
    while vy < 0:
        rho = rho0*np.exp(-y/H)
        v   = np.hypot(vx, vy)
        Fd  = 0.5 * rho * Cd * A * v**2
        axd = -Fd*(vx/v)/m if v>0 else 0
        ayd = -Fd*(vy/v)/m if v>0 else 0
        # thrust only vertical
        ax = axd
        ay = -g + ayd + thrust_acc
        vx += ax*dt
        vy += ay*dt
        x  += vx*dt+0.5*axd*dt*dt
        y  += vy*dt+0.5*ayd*dt*dt
        t  += dt
        hist['x'].append(x); hist['y'].append(y)
        hist['vx'].append(vx); hist['vy'].append(vy)
    return hist

# --- Main execution ---
sep_vx, sep_vy, sep_alt = 1500.0, 1100.0, 90000.0

# 1) Coast
h1 = coast_phase(sep_vx, sep_vy, sep_alt)
vx1, vy1, y1 = h1['vx'][-1], h1['vy'][-1], h1['y'][-1]

# 2) Boostback
x1, y1 = h1['x'][-1], h1['y'][-1]
h2 = boostback_burn(vx1, vy1, x1, y1)
vx2, vy2, y2 = h2['vx'][-1], h2['vy'][-1], h2['y'][-1]

# 3) Return coast to 60 km
entry_start_alt = 60000.0
x2, y2 = h2['x'][-1], h2['y'][-1]
h3 = final_coast(vx2, vy2, x2, y2, return_alt=entry_start_alt)
vx3, vy3, y3 = h3['vx'][-1], h3['vy'][-1], h3['y'][-1]

# 4) Entry burn
x3, y3 = h3['x'][-1], h3['y'][-1]
h4 = entry_burn(vx3, vy3, x3, y3)
vx4, vy4, y4 = h4['vx'][-1], h4['vy'][-1], h4['y'][-1]

# 5) Unpowered free-fall to 2 km
x4, y4 = h4['x'][-1], h4['y'][-1]
h5 = unpowered_descent(vx4, vy4, x4, y4)
vx5, vy5, x5, y5 = h5['vx'][-1], h5['vy'][-1], h5['x'][-1], h5['y'][-1]
print(f"Start Landing  → vx={vx5:.1f}, vy={vy5:.1f}, alt={y5:.1f} m")

# 6) Landing burn from end of unpowered descent
h6 = landing_burn(vx5, vy5, x5, y5)
vx6, vy6, x6, y6 = h6['vx'][-1], h6['vy'][-1], h6['x'][-1], h6['y'][-1]
print(f"Touchdown      → vx={vx6:.1f}, vy={vy6:.1f}, alt={y6:.1f} m")
h6 = landing_burn(vx5, vy5, x5, y5)
vx6, vy6, x6, y6 = h6['vx'][-1], h6['vy'][-1], h6['x'][-1], h6['y'][-1]

# 7) Print states
print(f"End Coast      → vx={vx1:.1f}, vy={vy1:.1f}, alt={y1:.1f} m")
print(f"End Boostback  → vx={vx2:.1f}, vy={vy2:.1f}, alt={y2:.1f} m")
print(f"Start Entry    → vx={vx3:.1f}, vy={vy3:.1f}, alt={y3:.1f} m")
print(f"End Entry      → vx={vx4:.1f}, vy={vy4:.1f}, alt={y4:.1f} m")
print(f"Start Landing  → vx={vx5:.1f}, vy={vy5:.1f}, alt={y5:.1f} m")
print(f"Touchdown      → vx={vx6:.1f}, vy={vy6:.1f}, alt={y6:.1f} m")

# 8) Plot full profile
X = np.hstack([h1['x'], h2['x'], h3['x'], h4['x'], h5['x'], h6['x']])
Y = np.hstack([h1['y'], h2['y'], h3['y'], h4['y'], h5['y'], h6['y']])
p1 = len(h1['x'])
p2 = p1 + len(h2['x'])
p3 = p2 + len(h3['x'])
p4 = p3 + len(h4['x'])
p5 = p4 + len(h5['x'])

plt.figure(figsize=(8,5))
X=X/1000
Y=Y/1000
plt.plot(X[:p1],    Y[:p1],    label='Coast')
plt.plot(X[p1:p2],  Y[p1:p2],  label='Boostback')
plt.plot(X[p2:p3],  Y[p2:p3],  label='Return to 60 km')
plt.plot(X[p3:p4],  Y[p3:p4],  label='Entry Burn')
plt.plot(X[p4:p5],  Y[p4:p5],  label='Free-fall to 2 km')
plt.plot(X[p5:],    Y[p5:],    label='Landing Burn')
plt.xlabel('Downrange (km)')
plt.ylabel('Altitude (km)')
plt.title('Full RTLS Mission Profile')
plt.legend(); plt.grid(True); plt.tight_layout(); plt.show()