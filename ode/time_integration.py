
import numpy as np

class TimeIntegration():

    def __init__(self, t0: float, tf: float, dt: float,
        dt_min: float, dt_max: float, atol: float, rtol: float,
        u0: NDArray[np.float64]):

        self.t = t0 + dt

        # self.t0 = t0
        self.tf = tf

        self.dt = dt
        self.dt_min = dt_min
        self.dt_max = dt_max
        self.atol = atol
        self.rtol = rtol

        self.dtp = dt

        # self.u0 = u0

        self.up = u0
        self.up2 = u0

        self.inc = 1

    def update(self, dt: float, u: NDArray[np.float64]):

        self.inc += 1

        self.t += dt

        self.dtp = self.dt
        self.dt = dt

        self.up2 = self.up
        self.up = u

    def print_info(self):

        print(f"i = {self.inc}: t = {self.t:.4f}, dt = {self.dt:.4f}, dtp = {self.dtp:.4f}")

    def update_time_step(self, u: NDArray[np.float64]) -> float:

        dt = self.dt
        t = self.t
        tf = self.tf

        scale = self.atol + self.rtol*np.maximum(np.abs(self.up), np.abs(u))
        error_diff = np.abs(u - self.up)

        rms_err = np.sqrt(np.mean((error_diff/scale)**2))

        exponent = 1.0/3.0 # 2nd order
        dt_new = 0.9*dt*((1.0/(rms_err + 1e-20))**exponent)

        dt_new = min(dt_new, 1.2*dt)
        dt_new = max(dt_new, 0.5*dt)

        dt = min(max(self.dt_min, dt_new), self.dt_max)

        if t+dt >= tf:
            dt = tf - t

        # print(f"Old dt: {self.dt:.4e} | New dt: {dt:.4e} | Normalized Error: {rms_err:.4e}")

        return dt
