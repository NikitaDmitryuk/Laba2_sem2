import re


def sweep_method(A, B, C, F, mu1, mu2, N):
    u = [0 for _ in range(N + 1)]
    a = [0 for _ in range(N + 1)]
    b = [0 for _ in range(N + 1)]
    u[0] = mu1
    u[N] = mu2
    b[1] = mu1
    a[1] = 0

    for j in range(N):
        a[j+1] = B / (C - a[j] * A)
        b[j+1] = (A * b[j] + F[j]) / (C - a[j] * A)

    for j in range(N, 1, -1):
        u[j-1] = a[j] * u[j] + b[j]

    return u


class DiffEquation:

    def __init__(self, type_of_method):
        self.u = []
        if type_of_method == 'explicit method':
            self.function_method = self.explicit_method
        elif type_of_method == 'implicit method':
            self.function_method = self.implicit_method

    def dsolve(self, x0, x1, t0, t1, dx, dt, border_conditions, a, source_function=lambda x, t: 0):
        N = round((x1 - x0) / dx)
        K = round((t1 - t0) / dt)
        self.u = []
        border_conditions = re.sub(r"[ ()=,u]", "", border_conditions).split(';')

        for k in range(K + 1):
            self.u.append([])
            for j in range(N + 1):
                self.u[k].append(0)

        for conditions in border_conditions:
            if conditions[0] == '0':
                for i in range(K+1):
                    self.u[i][0] = float(conditions[2])
            elif conditions[0] == 'l':
                for i in range(K+1):
                    self.u[i][N] = float(conditions[2])
            elif conditions[1] == '0':
                self.u[0] = [float(conditions[2]) for _ in range(N+1)]

        self.function_method(N, K, dx, dt, a, source_function)

        return self.u

    def explicit_method(self, N, K, dx, dt, a, source_function):
        for k in range(K):
            for j in range(1, N):
                u = self.u
                self.u[k + 1][j] = u[k][j] + a * a * dt / dx / dx * (u[k][j + 1] - 2 * u[k][j] + u[k][j - 1]) + \
                                   dt * source_function(j * dx, k * dt)

    def implicit_method(self, N, K, dx, dt, a, source_function):
        gamma = a * a * dt / (dx * dx)
        for k in range(K):
            F = [self.u[k][j] / gamma + a * a * dx * dx * source_function(j*dx, (k+1)*dt) for j in range(N + 1)]
            self.u[k+1] = sweep_method(1, 1, 2 + 1 / gamma, F, self.u[0][0], self.u[0][N], N)
