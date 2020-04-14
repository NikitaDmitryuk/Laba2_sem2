import re

# n - число уравнений (строк матрицы)
# C - диагональ, лежащая над главной (нумеруется: [0;n-2])
# B - главная диагональ матрицы A (нумеруется: [0;n-1])
# A - диагональ, лежащая под главной (нумеруется: [1;n-1])
# F - правая часть (столбец)
# x - решение, массив x будет содержать ответ


def tdma(A, C, B, F):

    alpha = [0]
    beta = [0]
    n = len(F)
    x = [0]*n

    for i in range(n-1):
        alpha.append(-C[i] / (A[i] * alpha[i] + B[i]))
        beta.append((F[i] - A[i] * beta[i]) / (A[i] * alpha[i] + B[i]))

    x[n-1] = (F[n - 1] - A[n - 2] * beta[n - 1]) / (B[n - 1] + A[n - 2] * alpha[n - 1])

    for i in reversed(range(n-1)):
        x[i] = alpha[i+1]*x[i+1] + beta[i+1]

    return x


def sweep_method(A, B, C, F, mu1, muN, N):
    u = [0 for _ in range(N + 1)]
    a = [0 for _ in range(N + 1)]
    b = [0 for _ in range(N + 1)]
    u[0] = mu1
    u[N] = muN
    b[1] = mu1

    for j in range(1, N):
        a[j + 1] = B / (C - a[j] * A)
        b[j + 1] = (A * b[j] + F[j]) / (C - a[j] * A)

    for j in range(N, 1, -1):
        u[j - 1] = a[j] * u[j] + b[j]

    return u


class DiffEquation:

    def __init__(self, type_of_method):
        self.u = []
        self.dx = None
        self.dt = None
        self.N = None
        self.K = None
        self.type_of_method = type_of_method

    def dsolve(self, x0, x1, t0, t1, dx, dt, border_conditions, a=1, betta=0, D=lambda x, t: 1,
               source_function=lambda x, t: 0):

        N = round((x1 - x0) / dx)
        K = round((t1 - t0) / dt)
        self.dx = dx
        self.dt = dt
        self.N = N
        self.K = K
        self.u = []
        # border_conditions = re.sub(r"[ ()=,u]", "", border_conditions)

        for k in range(K + 1):
            self.u.append([])
            for j in range(N + 1):
                self.u[k].append(0)

        for conditions in border_conditions:
            if conditions[0] == 'u(0, t)':
                for k in range(K + 1):
                    self.u[k][0] = conditions[1](0, k * dt)
            elif conditions[0] == 'u(l, t)':
                for k in range(K + 1):
                    self.u[k][N] = conditions[1](N * dx, k * dt)
            elif conditions[0] == 'u(x, 0)':
                self.u[0] = [conditions[1](j * dx, 0) for j in range(N + 1)]
            elif conditions[0].replace(' ', '') == 'b0,c0':
                b0, c0 = conditions[1]
            elif conditions[0].replace(' ', '') == 'bk,ck':
                bk, ck = conditions[1]

        if self.type_of_method == 'explicit method':
            self.explicit_method(a, source_function)
        elif self.type_of_method == 'implicit method':
            self.implicit_method(a, source_function)
        elif self.type_of_method == 'diffusion equation':
            self.diffusion_equation(D, betta, b0, c0, bk, ck)

        return self.u

    def explicit_method(self, a, source_function):
        dx = self.dx
        dt = self.dt
        N = self.N
        K = self.K

        for k in range(K):
            for j in range(1, N):
                u = self.u
                self.u[k + 1][j] = u[k][j] + a * a * dt / dx / dx * (u[k][j + 1] - 2 * u[k][j] + u[k][j - 1]) + \
                                   dt * source_function(j * dx, k * dt)

    def implicit_method(self, a, source_function):
        dx = self.dx
        dt = self.dt
        N = self.N
        K = self.K
        gamma = a * a * dt / (dx * dx)
        for k in range(K):
            F = [self.u[k][j] / gamma + a * a * dx * dx * source_function(j * dx, (k + 1) * dt) for j in range(N + 1)]
            self.u[k + 1] = sweep_method(1, 1, 2 + 1 / gamma, F, self.u[0][0], self.u[0][N], N)

    def diffusion_equation(self, D, betta, b0, c0, bk, ck):
        dx = self.dx
        dt = self.dt
        N = self.N
        K = self.K

        gamma = dt / dx / dx
        a = lambda j, m, n: D(j*dx + m * dx / 2, dt * (n + 1))  # m = +1, -1;

        A = [0.0]*(N + 1)
        B = [0.0]*(N + 1)
        C = [0.0]*N
        F = [0.0]*(N + 1)

        for k in range(K):
            B[0] = dx * betta / 2 - 1 / dx - dx / 2 / dt - b0((k + 1)*dt)
            C[0] = 1 / dx
            F[0] = c0((k + 1)*dt) - dx * self.u[k][0] / 2 / dt

            A[N] = - 1 / dx
            B[N] = 1 / dx + dx / 2 / dt - dx * betta / 2 - bk((k + 1)*dt)
            F[N] = ck((k + 1)*dt) + dx * self.u[k][N] / 2 / dt

            for j in range(1, N):
                A[j] = gamma * a(j, -1, k)
                C[j] = 1 + gamma * (a(j, 1, k) + a(j, -1, k)) - dt * betta
                B[j] = - gamma * a(j, 1, k)
                F[j] = self.u[k][j]

            self.u[k+1] = tdma(A=A, B=B, C=C, F=F)
