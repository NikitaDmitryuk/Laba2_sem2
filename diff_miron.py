import re

# n - число уравнений (строк матрицы)
# C - диагональ, лежащая над главной (нумеруется: [0;n-2])
# B - главная диагональ матрицы A (нумеруется: [0;n-1])
# A - диагональ, лежащая под главной (нумеруется: [1;n-1])
# F - правая часть (столбец)
# x - решение, массив x будет содержать ответ


def tdma_miron(A, C, B, F, un):

    alpha = [0]
    beta = [0]
    n = len(F)
    x = [0]*n

    for i in range(n-1):
        alpha.append(-C[i] / (A[i] * alpha[i] + B[i]))
        beta.append((F[i] - A[i] * beta[i]) / (A[i] * alpha[i] + B[i]))

    x[n-1] = un

    for i in reversed(range(n-1)):
        x[i] = alpha[i+1]*x[i+1] + beta[i+1]

    return x


class DiffEquationMiron:

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

        if self.type_of_method == 'explicit method':
            self.explicit_method_miron(source_function)
        elif self.type_of_method == 'implicit method':
            self.implicit_method_miron(source_function)

        return self.u

    def explicit_method_miron(self, source_function):
        dx = self.dx
        dt = self.dt
        N = self.N
        K = self.K

        for k in range(K):
            self.u[k + 1][0] = self.u[k][0] + 2 * dt / dx * (dx / 2 * source_function(0, k*dt) + (self.u[k][1] - self.u[k][0]) / dx)
            for j in range(1, N):
                u = self.u
                self.u[k + 1][j] = u[k][j] + dt / dx / dx * (u[k][j + 1] - 2 * u[k][j] + u[k][j - 1]) + \
                                   dt * source_function(j * dx, k * dt)

    def implicit_method_miron(self, source_function):
        dx = self.dx
        dt = self.dt
        N = self.N
        K = self.K
        gamma = dt / (dx * dx)
        A = [1] * (N + 1)
        B = [-2 - 1 / gamma] * (N + 1)
        C = [1] * (N + 1)

        B[0] = - 1 / dx - dx / 2 / dt
        C[0] = 1 / dx

        for k in range(K):
            F = [-self.u[k][j] / gamma - dx * dx * source_function(j * dx, (k + 1) * dt) for j in range(N + 1)]
            F[0] = -dx / 2 / dt * self.u[k][0] - dx / 2 * source_function(0, k*dt)
            self.u[k + 1] = tdma_miron(A=A, B=B, C=C, F=F, un=self.u[k][N])
