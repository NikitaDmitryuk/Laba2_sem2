import re


class DiffEquation:

    def __init__(self, type_of_circuit='explicit'):
        self.u = []
        if type_of_circuit == 'explicit':
            self.def_circuit = self.explicit_circuit

    def dsolve(self, x0, x1, t0, t1, dx, dt, border_conditions, source_function=lambda x, t: 0):
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

        self.def_circuit(N, K, dx, dt, source_function)

        return self.u

    def explicit_circuit(self, N, K, dx, dt, source_function):
        for k in range(K):
            for j in range(1, N):
                u = self.u
                self.u[k + 1][j] = u[k][j] + dt / dx / dx * (u[k][j + 1] - 2 * u[k][j] + u[k][j - 1]) + \
                                   dt * source_function(j * dx, k * dt)
