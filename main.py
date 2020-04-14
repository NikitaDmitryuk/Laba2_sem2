import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from diffequation import DiffEquation


def plot_tmp(u, title):
    t = np.arange(len(u)) * dt + t0
    x = np.arange(len(u[0])) * dx + x0
    ls = []
    for i in range(len(u)):
        ls += u[i]

    z = np.asarray(ls)
    X, T = np.meshgrid(x, t)
    Z = z.reshape(X.shape)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, T, Z, cmap=cm.coolwarm)
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_zlabel('T')
    plt.title(title)
    plt.savefig(title + '.png', bbox_inches='tight')


# dx = 0.1
# dt = 0.005
#
# x0 = 0
# x1 = 2 * np.pi
#
# t0 = 0
# t1 = 10

dx = 0.01
dt = 100

x0 = 0
x1 = 100

t0 = 0
t1 = 6000

D0 = 0.01
D = lambda x, t: D0 * (1 + x / x1)
k = lambda t: 0.1
betta = -1e-6


def main():
    # 'explicit method', 'implicit method', 'diffusion equation'

    # diff_equation = DiffEquation(type_of_method='explicit method')
    # u = diff_equation.dsolve(x0=x0, x1=x1,
    #                          t0=t0, t1=t1,
    #                          dx=dx, dt=dt,
    #                          border_conditions=[('u(x, 0)', lambda x, t: 1),
    #                                             ('u(0, t)', lambda x, t: 0),
    #                                             ('u(l, t)', lambda x, t: 0)],
    #                          source_function=lambda x, t: np.sin(x))
    # plot_tmp(u, 'Решение с помощью явной схемы')
    #

    # diff_equation = DiffEquation(type_of_method='implicit method')  # 'explicit method', 'implicit method'
    # u = diff_equation.dsolve(x0=x0, x1=x1,
    #                          t0=t0, t1=t1,
    #                          dx=dx, dt=dt,
    #                          border_conditions=[('u(x, 0)', lambda x, t: np.sin(x/2)),
    #                                             ('u(0, t)', lambda x, t: 0),
    #                                             ('u(l, t)', lambda x, t: 0)],
    #                          source_function=lambda x, t: np.sin(t))
    # plot_tmp(u, 'Решение с помощью неявной схемы')

    # dU(0, t) / dx = b0 * U(0, t) + c0(t)
    # dU(l, t) / dx = bk * U(l, t) + ck(t)
    diff_equation = DiffEquation(type_of_method='diffusion equation')
    u = diff_equation.dsolve(x0=x0, x1=x1,
                             t0=t0, t1=t1,
                             dx=dx, dt=dt,
                             betta=betta,
                             D=lambda x, t: 1,
                             border_conditions=[('u(x, 0)', lambda x, t: 1),
                                                ('a0', lambda t: 0),
                                                ('an', lambda t: -k(t) / D(x1, t))]
                             )
    plot_tmp(u, 'Решение уравнения диффузии')

    plt.show()


if __name__ == '__main__':
    main()
