
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from diff_miron import DiffEquationMiron


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
    plt.savefig(title + ' miron.png', bbox_inches='tight')


dx = 0.1
dt = 0.006

x0 = 0
x1 = np.pi / 2

t0 = 0
t1 = 3


def main():
    # 'explicit method', 'implicit method', 'diffusion equation'

    diff_equation = DiffEquationMiron(type_of_method='explicit method')
    u = diff_equation.dsolve(x0=x0, x1=x1,
                             t0=t0, t1=t1,
                             dx=dx, dt=dt,
                             border_conditions=[('u(x, 0)', lambda x, t: 2 * x),
                                                ('u(l, t)', lambda x, t: np.pi)],
                             source_function=lambda x, t: np.exp(t) * np.cos(x))
    plot_tmp(u, 'Решение с помощью явной схемы')

    # diff_equation = DiffEquationMiron(type_of_method='implicit method')  # 'explicit method', 'implicit method'
    # u = diff_equation.dsolve(x0=x0, x1=x1,
    #                          t0=t0, t1=t1,
    #                          dx=dx, dt=dt,
    #                          border_conditions=[('u(x, 0)', lambda x, t: 2*x),
    #                                             ('u(l, t)', lambda x, t: np.pi)],
    #                          source_function=lambda x, t: np.exp(t)*np.cos(x))
    # plot_tmp(u, 'Решение с помощью неявной схемы')

    plt.show()


if __name__ == '__main__':
    main()
