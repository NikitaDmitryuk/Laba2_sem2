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


dx = 0.1
dt = 0.005

x0 = 0
x1 = 2 * np.pi

t0 = 0
t1 = 2

a = 1


def main():
    diff_equation = DiffEquation(type_of_method='explicit method')  # 'explicit method', 'implicit method'
    u = diff_equation.dsolve(x0=x0, x1=x1,
                             t0=t0, t1=t1,
                             dx=dx, dt=dt,
                             border_conditions='u(x, 0)=1; u(0, t)=0; u(l, t)=0',
                             a=a,
                             source_function=lambda x, t: np.sin(x))
    plot_tmp(u, 'Решение с помощью явной схемы')

    diff_equation = DiffEquation(type_of_method='implicit method')  # 'explicit method', 'implicit method'
    u = diff_equation.dsolve(x0=x0, x1=x1,
                             t0=t0, t1=t1,
                             dx=dx, dt=dt,
                             border_conditions='u(x, 0)=1; u(0, t)=0; u(l, t)=0',
                             a=a,
                             source_function=lambda x, t: np.sin(x))
    plot_tmp(u, 'Решение с помощью неявной схемы')

    plt.show()


if __name__ == '__main__':
    main()
