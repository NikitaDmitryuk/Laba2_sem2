import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from diffequation import DiffEquation
from scipy import integrate


def write_file(U, name_file):
    with open(f'u_{name_file}.txt', 'w') as file:
        for u in U:
            file.write(' '.join(map(str, u))+'\n')


def read_file(name_file):
    data = []
    with open(f'u_{name_file}.txt', 'r') as f:
        for line in f:
            data.append([float(x) for x in line.split()])
        return data


def calc_potok(u):
    h_det1 = 0.01 * H1
    h_det2 = 0.1 * H1
    h_det3 = H1
    h_det4 = 10 * H1
    h_det5 = 20 * H1

    f = lambda t: k(t) * u[int(t / dt)][-1]
    F = integrate.quad(f, 0, t1)[0]

    f1 = lambda x, y: 1 / (2 * np.pi * ((H1 / 2 - x) ** 2 + (H2 / 2 - y) ** 2 + h_det1 ** 2))
    f2 = lambda x, y: 1 / (2 * np.pi * ((H1 / 2 - x) ** 2 + (H2 / 2 - y) ** 2 + h_det2 ** 2))
    f3 = lambda x, y: 1 / (2 * np.pi * ((H1 / 2 - x) ** 2 + (H2 / 2 - y) ** 2 + h_det3 ** 2))
    f4 = lambda x, y: 1 / (2 * np.pi * ((H1 / 2 - x) ** 2 + (H2 / 2 - y) ** 2 + h_det4 ** 2))
    f5 = lambda x, y: 1 / (2 * np.pi * ((H1 / 2 - x) ** 2 + (H2 / 2 - y) ** 2 + h_det5 ** 2))

    result1 = integrate.nquad(f1, [[0, H1], [0, H2]])
    result2 = integrate.nquad(f2, [[0, H1], [0, H2]])
    result3 = integrate.nquad(f3, [[0, H1], [0, H2]])
    result4 = integrate.nquad(f4, [[0, H1], [0, H2]])
    result5 = integrate.nquad(f5, [[0, H1], [0, H2]])

    print('T =', t1 // 3600)
    print('Поток в х1:', result1[0] * F)
    print('Поток в х2:', result2[0] * F)
    print('Поток в х3:', result3[0] * F)
    print('Поток в х4:', result4[0] * F)
    print('Поток в х5:', result5[0] * F)


def plot_tmp(u, title):
    t = (np.arange(len(u)) * dt + t0) / 3600
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
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$t$')
    ax.set_zlabel(r'$T$')
    plt.title(title)
    plt.savefig(title + ' nekit.png', bbox_inches='tight')


# dx = 0.01
# dt = 0.1
#
# x0 = 0
# x1 = 2*np.pi
#
# t0 = 0
# t1 = 2


dx = 0.01
dt = 10000


x0 = 0
x1 = 100

t0 = 0
t1 = 3600000

D0 = 0.01
D = lambda x, t: D0 * (1 + x / x1)
k = lambda t: 0.1

betta = -1e-6

H1 = 1e5
H2 = 2 * H1


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
    # plot_tmp(u, 'Решение с помощью явной схемы с нарушенным условием устойчивости')

    # diff_equation = DiffEquation(type_of_method='implicit method')  # 'explicit method', 'implicit method'
    # u = diff_equation.dsolve(x0=x0, x1=x1,
    #                          t0=t0, t1=t1,
    #                          dx=dx, dt=dt,
    #                          border_conditions=[('u(x, 0)', lambda x, t: 1),
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
                             D=D,
                             border_conditions=[('u(x, 0)', lambda x, t: 1),

                                                ('a0', lambda t: 0),
                                                ('an', lambda t: -k(t) / D(x1, t))]
                             )
    plot_tmp(u, r'Решение уравнения диффузии. $T={}$ часов.'.format(t1 // 3600))
    write_file(u, f't={t1//3600}')
    # u = read_file(f't={t1//3600}')

    calc_potok(u)

    plt.show()


if __name__ == '__main__':
    main()
