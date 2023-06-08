import numpy as np
from matplotlib import pyplot as plt
import scipy.stats as stats
import pandas as pd


def amount_of_transfer(giver, receiver, lamb=0.1, stochastic=True, capacity=1):
    new_lamb = lamb*(capacity-receiver)
    if stochastic:
        while True:
            amount = np.random.exponential(new_lamb)
            if amount <= capacity-receiver:
                break
        return min(giver, amount)
    return min(giver, capacity-receiver, new_lamb)


def check_transfer_rule():
    lis = []
    for ii in range(1000):
        lis.append(amount_of_transfer(1, 0))

    plt.hist(lis, density=True)
    plt.show()

    print(np.mean(lis))


def update(ants, total_given, lamb=0.15, stochastic=True, capacity=1):
    giver, receiver = np.random.randint(num_ants, size=2)
    amount = amount_of_transfer(ants[giver], ants[receiver], lamb, stochastic, capacity)
    ants[giver] -= amount
    ants[receiver] += amount
    total_given[giver] += amount
    return ants, total_given


def fit_gamma(data, floc=-0.01, ax=None):
    k, loc1, theta = stats.gamma.fit(data, floc=floc)
    if ax:
        x = np.linspace(0.01, max(data), 100)
        gamma = stats.gamma.pdf(x, k, loc1, theta)
        ax.plot(x, gamma, 'r')
        ax.set_title('total given | ' + f'k={k:.2f},'+r'$\theta$'+f'={theta:.2f}')
    return k, theta


def run(lamb, steps, capacity, num_ants, to_plot=True):
    ants = capacity*np.random.random(num_ants)/2
    total_given = np.zeros(num_ants)

    for ii in range(steps*num_ants):
        update(ants, total_given, lamb=lamb, capacity=capacity)

    if to_plot:
        fig, axs = plt.subplots(1, 2)
        axs[0].hist(total_given, density=True)
        k, theta = fit_gamma(total_given, ax=axs[0])
        axs[1].hist(ants, density=True)
        plt.show()
        return k, theta

    k, theta = fit_gamma(total_given, ax=None)
    return k, theta


def plot_results(results, parameter):
    df_k = pd.DataFrame(results['k'])
    df_theta = pd.DataFrame(results['theta'])
    k_mean = df_k.mean()
    k_std = df_k.std()
    theta_mean=df_theta.mean()
    theta_std=df_theta.std()

    plt.figure()
    plt.subplot(121)
    plt.plot(theta_mean)
    plt.title('theta')
    plt.xlabel(parameter)
    plt.subplot(122)
    plt.plot(k_mean)
    plt.title('k')
    plt.xlabel(parameter)


capacity = 1
steps = 10
num_ants = 100
results = {'k': {}, 'theta': {}}
for lamb in np.linspace(0.02, 1, 50):
    print(lamb)
    for param in ['k', 'theta']:
        results[param][lamb] = []
    for jj in range(50):
        # print(jj)

        k, theta = run(lamb, steps, capacity, num_ants, to_plot=False)  # (jj == 1))

        results['k'][lamb].append(k)
        results['theta'][lamb].append(theta)

plot_results(results, 'lambda')


results_steps = {'k': {}, 'theta': {}}
lamb = 0.15
for steps in np.linspace(1, 20, 20):
    steps=int(steps)
    print(steps)
    for param in ['k', 'theta']:
        results_steps[param][steps] = []
    for jj in range(50):
        # print(jj)

        k, theta = run(lamb, steps, capacity, num_ants, to_plot=False)  # (jj == 1))

        results_steps['k'][steps].append(k)
        results_steps['theta'][steps].append(theta)


plot_results(results_steps, '# steps per ant')

plt.show()

a = 1
