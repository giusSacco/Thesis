import numpy as np

def random_dir():
    vec = np.random.standard_normal(size=3)
    vec /= np.linalg.norm(vec)
    return vec

N = 15000

alpha_list = list(); beta_list = list(); gamma_list = list()
for i in range(N):
    alpha, beta, gamma = random_dir()
    alpha_list.append(alpha)
    beta_list.append(beta)
    gamma_list.append(gamma)

print(N, np.average(np.array(alpha_list)), np.average(np.array(beta_list)), np.average(np.array(gamma_list)))