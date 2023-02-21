from joblib import Parallel, delayed
from Flow import Flow
import numpy as np


def execute_flow(mu, T):
    flow = Flow(1e5, 1e-4, 6, 1000, mu, T, 2)
    flow.compute()
    flow.get_observables_for_all_positions()
    flow.save("phase_diagram")


if __name__ == "__main__":
    mu_array = np.arange(0.0, 0.8, 0.0125)
    T_array = np.arange(0.01, 0.7, 0.0125)
    job_list = []
    for mu in mu_array:
        for T in T_array:
            job_list.append(delayed(execute_flow)(mu, T))
    Parallel(n_jobs=-1)(job_list)
