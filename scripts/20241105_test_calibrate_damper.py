import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, fsolve

T1 = [2.11319, 2.11481, 2.11643, 2.11788, 2.11926, 2.12097, 2.12239, 2.12407, 2.12572, 2.12708, 2.13041]  # b_c = 0.404371245546549
A1 = [0.0773, 0.0618, 0.0519, 0.0508, 0.0455, 0.0341, 0.0336, 0.027, 0.0256, 0.0256, 0.0229]
T2 = [1.094, 1.09499, 1.0954, 1.09626, 1.09721, 1.09852, 1.09943, 1.10013, 1.10058, 1.10161, 1.10248, 1.10478, 1.10634]  # b_c = 0.4284206885466512
A2 = [0.0712, 0.0527, 0.0610, 0.0518, 0.0493, 0.0420, 0.0393, 0.0375, 0.0372, 0.0318, 0.0328, 0.0261, 0.0225]
T3 = [3.65089, 3.65183, 3.65269, 3.65331, 3.65394, 3.65480, 3.65566, 3.65723, 3.65871, 3.65950, 3.66036, 3.66129, 3.66771]  # b_lc = 0.3252700179281933
A3 = [0.072, 0.0656, 0.0548, 0.0492, 0.0611, 0.0513, 0.0473, 0.0469, 0.0439, 0.0295, 0.0255, 0.0270, 0.0196]

m_lcs = []
b_lcs = []
for T, A in zip([T1, T2, T3], [A1, A2, A3]):
    T = np.array(T) - T[0]
    A = np.array(A) - A[-1]

    f = lambda t, amp, tau: amp * np.exp(-t / tau)
    popt, pcov = curve_fit(f, T, A, p0=(A[0], 1 / 300), bounds=(0, (1, 1000)))
    print('amp', popt[0], 'tau', popt[1])
    tau = popt[1]

    k_lc = 18.017806799339144 * 1e3  # N/mm
    omega_lc = 2 * np.pi * 643.4036  # calcualted on 20241029 Daily Notes


    def equations(X):
        m_lc, b_lc = X
        return (omega_lc**2 - (k_lc / m_lc) + np.square(b_lc / 2 / m_lc),
                1 / tau - b_lc / 2 / m_lc)


    # m_lc = k_lc / np.square(omega_lc)  # N/mm --> N/m / (1/s^2)
    # b_lc = 2 * m_lc / popt[1]
    m_lc, b_lc = fsolve(equations, (1e-3, 0.5))
    print('k_lc = {}, m_lc = {}, b_lc = {}'.format(k_lc, m_lc, b_lc))

    m_lcs.append(m_lc)
    b_lcs.append(b_lc)
print("m_lc_avg", np.mean(m_lcs))
print("b_lc_avg", np.mean(b_lcs))
