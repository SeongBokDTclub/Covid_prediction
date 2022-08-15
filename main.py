import scipy.integrate
import numpy
import matplotlib.pyplot as plt

# Variable.
gamma = 0.07
R_ZERO = 10
beta = R_ZERO * gamma
INITIAL_S = 0.999
INITIAL_E = 0
INITIAL_I = 0.001
INITIAL_C = 0.001381
BIRTH_RATE = 0.00136354
DEATH_BY_TIME = 0.00106354
DEATH_BY_COVID = 0.001303421053
INCIDENCE_RATE = 0.5


def SIR(y, t, beta, gamma):
    S, I, R = y
    DIDT = beta * S * I - gamma * I
    DRDT = gamma * I
    DSDT = - beta * S * I

    return [DSDT, DIDT, DRDT]


def SEIR(y, t, beta, gamma):
    S, E, I, R = y
    DSDT = BIRTH_RATE - (beta * S * I + DEATH_BY_TIME) * S
    DEDT = beta * S * I - (DEATH_BY_TIME + INCIDENCE_RATE) * E
    DIDT = INCIDENCE_RATE * E - (DEATH_BY_TIME + gamma) * I
    DRDT = gamma * I - DEATH_BY_TIME * R

    return [DSDT, DEDT, DIDT, DRDT]


def SEIRS(y, t, beta, gamma):
    S, E, I, R = y

    C = INITIAL_C

    DSDT = BIRTH_RATE - (beta * I + DEATH_BY_TIME) * S + C * R
    DEDT = beta * S * I - (DEATH_BY_TIME + INCIDENCE_RATE) * E
    DIDT = INCIDENCE_RATE * E - (DEATH_BY_TIME + gamma) * I
    DRDT = gamma * I - (DEATH_BY_TIME + C) * R

    return [DSDT, DEDT, DIDT, DRDT]


t = numpy.linspace(0, 2000, 1000)
response = numpy.array(scipy.integrate.odeint(SIR, [INITIAL_S, INITIAL_I, 0], t, args=(beta, gamma)))
response_2 = numpy.array(
    scipy.integrate.odeint(SEIR, [INITIAL_S, INITIAL_E, INITIAL_I, 0], t, args=(beta, gamma)))
response_3 = numpy.array(
    scipy.integrate.odeint(SEIRS, [INITIAL_S, INITIAL_E, INITIAL_I, 0], t, args=(beta, gamma))
)

plt.figure(figsize=(10, 10))

plt.plot(t, response_3[:, 0], label="susceptible")
plt.plot(t, response_3[:, 1], label="exposed")
plt.plot(t, response_3[:, 2], label="infections")
plt.plot(t, response_3[:, 3], label="recoveries")


# plt.plot(t, response_2[:, 0], label="susceptible")
# plt.plot(t, response_2[:, 1], label="exposed")
# plt.plot(t, response_2[:, 2], label="infections")

#
# plt.plot(t, response[:,0], label="susceptible")
# plt.plot(t, response[:,1], label="infections")
# plt.plot(t, response[:,2], label="recoveries")

plt.legend()
plt.show()
