@np.vectorize
def R1(qtest):
    if(.05<=qtest<=2):
        return(.38 - .2 * np.log(qtest))
    elif(qtest>2):
        return(.426 / (1+qtest)**(1/3))
    else:
        return(0)

q = np.logspace(-1, 1)

Rcirc = (.5 - .227 * np.log(q) )**4 * (1+q)
plt.plot(q, R1(q), label='R1')
plt.plot(q, Rcirc, label = 'Rcirc')
plt.legend()
plt.show()
