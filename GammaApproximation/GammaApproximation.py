from sympy import *
#Code starts from here
def GammaAtGivenScale(m, x):
    gammaCdf = lowergamma(m,x)/gamma(m)
    #gammaCdf = integrate(exp(-x) * x**(m-1), (x, 0, x))/gamma(m)
    return gammaCdf

def GetMaxValue(f, x, useRound):
    fPrime = diff(f, x)

    solutions = None
    try:
        solutions = solveset(fPrime,x)
    except Exception as ex:
        print(ex)
        solutions = [0]

    max = float('-inf') 
    for solution in solutions:
        value = None
        if useRound:
            value = round(abs(f.subs(x, solution).evalf()),3)
        else:
            value = abs(f.subs(x, solution).evalf())
        if value > max:
            max = value
    return max

def GetThetaOneThetaTwoUpperBoundsForGammaApproximation(ms, useRound):
    x = symbols('x')
    mGammaDict = {}
    thetaOneThetaTwoUpperBounds = []
    for m in ms:
        gammaOfm = None
        if m in mGammaDict.keys():
            gammaOfm = mGammaDict[m]
        else:
            gammaOfm = GammaAtGivenScale(m, x)

        gammaOfmPlusOne = None
        if m+1 in mGammaDict.keys():
            gammaOfmPlusOne = mGammaDict[m+1]
        else:
            gammaOfmPlusOne = GammaAtGivenScale(m + 1, x)

        gammaOfmPlusTwo = None
        if m+2 in mGammaDict.keys():
            gammaOfmPlusTwo = mGammaDict[m+2]
        else:
            gammaOfmPlusTwo = GammaAtGivenScale(m + 2, x)

        gammaOfmPlusThree = None
        if m+3 in mGammaDict.keys():
            gammaOfmPlusThree = mGammaDict[m+3]
        else:
            gammaOfmPlusThree = GammaAtGivenScale(m + 3, x)

        gammaOfmPlusFour = None
        if m+4 in mGammaDict.keys():
            gammaOfmPlusFour = mGammaDict[m+4]
        else:
            gammaOfmPlusFour = GammaAtGivenScale(m + 4, x)

        thetaOne = -gammaOfm + 3*gammaOfmPlusOne - 3*gammaOfmPlusTwo + gammaOfmPlusThree
        thetaTwo = gammaOfm - 4*gammaOfmPlusOne + 6*gammaOfmPlusTwo - 4*gammaOfmPlusThree + gammaOfmPlusFour
        maxOfThetaOne = GetMaxValue(thetaOne, x, useRound)
        maxOfThetaTwo = GetMaxValue(thetaTwo, x, useRound)
        thetaOneThetaTwoUpperBounds.append((m, maxOfThetaOne, maxOfThetaTwo))
    
    return thetaOneThetaTwoUpperBounds

def CalculatePowerProfile(P_nkValues, P_nkPower, summationPower):
    sum = 0
    for P_nkValue in P_nkValues:
        sum = sum + (P_nkValue**P_nkPower)
    result = sum**summationPower
    return result

def CalculateCentralMoments(momentOrder, Pnks):
    if momentOrder > 4 or momentOrder < 1:
        raise ValueError(momentOrder)
    momentValue = None
    
    sumOfPnks = CalculatePowerProfile(Pnks, 1, 1)
    sumOfSquaredPnks = CalculatePowerProfile(Pnks, 2, 1)
    rho = sumOfPnks/sumOfSquaredPnks

    if momentOrder == 1:
        momentValue = rho * sumOfPnks
    elif momentOrder == 2:
        momentValue = (rho**2)*sumOfSquaredPnks
    elif momentOrder == 3:
        momentValue = 2*(rho**3) * CalculatePowerProfile(Pnks, 3, 1)
    else:
        momentValue = (rho**4)*(6*CalculatePowerProfile(Pnks, 4, 1) + 3*CalculatePowerProfile(Pnks, 2, 2))

    return momentValue

def CalculateAlphaConstantValue(m, Pnks):
    thirdCentralMoment = CalculateCentralMoments(3, Pnks)
    result = (thirdCentralMoment - 2*m)/factorial(3)
    return result

def CalculateGammaConstantValue(m, Pnks):
    thirdCentralMoment = CalculateCentralMoments(3, Pnks)
    fourthCentralMoment = CalculateCentralMoments(4, Pnks)
    result = (fourthCentralMoment - (12*thirdCentralMoment) - (3*(m**2)) + (18*m))/factorial(4)
    return result

def PrintSecondScenarioErrorUpperBound():
    Pnks = []
    for i in range(1,51):
        Pnks.append(10)#10dB
    for i in range(1,51):
        Pnks.append(10**(1.2))#12dB

    thetasUpperBounds = GetThetaOneThetaTwoUpperBoundsForGammaApproximation([96], False)
    for thetasUpperBound in thetasUpperBounds:
        print(f"Theta1: m = {thetasUpperBound[0]}, value = {thetasUpperBound[1]}")
        print(f"Theta2: m = {thetasUpperBound[0]}, value = {thetasUpperBound[2]}")
        alphaConst = CalculateAlphaConstantValue(thetasUpperBound[0],Pnks)
        print(f"Alpha Constant: m = {thetasUpperBound[0]}, value = {alphaConst}")
        gammaConst = CalculateGammaConstantValue(thetasUpperBound[0], Pnks)
        print(f"Gamma Constant: m = {thetasUpperBound[0]}, value = {gammaConst}")
        errorUpperBound = abs(alphaConst*thetasUpperBound[1] + gammaConst*thetasUpperBound[2])
        print(f"Error Upper Bound: m = {thetasUpperBound[0]}, value = {errorUpperBound}")

def PrintThirdScenarioErrorUpperBound():
    Pnks = [10]#10dB
    for i in range(2,101):
        Pnks.append(1)#0dB
    thetasUpperBounds = GetThetaOneThetaTwoUpperBoundsForGammaApproximation([60], False)
    for thetasUpperBound in thetasUpperBounds:
        print(f"Theta1: m = {thetasUpperBound[0]}, value = {thetasUpperBound[1]}")
        print(f"Theta2: m = {thetasUpperBound[0]}, value = {thetasUpperBound[2]}")
        alphaConst = CalculateAlphaConstantValue(thetasUpperBound[0],Pnks)
        print(f"Alpha Constant: m = {thetasUpperBound[0]}, value = {alphaConst}")
        gammaConst = CalculateGammaConstantValue(thetasUpperBound[0], Pnks)
        print(f"Gamma Constant: m = {thetasUpperBound[0]}, value = {gammaConst}")
        errorUpperBound = abs(alphaConst*thetasUpperBound[1] + gammaConst*thetasUpperBound[2])
        print(f"Error Upper Bound: m = {thetasUpperBound[0]}, value = {errorUpperBound}")

thetas = GetThetaOneThetaTwoUpperBoundsForGammaApproximation([1, 2, 4, 6, 8, 10], True)
for theta in thetas:
    print(f"Theta1: m = {theta[0]}, value = {theta[1]}")
    print(f"Theta2: m = {theta[0]}, value = {theta[2]}")

print("-"*25)
PrintSecondScenarioErrorUpperBound()
print("-"*25)
PrintThirdScenarioErrorUpperBound()

