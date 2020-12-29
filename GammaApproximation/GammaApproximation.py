from sympy import *
#Code starts from here
def GammaAtGivenScale(m, x):
    gammaCdf = lowergamma(m,x)/gamma(m)
    #gammaCdf = integrate(exp(-x) * x**(m-1), (x, 0, x))/gamma(m)
    return gammaCdf

def GetMaxValue(f, x):
    fPrime = diff(f, x)
    solutions = solve(fPrime,x)
    max = float('-inf') 
    for solution in solutions:
        value = round(abs(f.subs(x, solution).evalf()),3)
        if value > max:
            max = value
    return max

def PrintUpperBoundOfThetaOneAndThetaTwoForGammaApproximationErrorFormula():
    x = symbols('x')
    ms = [1,2,4,6,8]
    mGammaDict = {}
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
        maxOfThetaOne = GetMaxValue(thetaOne, x)
        maxOfThetaTwo = GetMaxValue(thetaTwo, x)
        print(f"Theta1: m = {m}, value = {maxOfThetaOne}")
        print(f"Theta2: m = {m}, value = {maxOfThetaTwo}\n\n")


PrintUpperBoundOfThetaOneAndThetaTwoForGammaApproximationErrorFormula()
