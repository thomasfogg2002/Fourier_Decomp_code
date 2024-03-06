
# GENERATE FOURIER DECOMPOSITION COEFFICIENTS AS FUNCTION OF s

import cmath
import numpy as np

order = 4   # Display coefficients up to 4th multipole
slen = 5  # Number of points along longitudinal coordinate
x = np.linspace(-20,20,41)
y = np.linspace(-20,20,41)
Bx = np.zeros((41,41,5))
By = np.zeros((41,41,5))




M = 24  # Number of measurements
r0 = 0.01  # Reference radius [m]

C_FD = []       # Real part
C_FD_imag = []  # Imaginary part

# Evaluate sums for each coefficient - done as function of longitudinal coordinate, s

for kref in range(slen):    # Generates the coeff's at each longitudinal point along the centre of the magnet

    Cn_values_real = []
    Cn_values_imag = []

    for n in range(1,order+1):

        sum_values = 0

        for m in range(M):

            theta_m = 2*m*np.pi/M

            xref,yref = [r0*np.cos(theta_m),r0*np.sin(theta_m)]

            if abs(xref)>max(x):    # Stops you increasing r0 beyond aperture
                break

            if xref>0:
                iref = max(np.where(x <= round(xref,5))[0])
            else:
                iref = min(np.where(x >= round(xref,5))[0])

            if yref>0:
                jref = max(np.where(y <= round(yref,5))[0])
            else:
                jref = min(np.where(y >= round(yref,5))[0])

            B_cart = complex(By[iref,jref,kref],Bx[iref,jref,kref])
            Bm = B_cart*cmath.rect(1,theta_m)

            value = (1/(M*r0**(n-1)))*Bm*cmath.rect(1,-(n)*theta_m)

            sum_values = sum_values + value

        Cn_values_real.append(sum_values.real)
        Cn_values_imag.append(sum_values.imag)

    C_FD.append(Cn_values_real)
    C_FD_imag.append(Cn_values_imag)

C_FD = np.transpose(C_FD)
C_FD_imag = np.transpose(C_FD_imag)

print('Dipole coeff\'s: C1(s) = '+str(C_FD[0]))
print('Dipole coeff\'s: C2(s) = '+str(C_FD[1]))
print('... up to C_order = C'+str(order)+' in this case')