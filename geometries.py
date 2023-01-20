import numpy as np

'''
Ideal geometries are generated in this file

Given the coordination number and an average distance from the central ion
to the first coordination sphere atoms (from the maximum of the rdf), the
geometries are created as numpy arrays
'''

def idealGeometries(coordNum, r):
    geometries = []


    if (coordNum == 3):
        # trigonal planar (TP)
        a = 2 * np.pi / 3
        TP = np.array([
            [r*np.cos(0*a), r*np.sin(0*a), 0.],
            [r*np.cos(1*a), r*np.sin(1*a), 0.],
            [r*np.cos(2*a), r*np.sin(2*a), 0.]
        ])
        geometries.append(['trigonal-planar', TP])




    elif (coordNum == 4):
        # square planar (SP)
        SP = np.array([
            [r, 0., 0.],
            [0., r, 0.],
            [-r, 0., 0.],
            [0., -r, 0.]
        ])
        geometries.append(['square-planar', SP])

        # tetrahedral (T)
        theta = np.pi - 2*np.arctan(np.sqrt(2))
        x = r * np.sin(theta)
        y = r * np.cos(theta)
        T = np.array([
            [0., r, 0.],
            [x, -y, 0.],
            [-0.5*x, -y, np.sqrt(3)/2*x],
            [-0.5*x, -y, -np.sqrt(3)/2*x]
        ])
        geometries.append(['tetrahedral', T])



    elif (coordNum == 5):
        # square pyramidal (SP)
        SP = np.array([
            [r, 0., 0.],
            [0., r, 0.],
            [-r, 0., 0.],
            [0., -r, 0.],
            [0., 0., r]
        ])
        geometries.append(['square-pyramidal', SP])

        # trigonal bipyramidal (TB)
        a = 2 * np.pi / 3
        TB = np.array([
            [0., 0., r],
            [0., 0., -r],
            [r*np.cos(0*a), r*np.sin(0*a), 0.],
            [r*np.cos(1*a), r*np.sin(1*a), 0.],
            [r*np.cos(2*a), r*np.sin(2*a), 0.]
        ])
        geometries.append(['trigonal-bipyramidal', TB])



    elif (coordNum == 6):
        # octahedral (O)
        O = np.array([
            [r, 0., 0.],
            [0., r, 0.],
            [0., 0., r],
            [-r, 0., 0.],
            [0., -r, 0.],
            [0., 0., -r]
        ])
        geometries.append(['octahedral', O])

        # Trigonal Prismatic (TP)
        theta = np.arctan(2 / np.sqrt(3))
        sin = np.sin(theta)
        cos = np.cos(theta)
        TP = np.array([
            [0., r*sin, r*cos],
            [0., r*sin, -r*cos],
            [-r*sin*np.sqrt(3)/2, -r*sin*0.5, r*cos],
            [-r*sin*np.sqrt(3)/2, -r*sin*0.5, -r*cos],
            [r*sin*np.sqrt(3)/2, -r*sin*0.5, r*cos],
            [r*sin*np.sqrt(3)/2, -r*sin*0.5, -r*cos]
        ])
        geometries.append(['trigonal-prismatic', TP])




    elif (coordNum == 7):
        # capped octahedral (CO)
        theta = np.atan(np.sqrt(0.5))
        a = (37 - 3*np.power(np.cos(theta), 2) + 12*np.sin(theta))
        b = r*(-14*np.sin(theta) - 84)
        c = 49*pow(r, 2)
        r1 = (-b - np.sqrt(b**2 - 4*a*c)) / (2*a)
        r2 = 7*r - 6*r1
        x = r1*np.cos(theta)
        y = r1*np.sin(theta)
        CO = np.array([
            [x, y, 0.],
            [-0.5*x, y, np.sqrt(3)/2],
            [-0.5*x, y, -np.sqrt(3)/2*x],
            [-x, -y, 0.],
            [0.5*x, -y, np.sqrt(3)/2*x],
            [0.5*x, -y, -np.sqrt(3)/2*x],
            [0., r2, 0.]
        ])
        geometries.append(['capped-octahedral', CO])

        # capped trigonal prismatic (CTP)
        theta = np.arctan(2 / np.sqrt(3))
        a = 9 + 3*np.sin(theta) - (5/7)
        b = (-24 - 4*np.sin(theta)) * r
        c = 16 * np.power(r, 2)
        r1 = (-b - np.sqrt(np.power(b, 2) - 4*a*c)) / (2*a)
        r2 = 4*r - 3*r1
        sin = np.sin(theta)
        cos = np.cos(theta)
        CTP = np.array([
            [0., r1*sin, r1*cos],
            [0., r1*sin, -r1*cos],
            [-r1*sin*np.sqrt(3)/2, -r1*sin*0.5, r1*cos],
            [-r1*sin*np.sqrt(3)/2, -r1*sin*0.5, -r1*cos],
            [r1*sin*np.sqrt(3)/2, -r1*sin*0.5, r1*cos],
            [r1*sin*np.sqrt(3)/2, -r1*sin*0.5, -r1*cos],
            [r2*np.sqrt(3)/2, r2*0.5, 0.]
        ])
        geometries.append(['capped-trigonal-prismatic', CTP])

        # pentagonal bipyramidal (BP)
        a = 2 * np.pi / 5
        PB = np.array([
            [0., 0., r1],
            [0., 0., -r1],
            [r*np.cos(0*a), r*np.sin(0*a), 0.],
            [r*np.cos(1*a), r*np.sin(1*a), 0.],
            [r*np.cos(2*a), r*np.sin(2*a), 0.],
            [r*np.cos(3*a), r*np.sin(3*a), 0.],
            [r*np.cos(4*a), r*np.sin(4*a), 0.]
        ])
        geometries.append(['pentagonal-bipyramidal', PB])




    elif (coordNum == 8):
        # cubic (C)
        x = r / np.sqrt(3)
        C = np.array([
            [x, x, x],
            [-x, x, x],
            [x, -x, x],
            [-x, -x, x],
            [x, x, -x],
            [-x, x, -x],
            [x, -x, -x],
            [-x, -x, -x]
        ])
        geometries.append(['cubic', C])

        # square antiprismatic (SA)
        theta = np.arctan(np.sqrt(np.sqrt(2) / 4))
        x = r*np.cos(theta)
        y = r*np.sin(theta)
        s = np.sqrt(2) / 2
        SA = np.array([
            [x, y, 0.],
            [-x, y, 0.0],
            [0., y, x],
            [0., y, -x],
            [s*x, -y, s*x],
            [-s*x, -y, s*x],
            [s*x, -y, -s*x],
            [-s*x, -y, -s*x]
        ])
        geometries.append(['square-antiprismatic', SA])

        # bicapped trigonal prismatic (BTP)
        theta = np.arctan(2 / np.sqrt(3))
        a = 9 + 3*np.sin(theta) - (5/7)
        b = (-24 - 4*np.sin(theta)) * r
        c = 16 * np.power(r, 2)
        r1 = (-b - np.sqrt(np.power(b, 2) - 4*a*c)) / (2*a)
        r2 = 4*r - 3*r1
        sin = np.sin(theta)
        cos = np.cos(theta)
        BTP = np.array([
            [0., r1*sin, r1*cos],
            [0., r1*sin, -r1*cos],
            [-r1*sin*np.sqrt(3)/2, -r1*sin*0.5, r1*cos],
            [-r1*sin*np.sqrt(3)/2, -r1*sin*0.5, -r1*cos],
            [r1*sin*np.sqrt(3)/2, -r1*sin*0.5, r1*cos],
            [r1*sin*np.sqrt(3)/2, -r1*sin*0.5, -r1*cos],
            [r2*np.sqrt(3)/2, r2*0.5, 0.],
            [-r2*np.sqrt(3)/2, r2*0.5, 0.]
        ])
        geometries.append(['bicapped-trigonal-prismatic', BTP])

        # dodecahedral (D)
        theta_A = 0.64315383 # angle in rads
        theta_B = 1.2123057 # angle in rads
        z_A = r * np.cos(theta_A)
        x_A = r * np.sin(theta_A)
        z_B = r * np.cos(theta_B)
        x_B = r * np.sin(theta_B)
        D = np.array([
            [x_A, 0., z_A],
            [-x_A, 0., z_A],
            [0., x_A, -z_A],
            [0., -x_A, -z_A],
            [x_B, 0., -z_B],
            [-x_B, 0., -z_B],
            [0., x_B, z_B],
            [0., -x_B, z_B]
        ])
        geometries.append(['dodecahedral', D])

        # hexagonal bipyramidal (HB)
        a = np.pi / 3
        HB = np.array([
            [0., 0., r],
            [0., 0., -r],
            [r*np.cos(0*a), r*np.sin(0*a), 0.],
            [r*np.cos(1*a), r*np.sin(1*a), 0.],
            [r*np.cos(2*a), r*np.sin(2*a), 0.],
            [r*np.cos(3*a), r*np.sin(3*a), 0.],
            [r*np.cos(4*a), r*np.sin(4*a), 0.],
            [r*np.cos(5*a), r*np.sin(5*a), 0.]
        ])
        geometries.append(['hexagonal-bipyramidal', HB])

        # trans bicapped octahedral (TBO)
        theta = np.arctan(np.sqrt(0.5))
        a = (10 - 3*np.power(np.cos(theta), 2) + 6*np.sin(theta))
        b = r*(-8*np.sin(theta) - 24)
        c = 16*np.power(r, 2)
        r1 = (-b - np.sqrt(b**2 - 4*a*c)) / (2*a)
        r2 = 4*r - 3*r1
        x = r1*np.cos(theta)
        y = r1*np.sin(theta)
        TBO = np.array([
            [x, y, 0.],
            [-0.5*x, y, np.sqrt(3)/2*x],
            [-0.5*x, y, -np.sqrt(3)/2*x],
            [-x, -y, 0.],
            [0.5*x, -y, np.sqrt(3)/2*x],
            [0.5*x, -y, -np.sqrt(3)/2*x],
            [0., r2, 0.],
            [0., -r2, 0.]
        ])
        geometries.append(['trans-bicapped-octahedral', TBO])




    elif (coordNum == 9):
        # capped square antiprismatic (CSA)
        theta = np.arctan(np.sqrt(np.sqrt(2)/4))
        a = (65 - 2*np.power(np.cos(theta), 2) + 16*np.sin(theta))
        b = (-18*r*np.sin(theta) - 144*r)
        c = 81*np.power(r, 2)
        r1 = (-b - np.sqrt(b**2 - 4*a*c)) / (2*a)
        r2 = 9*r - 8*r1
        x = r1*np.cos(theta)
        y = r1*np.sin(theta)
        s = np.sqrt(2)/2
        CSA = np.array([
            [x, y, 0.],
            [-x, y, 0.],
            [0., y, x],
            [0., y, -x],
            [s*x, -y, s*x],
            [-s*x, -y, s*x],
            [s*x, -y, -s*x],
            [-s*x, -y, -s*x],
            [0., r2, 0.]
        ])
        geometries.append(['capped-square-antiprismatic', CSA])

        # capped square (CS)
        theta = np.arctan(np.sqrt(2) / 2)
        a = (65 - 2*np.power(np.cos(theta), 2) + 16*np.sin(theta))
        b = (-18*r*np.sin(theta) - 144*r)
        c = 81*np.power(r, 2)
        r1 = (-b - np.sqrt(b**2 - 4*a*c)) / (2*a)
        r2 = 9*r - 8*r1
        x = r1*np.cos(theta)
        y = r1*np.sin(theta)
        CS = np.array([
            [x, y, 0.],
            [-x, y, 0.],
            [0., y, x],
            [0., y, -x],
            [x, -y, 0.],
            [-x, -y, 0.],
            [0., -y, x],
            [0., -y, -x],
            [0., r2, 0.]
        ])
        geometries.append(['capped-square', CS])

        # tricapped trigonal prismatic (TTP)
        theta = np.arctan(2/np.sqrt(3))
        a = 4 + 2*np.sin(theta) - (5/7)
        b = (-12 - 3*np.sin(theta)) * r
        c = 9 * np.power(r, 2)
        r1 = (-b - np.sqrt(np.power(b, 2) - 4*a*c)) / (2*a)
        r2 = 3*r - 2*r1
        sin = np.sin(theta)
        cos = np.cos(theta)
        TTP = np.array([
            [0., r1*sin, r1*cos],
            [0., r1*sin, -r1*cos],
            [-r1*sin*np.sqrt(3)/2, -r1*sin*0.5, r1*cos],
            [-r1*sin*np.sqrt(3)/2, -r1*sin*0.5, -r1*cos],
            [r1*sin*np.sqrt(3)/2, -r1*sin*0.5, r1*cos],
            [r1*sin*np.sqrt(3)/2, -r1*sin*0.5, -r1*cos],
            [r2*np.sqrt(3)/2, r2*0.5, 0.],
            [-r2*np.sqrt(3)/2, r2*0.5, 0.],
            [0., -r2, 0.]
        ])
        geometries.append(['tricapped-trigonal-prismatic', TTP])


    

    elif (coordNum == 10):
        # bicapped square antiprismatic (BSA)
        theta = np.arctan(np.sqrt(np.sqrt(2)/4))
        a = (17 - 2*np.power(np.cos(theta), 2) + 8*np.sin(theta))
        b = (-10*r*np.sin(theta) - 40*r)
        c = 25*np.power(r, 2)
        r1 = (-b - np.sqrt(b**2 - 4*a*c)) / (2*a)
        r2 = 5*r - 4*r1
        x = r1*np.cos(theta)
        y = r1*np.sin(theta)
        s = np.sqrt(2)/2
        BSA = np.array([
            [x, y, 0.],
            [-x, y, 0.],
            [0., y, x],
            [0., y, -x],
            [s*x, -y, s*x],
            [-s*x, -y, s*x],
            [s*x, -y, -s*x],
            [-s*x, -y, -s*x],
            [0., r2, 0.],
            [0., -r2, 0.]
        ])
        geometries.append(['bicapped-square-antiprismatic', BSA])

        # bicapped square (BS)
        theta = np.arctan(np.sqrt(2)/2)
        a = (17 - 2*np.power(np.cos(theta), 2) + 8*np.sin(theta))
        b = (-10*r*np.sin(theta) - 40*r)
        c = 25*np.power(r, 2)
        r1 = (-b - np.sqrt(b**2 - 4*a*c)) / (2*a)
        r2 = 5*r - 4*r1
        x = r1*np.cos(theta)
        y = r1*np.sin(theta)
        BS = np.array([
            [x, y, 0.],
            [-x, y, 0.],
            [0., y, x],
            [0., y, -x],
            [x, -y, 0.],
            [-x, -y, 0.],
            [0., -y, x],
            [0., -y, -x],
            [0., r2, 0.],
            [0., -r2, 0.]
        ])
        geometries.append(['bicapped-square', BS])

        # octagonal bipyramidal (OB)
        a = np.pi / 4
        OB = np.array([
            [0., 0., r],
            [0., 0., -r],
            [r*np.cos(0*a), r*np.sin(0*a), 0.],
            [r*np.cos(1*a), r*np.sin(1*a), 0.],
            [r*np.cos(2*a), r*np.sin(2*a), 0.],
            [r*np.cos(3*a), r*np.sin(3*a), 0.],
            [r*np.cos(4*a), r*np.sin(4*a), 0.],
            [r*np.cos(5*a), r*np.sin(5*a), 0.],
            [r*np.cos(6*a), r*np.sin(6*a), 0.],
            [r*np.cos(7*a), r*np.sin(7*a), 0.]
        ])
        geometries.append(['octagonal-bipyramidal', OB])

        # pentagonal antiprismatic (PA)
        a = np.arctan( np.sqrt((1 - np.cos(2*np.pi / 5)) / 2) )
        i = 2*np.pi / 5
        PA = np.array([
            [r*np.cos(a), 0., r*np.sin(a)],
            [-r*np.cos(a), 0., -r*np.sin(a)],
            [r*np.cos(a)*np.cos(i), r*np.cos(a)*np.sin(i), r*np.sin(a)],
            [-r*np.cos(a)*np.cos(i), r*np.cos(a)*np.sin(i), -r*np.sin(a)],
            [r*np.cos(a)*np.cos(i), -r*np.cos(a)*np.sin(i), r*np.sin(a)],
            [-r*np.cos(a)*np.cos(i), -r*np.cos(a)*np.sin(i), -r*np.sin(a)],
            [r*np.cos(a)*np.cos(2*i), r*np.cos(a)*np.sin(2*i), r*np.sin(a)],
            [-r*np.cos(a)*np.cos(2*i), r*np.cos(a)*np.sin(2*i), -r*np.sin(a)],
            [r*np.cos(a)*np.cos(2*i), -r*np.cos(a)*np.sin(2*i), r*np.sin(a)],
            [-r*np.cos(a)*np.cos(2*i), -r*np.cos(a)*np.sin(2*i), -r*np.sin(a)]
        ])
        geometries.append(['pentagonal-antiprismatic', PA])

        # pentagonal prismatic (PP)
        a = np.arctan( np.sqrt((1 - np.cos(2*np.pi / 5)) / 2) )
        i = 2*np.pi / 5
        PP = np.array([
            [r*np.cos(a), 0., r*np.sin(a)],
            [r*np.cos(a), 0., -r*np.sin(a)],
            [r*np.cos(a)*np.cos(i), r*np.cos(a)*np.sin(i), r*np.sin(a)],
            [r*np.cos(a)*np.cos(i), r*np.cos(a)*np.sin(i), -r*np.sin(a)],
            [r*np.cos(a)*np.cos(i), -r*np.cos(a)*np.sin(i), r*np.sin(a)],
            [r*np.cos(a)*np.cos(i), -r*np.cos(a)*np.sin(i), -r*np.sin(a)],
            [r*np.cos(a)*np.cos(2*i), r*np.cos(a)*np.sin(2*i), r*np.sin(a)],
            [r*np.cos(a)*np.cos(2*i), -r*np.cos(a)*np.sin(2*i), r*np.sin(a)],
            [r*np.cos(a)*np.cos(2*i), r*np.cos(a)*np.sin(2*i), -r*np.sin(a)],
            [r*np.cos(a)*np.cos(2*i), -r*np.cos(a)*np.sin(2*i), -r*np.sin(a)]
        ])
        geometries.append(['pentagonal-prismatic', PP])




    elif (coordNum == 12):
        # Icosahedral (I)
        a = r / np.sin(2 * np.pi / 5) / 2
        phi = (1 + np.sqrt(5)) / 2
        I = np.array([
            [0., a, a*phi],
            [0., -a, a*phi],
            [0., a, -a*phi],
            [0., -a, -a*phi],
            [a, a*phi, 0.],
            [-a, a*phi, 0.],
            [a, -a*phi, 0.],
            [-a, -a*phi, 0.],
            [a*phi, 0., a],
            [-a*phi, 0., a],
            [a*phi, 0., -a],
            [-a*phi, 0., -a]
        ])
        geometries.append(['icosahedral', I])



    return geometries
