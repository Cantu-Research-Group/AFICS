import numpy as np
import csv
import geometries
import itertools


class Trajectory ():
    def __init__(
                    self,
                    atoms=[],
                    atomNum=0,
                    frameCount=0,
                    ionID='',
                    ionNum=0,
                    whichIon=1,
                    elements=[],
                    boxSize=10.0,
                    rdf=[],
                    rdfIntegral=[],
                    adf=[],
                    binSize=0.05,
                    dist=0.0,
                    coorNum=0,
                    maxR=0.0,
                    missingFrames=[],
                    RMSDs=[],
                    framesForRMSD=100,
                    idealGeos=[],
                    startFrame=None,
                    endFrame=100000000
                    ):
        self.atoms = atoms              # array of atoms and their respective coordinates
        self.atomNum = atomNum          # number of atoms per frame
        self.frameCount = frameCount    # number of frames in the trajectory
        self.ionID = ionID              # string of the name of the central ion for rdf, ie. 'Lu' (needs to changed to accomadate more than one of the same ion)
        self.ionNum = ionNum            # number of ions with the ion name, ionID
        self.whichIon = whichIon        # if there are multiple of the same ion, this specifies which one to center the rdf around (not indexed at 0), by default gets 1st one
        self.elements = elements        # elements to include in the rdf/adf ie. 'O', 'N'
        self.boxSize = boxSize          # size of the cell
        self.rdf = rdf                  # array for rdf bins used in histogram
        self.rdfIntegral = rdfIntegral  # array for integral of rdf used in histogram
        self.adf = adf
        self.binSize = binSize          # size of bin in histogram
        self.dist = dist                # average distance in first sphere, ie max of rdf
        self.coorNum = coorNum          # coordination number
        self.maxR = maxR                # maximum distance of atoms within the first coordination sphere
        self.thresholdAtoms = []        # atoms within the first coordination sphere
        self.missingFrames = missingFrames  # index of frames that were not able to be collected
        self.RMSDs = RMSDs              # array holding RMSD values
        self.framesForRMSD = framesForRMSD  # number of frames each optimal permutation is used for, higher numbers make the program faster
        self.idealGeos = idealGeos      # array containing the ideal geometries
        self.startFrame = startFrame
        self.endFrame = endFrame

    def getRMSDs(self):
        # [ ['geometry name', [[x,y,z],[x,y,z],...], ... ]
        print ('Calculating RMSDs... (this may take a while)')
        geoNames = []
        for ideal in self.idealGeos:
            geoNames.append(ideal[0])
        self.RMSDs.append(geoNames)
        for i, frame in enumerate(self.thresholdAtoms):
            if (i % self.framesForRMSD == 0):
                self.idealGeos = reorderIdealGeometries(frame, self.idealGeos)
            tempRMSDs = []
            for ideal in self.idealGeos:
                tempRMSDs.append(str(kabschRMSD(frame, ideal[1])))
            self.RMSDs.append(tempRMSDs.copy())
        print ('RMSDs done')

    def getAtoms(self, filename):
        '''
        Retrieves atoms from an .xyz file and saves all the elements and their respective
        coordinates per frame as an array in the form:
        [ [ ['element', [x, y, z]], ['element', [x, y, z]], ... ]... ] 
        '''
        with open(filename, 'r') as f:
            frame = []
            if self.startFrame is not None and self.endFrame: # if startFrame and endFrame are numbers
                if self.startFrame < 1:
                    print ('Starting frame must be greater than 0')
                    return None
                if self.startFrame < self.endFrame:
                    if self.endFrame != 100000000:
                        print ('Collecting atoms from frame '+str(self.startFrame)+' to frame '+str(self.endFrame)+'...')
                    else:
                        print ('Collecting atoms starting from frame'+str(self.startFrame)+'...')
                    framesCounter = 0
                    for i, row in enumerate(f):
                        row = row.split()
                        if i == 0:
                                self.atomNum = int(row[0])
                        if (i % (self.atomNum + 2) == 0):
                            if row: # ensures its not an empty line at the end of the file
                                framesCounter += 1
                        if (framesCounter >= self.startFrame and framesCounter <= (self.endFrame + 1)):
                            if (i % (self.atomNum + 2) != 0) and (i % (self.atomNum + 2) != 1): # + 2 for the atom count and comment line of .xyz
                                frame.append([row[0], [float(row[1]), float(row[2]), float(row[3])]])
                            if (i % (self.atomNum + 2) == 0) and (self.frameCount != 0):
                                self.atoms.append(frame.copy())
                                frame.clear()
                            if (i % (self.atomNum + 2) == 0):
                                if row:
                                    self.frameCount += 1
                else:
                    print ('End frame must be a larger number than startframe')
            else:
                print ('Collecting all atoms...')
                for i, row in enumerate(f):
                    row = row.split()
                    if i == 0:
                        self.atomNum = int(row[0])
                    if (i % (self.atomNum + 2) != 0) and (i % (self.atomNum + 2) != 1): # + 2 for the atom count and comment line of .xyz
                        frame.append([row[0], [float(row[1]), float(row[2]), float(row[3])]])
                    if (i % (self.atomNum + 2) == 0) and (i != 0):
                        self.atoms.append(frame.copy())
                        frame.clear()
                    if (i % (self.atomNum + 2) == 0):
                        if row:
                            self.frameCount += 1

        print ('Atoms collected')
    
    def getThresholdAtoms(self):
        '''
        Retrieves atoms in the first coordinate sphere of the ion
        Format:
        [ [ [x, y, z], [x, y, z], ... ], ... ]
        '''
        for i, frame in enumerate(self.atoms):
            frameAtoms = []
            ionCoords = []
            ionCount = 0
            for atom in frame:
                if atom[0] == self.ionID:
                    ionCount += 1
                    if self.whichIon == ionCount:
                        ionCoords = atom[1]
                        break
            for atom in frame:
                if (atom[0] in self.elements) and (dist_eq(atom[1], ionCoords) <= self.maxR):
                    frameAtoms.append(atom[1])
            if (len(frameAtoms) == self.coorNum):
                self.thresholdAtoms.append(frameAtoms.copy())
                frameAtoms.clear()
            else:
                frameAtoms.clear()
                self.missingFrames.append(i + self.startFrame)
                print ('Frame '+str(i + self.startFrame)+' was unable to be collected')
        self.thresholdAtoms = np.array(self.thresholdAtoms)
        print ('First sphere atoms complete')

    def getDist(self):
        '''
        Calculates the maximum of the RDF to be used as the average distance
        in generating ideal structures
        '''
        max = 0.0
        maxIndex = 0
        count = 0
        breakCount = 0
        for i, value in enumerate(self.rdf):		
            if (value <= max) and (max != 0.0):
                breakCount += 1
            if (value > max):
                max = value
                maxIndex = i
            if (breakCount > 5): # 5 seems to work well for break count with binSize 0.02 or larger, 10 works for binSize 0.01
                self.dist = maxIndex*self.binSize
                break
            count += 1

    def getMaxR(self):
        '''
        Calculates the maximum radius for determining atoms within the first coordination sphere
        '''
        max = 0.0
        breakCount = 0
        for i, value in enumerate(self.rdfIntegral):
            if (i != 0):
                if (value - self.rdfIntegral[i - 1] < 0.02) and (value > 0):
                    breakCount += 1
            if (value > max):
                max = value
            if (breakCount >= 3):
                self.maxR = self.binSize * (i + 1) # plus 1 to get the outside of the bin
                break

    def getCN(self):
        '''
        Calculates coordination number
        '''
        max = 0.0
        breakCount = 0
        for i, value in enumerate(self.rdfIntegral):
            if (i != 0):
                if (value - self.rdfIntegral[i - 1] < 0.02) and (value > 0):
                    breakCount += 1
            if (value > max):
                max = value
            if (breakCount > 3):
                self.coorNum = round(max)
                break

    def printRDF(self, filename):
        '''
        Outputs RDF and RDF integral as .csv files
        '''
        RDFHeader = ['Distance (A)', 'Count']
        RDFIntHeader = ['Distance (A)', 'Integral']
        RDFout = []
        RDFIntout = []
        for i, value in enumerate(self.rdf):
            RDFout.append([(i * self.binSize), value])
        for i, value in enumerate(self.rdfIntegral):
            RDFIntout.append([(i * self.binSize), value])
        with open(filename+'-rdf.csv', 'w') as rdffile:
            rdfwriter = csv.writer(rdffile)
            rdfwriter.writerow(RDFHeader)
            rdfwriter.writerows(RDFout)
        with open(filename+'-rdf-integral.csv', 'w') as rdffile:
            rdfwriter = csv.writer(rdffile)
            rdfwriter.writerow(RDFIntHeader)
            rdfwriter.writerows(RDFIntout)

    def printADF(self, filename):
        '''
        Outputs RDF as .csv files
        '''
        ADFdict = {}
        for i in range(180):
            ADFdict[i] = 0 # initializes 180 entries for adf histogram of bin size 1 degree
        for value in self.adf:
            ADFdict[int(np.floor(value))] += 1
        with open(filename+'-adf.csv', 'w') as adffile:
            adfwriter = csv.writer(adffile)
            adfwriter.writerow(['Angle', 'Count'])
            for key in ADFdict:
                adfwriter.writerow([key, ADFdict[key]])

    def printRMSDs(self, filename):
        '''
        Outputs RMSDs as a .csv file
        '''
        with open(filename+'-RMSDs.csv', 'w') as rmsdfile:
            rmsdfwriter = csv.writer(rmsdfile)
            rmsdfwriter.writerow(self.RMSDs[0])
            RMSDs_float = np.array(self.RMSDs[1:], dtype=float)
            avg = np.mean(RMSDs_float, axis=0)
            stdDev = np.std(RMSDs_float, axis=0)
            avgHeader = []
            stdDevHeader = []
            for value in avg:
                avgHeader.append('Avg: '+str('{:.3f}'.format(value)))
            for value in stdDev:
                stdDevHeader.append('StdDev: '+str('{:.3f}'.format(value)))
            rmsdfwriter.writerow(avgHeader)
            rmsdfwriter.writerow(stdDevHeader)
            rmsdfwriter.writerows(self.RMSDs[1:])

    def outputIdealGeometries(self, subfolder=''):
        if len(subfolder):
            subfolder += '/'
        for ideal in self.idealGeos:
            header = [
                [str(int(self.coorNum + 1))], # add 1 for the central 'ion'
                ['CN='+str(int(self.coorNum))+': '+str(ideal[0])] # comment line of xyz says the geometry name
            ]
            with open ('./'+str(subfolder)+str(ideal[0])+'.xyz', 'w') as file:
                csvwriter = csv.writer(file)
                csvwriter.writerows(header)
                csvwriter.writerow(['O 0.0 0.0 0.0'])
                for row in ideal[1]:
                    csvwriter.writerow(['H '+str(row[0])+' '+str(row[1])+' '+str(row[2])])

    def getIonNum(self):
        '''
        Finds the total number of ions specified by the ionID per frame
        '''
        for frame in self.atoms:
            for atom in frame:
                if atom[0] == self.ionID:
                    self.ionNum += 1
            break

    def getWhichIon(self):
        ionsCoords = []
        for frame in self.atoms:
            for atom in frame:
                if atom[0] == self.ionID:
                    ionsCoords.append(atom[1])
        print ('There are '+str(self.ionNum)+' '+str(self.ionID)+' ions.')
        print ('Enter which one you would like to use:')
        for i in range(self.ionNum):
            print (str(i + 1)+' - '+str(self.ionID)+' at ('+str(ionsCoords[i][0])+', '+str(ionsCoords[i][1])+', '+str(ionsCoords[i][2])+')')
        while (True):
            usr_input = input ('Enter which ion to use: ')
            try:
                if (int(usr_input) > 0 and int(usr_input) <= self.ionNum):
                    self.whichIon = int(usr_input)
                    break
                else:
                    print ('Enter an integer between 1 and '+str(self.ionNum))
            except:
                print ('Enter an integer')


    def getRDF(self):
        """
        Obtains the RDF and integral of the RDF using bins of 0.05 A width
        """
        print ('Calculating RDF... (this may take a while)')
        bin = self.binSize # bin size is not adjustable by user and is set at 0.05 A
        r_max = 6.0 # the max the rdf will go to is 6 A
        data_num = int(r_max / bin)

        vol = 4.0/3.0 * np.pi * np.power(r_max, 3)
        density_A = 1.0 / vol

        distances = []

        for frame in self.atoms:
            frame_dist = []
            ion_coords = []
            # get ion coordinates
            ionCount = 0
            for atom in frame:
                if atom[0] == self.ionID:
                    ionCount += 1
                    if self.whichIon == ionCount:
                        ion_coords = atom[1]
                        break
            for atom in frame:
                if atom[0] in self.elements:
                    frame_dist.append(dist_eq(ion_coords, atom[1]))
            distances.append(frame_dist.copy())
            frame_dist.clear()
        atom_B_tol_num = []
        for frame in distances:
            tmp_dist_count = 0
            for distance in frame:
                if distance <= r_max:
                    tmp_dist_count += 1
            atom_B_tol_num.append(tmp_dist_count)

        rdf_values = []
        for i, frame in enumerate(distances):
            tmp_rdf_values = []
            for binNum in range(data_num):
                num_B_r = 0.0
                for distance in frame:
                    if (distance >= binNum*bin) and (distance < (binNum + 1)*bin):
                        num_B_r += 1
                tmp_rdf_values.append(1.0 / (density_A*atom_B_tol_num[i])*num_B_r / (4.0*np.pi*(binNum + 1)**2 * bin**3))
            rdf_values.append(tmp_rdf_values.copy())
            tmp_rdf_values.clear()

        rdfint_values = []
        for i, frame in enumerate(rdf_values):
            density_B = atom_B_tol_num[i] / vol
            tmp_rdfint_values = []
            for j in range(data_num):
                sum_int = 0.0
                for k in range(j):
                    sum_int += 4*np.pi*density_B*((k + 1)*bin)**2 * frame[k]*bin
                tmp_rdfint_values.append(sum_int)
            rdfint_values.append(tmp_rdfint_values.copy())
            tmp_rdfint_values.clear()


        for i in range(data_num):
            sum = 0.0
            for j in range(len(rdfint_values)):
                sum += rdfint_values[j][i]
            self.rdfIntegral.append(sum / self.frameCount)
        
        for i in range(data_num):
            sum = 0.0
            for j in range(len(rdf_values)):
                sum += rdf_values[j][i]
            self.rdf.append(sum / self.frameCount)
        print ('RDF complete')

    def getADF(self):
        '''
        Calculates the ADF for the treshold atoms in each frame
        '''
        for frame in self.thresholdAtoms:
            # center each frame using centroid
            frame -= centroid(frame)
            for i in range(len(frame) - 1):
                for j in range(i + 1, len(frame)):
                    self.adf.append(angle_eq(frame[i], frame[j]))
        print ('ADF complete')

    def getIdealGeos(self):
        self.idealGeos = geometries.idealGeometries(self.coorNum, self.dist)

    def checkWithUser(self):
        print ('Do the following look OK?')
        while (True):
            print ('\tMax of RDF: '+str('{:.2f}'.format(self.dist)))
            YN_inp = input ('Y/N?: ')
            try:
                if (YN_inp.upper() == 'Y' or YN_inp.upper() == 'YES'):
                    break
                elif (YN_inp.upper() == 'N' or YN_inp.upper() == 'NO'):
                    r_inp = input ('Enter value for the max of the RDF: ')
                    try:
                        self.dist = float(r_inp)
                        break
                    except:
                        print ('Enter a floating point number')
                else:
                    print ('Enter Y or N')
            except:
                print ('Enter Y or N')
        while (True):
            print ('\tThreshold of 1st coordination sphere: '+str('{:.2f}'.format(self.maxR)))
            YN_inp = input ('Y/N?: ')
            try:
                if (YN_inp.upper() == 'Y' or YN_inp.upper() == 'YES'):
                    break
                elif (YN_inp.upper() == 'N' or YN_inp.upper() == 'NO'):
                    r_inp = input ('Enter value for the threshold of the coordination sphere: ')
                    try:
                        self.maxR = float(r_inp)
                        break
                    except:
                        print ('Enter a floating point number')
                else:
                    print ('Enter Y or N')
            except:
                print ('Enter Y or N')
        while (True):
            print ('\tCoordination number: '+str(int(self.coorNum)))
            YN_inp = input ('Y/N?: ')
            try:
                if (YN_inp.upper() == 'Y' or YN_inp.upper() == 'YES'):
                    break
                elif (YN_inp.upper() == 'N' or YN_inp.upper() == 'NO'):
                    r_inp = input ('Enter value for the coordination number: ')
                    try:
                        self.coorNum = int(r_inp)
                        break
                    except:
                        print ('Enter an integer number')
                else:
                    print ('Enter Y or N')
            except:
                print ('Enter Y or N')


def dist_eq(c1, c2):
    return np.sqrt( (c1[0] - c2[0])**2 + (c1[1] - c2[1])**2 + (c1[2] - c2[2])**2 )

def angle_eq(c1, c2):
    # make sure to center coordinates around origin first
    a = dist_eq(c1, [0., 0., 0.])
    b = dist_eq(c2, [0., 0., 0.])
    c = dist_eq(c1, c2)
    return np.arccos( (a**2 + b**2 - c**2) / (2*a*b) ) * (180. / np.pi)

def centroid(arr):
    # a - 3D array of coordinates
    return np.mean(arr, axis=0)

def kabschRMSD(A, B):
    # kabsch algorithm for rotating/translating set of coordinates A onto B
    A = A - centroid(A)
    B = B - centroid(B)
    A = kabschRotate(A, B)
    return rmsd(A, B)

def kabschRotate(A, B):
    U = kabsch(A, B)
    A = np.dot(A, U)
    return A

def kabsch(A, B):
    C = np.dot(np.transpose(A), B)
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    U = np.dot(V, W)
    return U

def rmsd(A, B):
    diff = np.array(A) - np.array(B)
    N = len(A)
    return np.sqrt((diff * diff).sum() / N)

def reorderIdealGeometries(real, ideals):
    minIdeals = ideals
    for i, ideal in enumerate(ideals):
        min = 999.
        for permutation in itertools.permutations(ideal[1]):
            checkPerm = kabschRMSD(real, permutation)
            if (checkPerm < min):
                min = checkPerm
                minIdeals[i][1] = np.array(permutation)
    return minIdeals

