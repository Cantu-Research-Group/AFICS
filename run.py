import rmsd

outputFile = 'TEST2' # enter output file name here (ex. 'La-3')

traj = rmsd.Trajectory(ionID='Lu', elements=['O','N'], boxSize=10., framesForRMSD=500, binSize=0.02, startFrame=100, endFrame=200)
traj.getAtoms('./aimd-ph11-pos-1.xyz') # input .xyz filename here
traj.getIonNum()
if (traj.ionNum > 1): # get user input if there is more than one of the ion
    traj.getWhichIon()
traj.getRDF()
traj.getDist()
traj.getMaxR()
traj.getCN()
traj.printRDF(outputFile)
traj.checkWithUser() # comment this out if you don't want to be prompted about max and threshold values, for automated running, no need for "Y"
traj.getThresholdAtoms()
traj.getADF()
traj.printADF(outputFile)
traj.getIdealGeos()
traj.getRMSDs()
traj.printRMSDs(outputFile)
traj.outputIdealGeometries('') # if there is a subfolder you would like to save to, enter it in the ''
