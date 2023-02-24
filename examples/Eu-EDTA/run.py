import rmsd

outputFile = 'Eu-results'
traj = rmsd.Trajectory(ionID='Eu', elements=['O','N'], boxSize=17.5, framesForRMSD=1000, binSize=0.02, startFrame=1, endFrame=10000)
traj.getAtoms('./Eu-EDTA.xyz') 

traj.getIonNum()
if (traj.ionNum > 1): 
    traj.getWhichIon()
traj.getRDF()
traj.getDist()
traj.getMaxR()
traj.printRDF(outputFile)
traj.checkWithUser() 
traj.getThresholdAtoms()
traj.getADF()
traj.printADF(outputFile)
traj.getIdealGeos()
traj.getRMSDs()
traj.printRMSDs(outputFile)
traj.outputIdealGeometries('') 
