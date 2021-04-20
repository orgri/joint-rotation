fileName = 'test.txt'
fileDirectory = '/home/{}/Dropbox/Common_Folder/Testing/current_results/'  # must exist on system
intelegentStop = True
amplitude = [0.1, 0.2, 0.3, 0.4]
numberOfPeriods = 1000              	# number of periods between energy comparisons for intelegent stop
maxNumberOfPeriods = 1e6            	#
numberOfExtremums = 2000		#
skipPeriods = 0                     	# skip periods before start printing with printStep to log.txt
printStep = 1000.0                  	# print: time=printStep*period/timeStep, bigger value - less print
timeStep = 2000                     	# dt = period / timeStep
acceptableDiff = 1e-5
variable = 'Frequency'              	# Amplitude, Frequency
variableMin = 0.05
variableMax = 1.21
variableStep = 0.025
alpha = 0.05
magnetization = 338
viscosity = 0.0001
smallTheta = 0.1
smallPhi = 0.01
bigTheta = 0.7
bigPhi = 0.2
linearPolarization = False  		# if false then circular polarization: hx*cos(omega*t)+hy*sin(omega*t)
task = 'SimpleRun'  			# SimpleRun, PowerLoss, currentlu use only SimpleRun
numericalMethod = 'RK4'  		# Euler, RK2, RK4
debug = False
frequency = 0.2
outputPrecision = 15
#===================================================================================================

import subprocess

configTemplate = '''
fileName = {}
fileDirectory = {}
intelegentStop = {}
amplitude = {}
numberOfPeriods = {}
maxNumberOfPeriods = {}
numberOfExtremums = {}
skipPeriods = {}
printStep = {}
timeStep = {}
acceptableDiff = {}
variable = {}
variableMin = {}
variableMax = {}
variableStep = {}
alpha = {}
magnetization = {}
viscosity = {}
smallTheta = {}
smallPhi = {}
bigTheta = {}
bigPhi = {}
linearPolarization = {}
task = {}
numericalMethod	= {}
debug = {}
frequency = {}
outputPrecision = {}
'''
firstPCUser = 'vlad'
secondPCUser = 'cuda'
workDir = '/home/{}/Dropbox/Common_Folder/Testing/{}'
sync = ''

def writeConfig(output, file, fileDir, ampl):
    with open(output, 'w') as f:
        f.write(configTemplate.format(fileName, fileDir, str(intelegentStop).lower(),
                            ampl, numberOfPeriods, maxNumberOfPeriods,
                            numberOfExtremums, skipPeriods, printStep,
                            timeStep, acceptableDiff, variable,
                            variableMin, variableMax, variableStep,
                            alpha, magnetization, viscosity,
                            smallTheta, smallPhi, bigTheta, bigPhi,
                            str(linearPolarization).lower(), task, numericalMethod,
                            str(debug).lower(), frequency, outputPrecision))

for h in amplitude:
    fileName = "h={}.txt".format(h)
    idx = amplitude.index(h)
    inputFileName = "input_{}.cfg".format(idx)

    if idx < (len(amplitude) / 2):
        runFile = workDir.format(firstPCUser, 'Joint_Rotation')
        configDir = workDir.format(firstPCUser, inputFileName)
        outputDir = fileDirectory.format(firstPCUser)
        cmd = "{} -f {} &".format(runFile, configDir)

        writeConfig(configDir, fileName, outputDir, h)

        p = subprocess.Popen(cmd, shell=True)			#run Joint_Rotation on first PC (192.168.5.33)

    else:
        runFile = workDir.format(secondPCUser, 'Joint_Rotation')
        configDir = workDir.format(secondPCUser, inputFileName)
        inputHost = workDir.format(firstPCUser, inputFileName)
        outputDir = fileDirectory.format(secondPCUser)
        cmd = "ssh -t cuda@192.168.5.39 -p 777 '{} -f {}'".format(runFile, configDir)
        sync = "rsync -aP -e 'ssh -p 777' {} cuda@192.168.5.39:{}".format(inputHost, configDir)

        writeConfig(inputHost, fileName, outputDir, h)

        p = subprocess.Popen(sync, shell=True)			# sync config file on second PC (192.168.5.39)
        p.wait()

        p = subprocess.Popen(cmd, shell=True)			#run Joint_Rotation on second PC (192.168.5.39) 

#p=subprocess.Popen(cmd, shell=True)
#ssh -t cuda@192.168.5.39 -p 777 'cat .ssh/authorized_keys &'
#rsync -aP -e 'ssh -p 777' ./input.cfg  cuda@192.168.5.39:Dropbox/Common_Folder/Testing/input.cfg


def printFile(file, ampl):
    print(configTemplate.format(fileName, fileDirectory, str(intelegentStop).lower(),
                            ampl, numberOfPeriods, maxNumberOfPeriods,
                            numberOfExtremums, skipPeriods, printStep,
                            timeStep, acceptableDiff, variable,
                            variableMin, variableMax, variableStep,
                            alpha, magnetization, viscosity,
                            smallTheta, smallPhi, bigTheta, bigPhi,
                            str(linearPolarization).lower(), task, numericalMethod,
                            str(debug).lower(), frequency, outputPrecision))
