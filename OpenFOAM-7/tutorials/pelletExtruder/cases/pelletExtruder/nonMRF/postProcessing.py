import pandas as pd

#Convert forces to CSV
read_file = pd.read_csv (r'postProcessing/forcesRotor/0/forces2.tmp')
read_file.to_csv (r'postProcessing/forcesRotor/0/forces.csv', index=None)

#Convert barrel shear stress to CSV
read_file = pd.read_csv (r'postProcessing/barrelShearStress/0/wallShearStress2.tmp')
read_file.to_csv (r'postProcessing/barrelShearStress/0/wallShearStress2.csv', index=None)

#Convert screw shear stress to CSV
read_file = pd.read_csv (r'postProcessing/screwShearStress/0/wallShearStress2.tmp')
read_file.to_csv (r'postProcessing/screwShearStress/0/wallShearStress2.csv', index=None)

#Convert delta_p stress to CSV
read_file = pd.read_csv (r'postProcessing/delta_p/0/fieldValueDelta.tmp')
read_file.to_csv (r'postProcessing/delta_p/0/fieldValueDelta.csv', index=None)

#Convert delta_p.region1 to CSV
read_file = pd.read_csv (r'postProcessing/delta_p.region1/0/surfaceFieldValue.tmp')
read_file.to_csv (r'postProcessing/delta_p.region1/0/surfaceFieldValue.csv', index=None)

#Convert delta_p.region2 to CSV
read_file = pd.read_csv (r'postProcessing/delta_p.region2/0/surfaceFieldValue.tmp')
read_file.to_csv (r'postProcessing/delta_p.region2/0/surfaceFieldValue.csv', index=None)

#Convert massFlowRateInlet to CSV
read_file = pd.read_csv (r'postProcessing/massFlowRateInlet/0/surfaceFieldValue.tmp')
read_file.to_csv (r'postProcessing/massFlowRateInlet/0/surfaceFieldValue.csv', index=None)

#Convert massFlowRateOutlet to CSV
read_file = pd.read_csv (r'postProcessing/massFlowRateOutlet/0/surfaceFieldValue.tmp')
read_file.to_csv (r'postProcessing/massFlowRateOutlet/0/surfaceFieldValue.csv', index=None)

#Convert minMax to CSV
read_file = pd.read_csv (r'postProcessing/minMax/0/volFieldValue.tmp')
read_file.to_csv (r'postProcessing/minMax/0/volFieldValue.csv', index=None)

#Convert wallHeatFlux to CSV
read_file = pd.read_csv (r'postProcessing/wallHeatFlux/0/wallHeatFlux.tmp')
read_file.to_csv (r'postProcessing/wallHeatFlux/0/wallHeatFlux.csv', index=None)
