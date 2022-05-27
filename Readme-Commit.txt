This is a new version.

Here we read the data from the Raw file only. We do not iterate through all the files. Also, the "Metadata" is taken from the raw file instead from the file names

It works taking Normal rotation and Inverse rotation of the mole. It expects both rotations

We have fixed the ZeroLoc functionality

We have improved the drift correction

We have added a circular shift in the abs, cmd, time data, needed from a bug in FFMM 

We have fixed the feeddown correction

The x,y values are expressed in mm

The Angle is expressed in mrad

The angle is the angle obtained during angle rotation (Bmain only one component) and we have added the level coil angle calibration
