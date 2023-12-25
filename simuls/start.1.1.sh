# time internal unit: 1 second
# length internal unit: 1 mm
# period: 1 day
# number of hoses: 5
# radius: 6 mm
../simulWaterHoseWall/bin/Release/simulWaterHoseWall -v --lengthU=0.001 --timeU=1 --width=0.250 --deep=2. --period=86400 --earthLambda=0.6 --earthVHeatCap=1.8E6 --waterVHeatCap=4.18E6 --hoseNb=5 --hoseRadius=0.006 --hose1Deep=0.100 --hoseInterv=0.433 res.1.1
