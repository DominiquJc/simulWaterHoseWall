# time internal unit: 10 seconds
# length internal unit: 2 mm
# period: 9 days
# number of hoses: 6
# radius: 6 mm
../simulWaterHoseWall/bin/Release/simulWaterHoseWall -v --lengthU=0.002 --timeU=10 --width=0.208 --deep=2. --period=777600 --earthLambda=0.6 --earthVHeatCap=1.8E6 --waterVHeatCap=4.18E6 --hoseNb=6 --hoseRadius=0.006 --hose1Deep=0.100 --hoseInterv=0.360 res.3.2.6
