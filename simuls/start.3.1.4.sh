# time internal unit: 10 seconds
# length internal unit: 2 mm
# period: 4 days
# number of hoses: 12
# radius: 4 mm
../simulWaterHoseWall/bin/Release/simulWaterHoseWall -v --lengthU=0.002 --timeU=10 --width=0.094 --deep=2. --period=345600 --earthLambda=0.6 --earthVHeatCap=1.8E6 --waterVHeatCap=4.18E6 --hoseNb=12 --hoseRadius=0.004 --hose1Deep=0.100 --hoseInterv=0.164 res.3.1.4
