# What is it?
This program performs a thermal simulation of a thick wall, several time crossed by hoses where water circulates in a closed loop.  
It provides simulation results, in amplitude and phase shift, including flux density (thermal power by surface), and temperature at various locations.  
What is expected from water circulation: increase the heat flow with the room in response to a temperature change. Or else, increase the apparent thermal inertia.  

It is a little program, unpretentious, which may be used as source, for example to add temperature monitoring points. It is developed in C language with codeblocks.

# Who is it for?
People interested in the thermal aspect in alternative constructions (such as earthships).

# What is its goal?
This program allows the **study** of the thermal behavior of this device, for various periods, and various wall and hoses configurations.

# How use it?
You can install it locally.  
But you can also keep its source form.

# Use as source
- download simulWaterHoseWall-X.x.tgz
- expand simulWaterHoseWall-X.x.tgz
- open simulWaterHoseWall-X.x in your browser
- if codeblocks is installed, click on the cbp file
- else: ./configure, make, and use your development tools

# Local installation
- download simulWaterHoseWall-X.x.tgz
- expand simulWaterHoseWall-X.x.tgz
- cd simulWaterHoseWall-X.x
- ./configure --prefix=$HOME
- make
- make install

# Use
- create a work directory, for example "simuls"
- cd this directory
- create lauching shells based on the files in examples directory
- give execution rights to these launchers
- to lauch, type "./*launcher_name*" in the console.
- otherwise, without execution rights, type "sh *launcher_name*" in the console.

