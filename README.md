# Bead-Immunoassay-Simulation
Simulation of a bead-based 'sandwich' immunassay.

## Requirements
[Optional] If you want to save animation into a mp4 video file, the application 'FFmpeg' must downloaded.
'FFmpeg' can be downloaded in the following link: https://ffmpeg.org/download.html#build-windows
The code can also save as a 'gif' file if indicated.

## 1. Introduction to Bead-Based 'Sandwich' Immunoassay

A 'sandwich' immunoassay is an assay to determine the concentration of a target antigen/analyte.
It is composed of three main components:  
1. solid substrate coated with antibodies (capture antibodies) that bind to the target antigen
2. target antigen                        
3. antibodies (detector antibodies) conjugated with a signal inducing agent (ex. fluorophore)

When all components are mixed together, they form a 'sandwich' complex where both antibodies (1 and 2) bind to the target antigen.
Since the capture anibodies are anchored to a solid substrate, the unbound components can be washed away and only the bound set of particles will remain.
Among the remaining bound particles, only the bound sets that have formed a complete 'sandwich' complex would be able to be read as a signal.
The more antigens you have, more 'sandwich' complexes will form and the signals will be stronger.

A bead-based 'sandwich' immunoasssay is when the solid substrate, the capture antibodies are bound to, are beads. Since beads have mobility in solution, they can move around the solution and bind with the other component.

This code shows a simluation of the components as circular particles that move in certain velocities in a 2 dimentional sapce. Through the simulation we can see how fast or how much the 'sandwich' complex is formed depending on the number of each component as well as the bead size.

## 2. Running the Simulation
### Initial setup of parameters
The following parameters can be set to run the simulation.
Component particles:
Bead, antigen, dAb-fluorophore (dAb)
1. Number of each component particle:
    bead_num, antigen_num, dAb_num
2. Radius of bead particle:
    bead_rad (Antigen and dAb-fluorophore radius is fixed.)
3. Color of each component particle when unbound and bound, in order:
    bead_color, antigen_color, dAb_color
    
Such parameters are organized in a list named 'components'.

Other:
- frame_sz: size of the frame of simulation
- speed_calibration_constant: multiplied to initial velocity of the particles
                                    
If particles are too crowded, errors occur. Preset parameters in the code are optimal for demo.

### Running and saving the simulation animation
The simulation is shown as an animation using the matplotlib.animation.Funcanimation class.
You can specfiy the following parameters in this line:

'''
anim = animation.FuncAnimation(fig2, animate, frames=9000, interval=2, blit=True, repeat=False, save_count=10000
'''

- frames: number of frames the animation will be run for
- interval: time delay between frames in milisecond
- repeat: If True, the animation will run continuously. If False, the animation will stop after the number of frames.

For saving and showing the animation in a figure, you can set the parameters in this line:

'''
save_and_show_animation(anim, save=True, fps=50, filetype='mp4', filename='Bead_Immunoassay_Sim.mp4')
'''

- save: If save=True, the animation will be saved as a video file with the specified file type (filetype) and name (filename). Animation will show on a figure window after the video file is created. This takes a long time depending on the frame number specified in the animation.FuncAnimation function (frames). For saving in mp4 files, ffmpeg must be installed and its file path should be directed through matplotlib.rcParams['animation.ffmpeg_path']. If save=False, the animation will not be saved and show on figure window right away.

                                    


