# Software Carpentry Final Project - Bead Immunoassay Simulation
# Name: Jeongyun Kim (ID:jkim667)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.patches import Rectangle
from matplotlib import animation
from itertools import combinations
import matplotlib
matplotlib.rcParams['animation.ffmpeg_path'] = r'C:\Users\Jeongyun\Desktop\ffmpeg\bin\ffmpeg.exe'


def overlap_test(p1, p2):
    """
    A function to see if a particle and another particle overlaps.

    Parameters
    ----------
    p1 : patches.Circle
        Investigating particle.
    p2 : patches.Circle
        Comparing particle.

    Returns
    -------
    Boolean
    True : If particles overlap.
    False: If particles do not overlap.

    """
    return np.hypot(*(p1.center - p2.center)) < p1.radius + p2.radius


def init_particle_setup():
    """
    Initializes the properties of the particles.
    Postition of particles are randomly selected without overlapping with each
    other.
    Velocity of particles are set to have a fixed amplitude that is
    porportional to the inverse of the particle area and a randomly selected
    direction.
    A status vector is set to reflect the type of particle, color and bound
    status.

    Returns
    -------
    particles : list of matplotlib.patches.Circle
        Contains the center coordinates, radius and color of the particle.
        Can be drawn as pathces in a axes.
    velocities : list of lists
        Velocity of the particles are stored.
    status : list of lists
        For each particle, the following status of the particle is stored:
            [particle type, color, bound status]
            particle type: bead (0), antigen (1) and dAb-fluorophore (2)
            color: RGB code for current bound status
            bound status:
                bead: index of all bound antigens
                antigen: index of bound bead or dAb-fluorophore
                dAb-fluorophore: index of bound antigen and bead

    """
    fig, ax = plt.subplots(figsize=(10,10))

    ax.set_xlim(0,frame_sz)
    ax.set_ylim(0,frame_sz)

    particles = []
    velocities = []
    status = []

    for i in range(len(components)):
        num = components[i][0]
        rad = components[i][1]
        color = components[i][2]

        n = 0
        while n < num:
            center = np.random.randint(rad, frame_sz-rad, 2)
            curr_particle = Circle(center, rad, color=color[0])
            if particles == []:
                particles.append(curr_particle)
                status.append([i, color[0], [len(particles)-1]])
                n += 1
            else:
                for particle in particles:
                    if overlap_test(curr_particle, particle):
                        break
                else:
                    particles.append(curr_particle)
                    status.append([i, color[0], [len(particles)-1]])
                    n += 1

    for particle in particles:
        ax.add_patch(particle)
        velocity_amp = speed_calibration_constant*1/rad**2
        velocity_ang = 2*np.pi*np.random.random()
        velocities.append(
            [velocity_amp*np.cos(velocity_ang), velocity_amp*np.sin(velocity_ang)]
            )
    plt.show()    
    return particles, velocities, status


def next_step():
    """
    Updating particle properties after next step of the particles.
    particle_collision():
        Velocity update and position adjustment after all collisions.
    boundary_collision():
        Velocity update and position adjustment after colliding to the frame
        boundary.

    Returns
    -------
    particles : list of matplotlib.patches.Circle
        Updated next step of particles.

    """
    particle_collision()

    for n in range(len(particles)):  # position update for next step

        boundary_collision(particles[n].center, n)
        n_center = particles[n].center+np.array(velocities[n])
        n_particle = Circle(n_center, Circle.get_radius(particles[n]), color=status[n][1])
        particles[n] = n_particle

    return particles


def complex_count_update():
    """
    Count all 'sandwich' complex formed.
    When 'sandwich' complex is formed, the length of the dAb_fluorophore bound
    status will be 3.

    Returns
    -------
    complex_count : int
        The number of 'sandwich' complex formed.

    """
    complex_count = 0
    dAb_index = list(range(int(components[0][0]+components[1][0]), len(particles)))
    for i in dAb_index:
        if len(status[i][2]) == 3:
            complex_count += 1
    return complex_count


def particle_collision():
    """
    Handles particle pairs that collide to each other.
    Depending on the type of particles and their status, different collision
    updates are made.
    
    """
    particle_pairs = combinations(range(len(particles)), 2)

    for i, j in particle_pairs:
        if overlap_test(particles[i], particles[j]):

            overlap_distance_modification(i, j)

            p1_type, p2_type = status[i][0], status[j][0]

            # If particles are same type of components, elastic collision
            if (p1_type == p2_type) or (p1_type == 0 and p2_type == 2) or (p1_type == 2 and p2_type == 0):
                particle_elastic_collision(i, j)

            # If particle pair is bead and antigen
            if p1_type + p2_type == 1:
                particle_binding1(i, j)

            # If particle pair is antigen and dAb-fluorophore
            if p1_type + p2_type == 3:
                particle_binding2(i, j)
            
            # overlap_distance_modification(i, j)

                

def complete_bound_set(i):
    """
    Set of particles that are bound together and where particle i is included.

    Parameters
    ----------
    i : int
        Particle index.

    Returns
    -------
    list
    
    """    
    
    bound_set = []
    for item1 in status[i][2]:
        for item2 in status[item1][2]:
            for item3 in status[item2][2]:
                for item4 in status[item3][2]:
                    bound_set.append(item1)
                    bound_set.append(item2)
                    bound_set.append(item3)
                    bound_set.append(item4)
    return [*set(bound_set)] # Remove any overlapping index


def overlap_distance_modification(i, j):
    """
    For particles that overlap, the overlapping distance is modified so that
    the particles are contacting at each others border.
    This fixes the error when the particles are trapped infinitely within their
    particle border.

    Parameters
    ----------
    i : int
        Particle index i
    j : int
        Particle index j

    """
    distance = np.linalg.norm(particles[i].center-particles[j].center) # distance between two particle centers
    overlap_distance = particles[i].radius + particles[j].radius - distance # distance overlapping
    
    m_distance = overlap_distance/distance*(particles[j].center-particles[i].center) # modified distance change
     
    bound_set = complete_bound_set(j)

    for item in bound_set:
        if not item == i:
            particles[item].center = particles[item].center+m_distance


def total_mass(i):
    
    total_mass = 0
    for i in status[i][2]:
        total_mass += np.pi*(particles[i].radius)**2
    return total_mass


def elastic_collision_update(p1, p2, p1_mass, p2_mass, p1_vel, p2_vel):
    """
    Update velocity of particles that are to bounce off of each other.
    (Elastic collision)
    
    Parameters
    ----------
    i : int
        Particle index i
    j : int
        Particle index j

    """
    
    n_p1_vel = (p1_mass-p2_mass)/(p1_mass+p2_mass)*p1_vel+2*p2_mass/(p1_mass+p2_mass)*p2_vel
    n_p2_vel = 2*p1_mass/(p1_mass+p2_mass)*p1_vel-(p1_mass-p2_mass)/(p1_mass+p2_mass)*p2_vel

    velocities[p1] = n_p1_vel
    velocities[p2] = n_p2_vel


def particle_binding_update(p1, p2, p1_mass, p2_mass, p1_vel, p2_vel):
    """
    Update velocity of particles that are to be bound to each other. 
    (Inelastic collision)
    
    Parameters
    ----------
    i : int
        Particle index i
    j : int
        Particle index j

    """
    
    n_vel = (p1_mass*p1_vel+p2_mass*p2_vel)/(p1_mass+p2_mass)

    velocities[p1] = n_vel
    velocities[p2] = n_vel
    
    # Update bound status in 'status' vector.
    if p2 not in status[p1][2]:
        status[p1][2].append(p2)
    if p1 not in status[p2][2]:
        status[p2][2].append(p1)

    # Change colors when bound to another particle.
    status[p1][1] = components[status[p1][0]][2][1]
    status[p2][1] = components[status[p2][0]][2][1]
    

def syncronize_bound_particles(i, j):
    """
    Update velocity of all particles that are bound to each particle (i and j)
    to the same velocity as each particle.
    
    Parameters
    ----------
    i : int
        Particle index i
    j : int
        Particle index j

    """
    
    i_bound_set = complete_bound_set(i)
    for index1 in i_bound_set:
        velocities[index1] = velocities[i]
    
    j_bound_set = complete_bound_set(j)
    for index1 in j_bound_set:
        velocities[index1] = velocities[j]


def particle_elastic_collision(i, j):
    """
    Collision update When particles are the same component type or if the pair
    is bead and dAb-fluorophore.
    Velocity is updated through elastic collision.
    
    All particles that are bound to each particle is also updated to the same
    velocity to move together at the next step.

    Parameters
    ----------
    i : int
        Particle index i
    j : int
        Particle index j

    """
    p1_mass, p2_mass = total_mass(i), total_mass(j)
    p1_vel, p2_vel = np.array(velocities[i]), np.array(velocities[j])
    
    elastic_collision_update(i, j, p1_mass, p2_mass, p1_vel, p2_vel)
                
    syncronize_bound_particles(i, j)


def particle_binding1(i, j):
    """
    Collision update when particles are bead (i) and antigen (j).
    If antigen is not bound to anything ot just bound to a dAb-fluorophore,
    binding occurs and velocity is updated as inelastic collision. Otherwise,
    binding does not occur and velocity is updated as elastic collision.

    All particles that are bound each particle is also updated to the same
    velocity to move together at the next step.

    Parameters
    ----------
    i : int
        Particle index i
    j : int
        Particle index j

    """  
    p1_mass, p2_mass = total_mass(i), total_mass(j)
    p1_vel, p2_vel = np.array(velocities[i]), np.array(velocities[j])

    dAb_index = list(range(int(components[0][0]+components[1][0]), len(particles)))
    # If antigen (j) is not bound to anything or just bound to dAb-fluorophore
    if (len(status[j][2]) == 1) or ((len(status[j][2]) == 2) and (status[j][2][1] in dAb_index)):       
        particle_binding_update(i, j, p1_mass, p2_mass, p1_vel, p2_vel)
    else:        
        elastic_collision_update(i, j, p1_mass, p2_mass, p1_vel, p2_vel)
        
    syncronize_bound_particles(i, j)

        
def particle_binding2(i,j): # Antigen (i) and dAb-Fluorophore (j) binding
    """
    Collision update when particles are antigen (i) and dAb-fluorophore (j).
    If the antigen is not bound to anything or just bound to a bead, and the
    dAb-fluorophore is not bound to anything binding occurs and velocity is
    updated as inelastic collision. Otherwise, binding does not occur and
    velocity is updated as elastic collision.
    
    All particles that are bound each particle is also updated to the same
    velocity to move together at the next step.

    Parameters
    ----------
    i : int
        Particle index i
    j : int
        Particle index j

    """    
    p1_mass, p2_mass = total_mass(i), total_mass(j)
    p1_vel, p2_vel = np.array(velocities[i]), np.array(velocities[j])

    bead_index = list(range(components[0][0]))
    # If antigen is not bound to anything or just bound to bead
    if (len(status[i][2]) == 1) or (status[i][2][1] in bead_index):
        if (len(status[j][2])) == 1:        
            particle_binding_update(i, j, p1_mass, p2_mass, p1_vel, p2_vel)
        else:            
            elastic_collision_update(i, j, p1_mass, p2_mass, p1_vel, p2_vel)
    else:        
        elastic_collision_update(i, j, p1_mass, p2_mass, p1_vel, p2_vel)
    
    syncronize_bound_particles(i, j)

    # Update bead if it is bound to the antigen the dAb is bound to.
    if len(status[j][2]) == 2:
        for index3 in bead_index:
            if status[j][2][1] in status[index3][2]:
                status[j][2].append(index3)


def boundary_collision(n_center,n):
    """
    When particle is out of the bounds of the frame, velocity is updated to
    bounce off the wall and the center is updated so the the particle is within
    the borders of the frame.

    All particles that are bound each particle is also updated to the same
    velocity to move together at the next step.

    Parameters
    ----------
    n_center :
        Corrdinate for particle center.
    n : int
        Particle index.

    """
    if n_center[0] - particles[n].radius < 0:
        n_center[0] = particles[n].radius
        velocities[n][0] = -velocities[n][0]
    if n_center[0] + particles[n].radius > frame_sz:
        n_center[0] = frame_sz-particles[n].radius
        velocities[n][0] = -velocities[n][0]
    if n_center[1] - particles[n].radius < 0:
        n_center[1] = particles[n].radius
        velocities[n][1] = -velocities[n][1]
    if n_center[1] + particles[n].radius > frame_sz:
        n_center[1] = frame_sz-particles[n].radius
        velocities[n][1] = -velocities[n][1]

    bound_set = complete_bound_set(n)
    for index in bound_set:  # Need to update
        velocities[index] = velocities[n]


def animate(i):
    """
    Animation function for animation generation using
    matplotlib.animation.FuncAnimation.

    """    
    white_square = Rectangle((0,0),frame_sz,frame_sz,color=(1,1,1))
    ax2.add_patch(white_square)    
    next_step()
    for particle in particles:
       ax2.add_patch(particle)

    return particles

def save_and_show_animation(anim, save, fps, filetype, filename):
    """
    Animation saving and showing function for generated animation.

    """
    if save:
        if filetype == 'mp4':
            writer = animation.FFMpegWriter(fps)
        anim.save(filename)
    
    plt.show()
    
  
if __name__ == '__main__':
    '''
    ### Parameters for simualtion ###
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
        frame_sz: size of the frame of simulation
        speed_calibration_constant: multiplied to initial velocity of the
                                    particles
                                    
    *If particles are too crowded, errors occur.
    Preset parameters below are optimal for demo.
        
    '''
    bead_num = 1
    bead_rad = 75
    bead_color = [(0.45,0.65,0.78), (0.19,0.35,0.52)]

    antigen_num = 20
    antigen_color = [(0.99,0.56,0.20), (0.87,0.12,0.10)]

    dAb_num = 40
    dAb_color = [(0.5,0.5,0.5), (0.03,1,0.03)]

    components = [
        [bead_num, bead_rad, bead_color],
        [antigen_num, 2.5, antigen_color],
        [dAb_num, 7, dAb_color]
        ]   

    frame_sz = 400

    speed_calibration_constant = 10
    
    # Initialize particle setup
    particles, velocities, status = init_particle_setup()
    
    # Initialize figure window, the animation will show on.
    fig2, ax2 = plt.subplots(figsize=(10,10))
    ax2.set_xlim(0,frame_sz)
    ax2.set_ylim(0,frame_sz)
    
    # Run animation (update particles)
    anim = animation.FuncAnimation(
        fig2, animate, frames=200, interval=2, blit=True, repeat=True, save_count=200
        )
    
    '''
    ### Save and/or show animation ###
    If save=True, the animation will be saved as a video file with the
    specified file type (filetype) and name (filename). Animation will show on
    a figure window after the video file is created. This takes a long time
    depending on the frame number specified in the animation.FuncAnimation
    function (frames).
    
    For saving in mp4 files, ffmpeg must be installed and its file path should
    be directed through matplotlib.rcParams['animation.ffmpeg_path'].
    
    If save=False, the animation will not be saved and show on figure window
    right away.
    
    '''
    save_and_show_animation(
        anim, save=False, fps=5, filetype='mp4', filename='Bead_immunoassay.mp4'
        )
    
    # Counting the number of 'sandwich' complexes formed after simulation.
    complex_count = complex_count_update()
    




