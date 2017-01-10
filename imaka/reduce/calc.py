import numpy as np
from scipy import ndimage
import math

def calc_rotation():
    sep = 15.53 # mm in front focal plane

    imaka_magnify = 1.33

    sci_scale = 7.07 # asec / mmm

    sep_in_asec = imaka_magnify * sep * sci_scale

    # Calculate offsets from rotation.
    pos_orig = np.array([sep_in_asec, 0.0])

    angle = np.radians(1.0)

    rot_matrix = np.array([[np.cos(angle), -np.sin(angle)],
                           [np.sin(angle), np.cos(angle)]])

    pos_rot = np.dot(rot_matrix,  pos_orig)

    delta_vec = pos_rot - pos_orig

    print("Original Pos:       E = {0:6.2f}   N = {1:6.2f}".format(pos_orig[0], pos_orig[1]))
    print("Rotated Pos:        E = {0:6.2f}   N = {1:6.2f}".format(pos_rot[0], pos_rot[1]))
    print("Delta (Rot - Orig): E = {0:6.2f}   N = {1:6.2f}".format(delta_vec[0], delta_vec[1]))


def get_pa_scale():
    # e_pos = np.array([0.0, 0.0])
    # n_pos = np.array([10.0, -10.0])
    # x_pos = np.array([614.0, 952.0])
    # y_pos = np.array([4723.0, 5377.])

    # e_pos = np.array([6.0, -3.0])
    # n_pos = np.array([0.0, 0.0])
    # x_pos = np.array([381.4, 172.5])
    # y_pos = np.array([329.0, 374.6])

    # e_pos = np.array([6.0, -3.0])
    # n_pos = np.array([0.0, 0.0])
    # x_pos = np.array([381.4, 172.5])
    # y_pos = np.array([329.0, 374.6])

    # E-position of the science camera (2017-01-09 HST)
    # e_pos = np.array([15.0, -30.0])  # OR 0 60, OR 15 0
    # n_pos = np.array([0.0, 0.0])
    # x_pos = np.array([1770.0, 1650.0])
    # y_pos = np.array([2023.0, 1815.0])

    # E-position of the science camera (2017-01-09 HST)
    e_pos = np.array([0.0, -28.0])  # OR 0 60, OR 15 0
    n_pos = np.array([0.0, 0.0])
    x_pos = np.array([1495.0, 1382.0])
    y_pos = np.array([1985.0, 1796.0])
    
    scale, angle = two_point_solve(e_pos[0], n_pos[0], x_pos[0], y_pos[0],
                                   e_pos[1], n_pos[1], x_pos[1], y_pos[1])

    return


# Want center at (fiber center):
# x = 1660
# y = 1745
# Star coming into:
# x = 1522 (1525)
# y = 2027 
# Star coming into:
# x = 1482
# y = 1965
# After AO -28 25
# x = 1520
# y = 1688
# 

def two_point_solve(e1, n1, x1, y1, e2, n2, x2, y2):
    """
    Calculate the scale and angle from two measured points.
    dE = scale * (dX cos a - dY sin a)
    dN = scale * (dX sin a + dY cos a)

    scale = sqrt( (dE**2 + dN**2) / (dX**2 + dY**2) )
    angle = arctan( (dN * dX - dE * dY) / (dE * dX + dN * dY) )
    """
    de = e1 - e2
    dn = n1 - n2
    dx = x1 - x2
    dy = y1 - y2
    
    scale = ((de**2 + dn**2) / (dx**2 + dy**2))**0.5

    angle = np.arctan2( dn*dx - de*dy, de*dx + dn*dy )

    if hasattr(scale, '__iter__'):
        print('Scale = {0:.4f} arcsec / pixel'.format(scale[0]))
        print('Angle = {0:.3f} degrees'.format(np.degrees(angle[0])))
    else:
        print('Scale = {0:.4f} arcsec / pixel'.format(scale))
        print('Angle = {0:.3f} degrees'.format(np.degrees(angle)))

    return scale, angle

# These are numbers reported from scale_from_fiber_positions()
# camera_angle = {'N':76.33, 'E':-15.68}

# These are numbers reported from get_pa_scale() on 2017-01-17
camera_angle = {'N':76.33, 'E': -60.0}
# camera_angle = {'N':76.33, 'E': -140.0}

def xy_to_en(x, y, camera_pa='E'):
    """
    Calculate the East North offsets for a desired X Y offset.

    Note: x and y should be CURRENT - DESIRED
    """
    scale = 0.127 # arcsec / pixel
    angle = camera_angle[camera_pa] # degrees

    cosa = np.cos(np.radians(angle))
    sina = np.sin(np.radians(angle))

    east = scale * (x * cosa - y * sina)
    north = -scale * (x * sina + y * cosa)

    return east, north
    

def en_to_xy(e, n, camera_pa='E'):
    """
    Calculate the East North offsets for a desired X Y offset.
    """
    scale = 0.127 # arcsec / pixel
    angle = camera_angle[camera_pa] # degrees

    cosa = np.cos(np.radians(angle))
    sina = np.sin(np.radians(angle))

    x = (1.0 / scale) * (e * cosa - n * sina)
    y = (1.0 / -scale) * (e * sina + n * cosa)

    return x, y
    

    
# Designed focal plane scale 0.1"  / 15. micron

def scale_from_fiber_positions():
    """
    Analysis of the fiber positions
    """
    

    # Camera Position N - 2016-11-15
    # fiber_en_1 = np.array([0.0, 0.0])
    # fiber_en_2 = np.array([0.0, 12.1858])  # (supposedly N offset)
    # fiber_en_3 = np.array([-15.53, 0.0])  # (supposedly E offset... note the sign flip)
    # The sign flip is necessary to make the two coordinate have the same "handed-ness".
    
    # fiber_xy_1 = np.array([1788.502, 1930.594])
    # fiber_xy_2 = np.array([4211.805, 2553.25])
    # fiber_xy_3 = np.array([917.762, 5300.377])

    # Position E - 2016-11-17
    fiber_en_1 = np.array([0.0, 0.0])
    fiber_en_2 = np.array([0.0, 12.1858])  # (supposedly N offset)
    fiber_en_3 = np.array([15.53, 0.0])  # (supposedly E offset... note the sign flip)
    # The sign flip is necessary to make the two coordinate have the same "handed-ness".
    
    fiber_xy_1 = np.array([1915.28, 2197.8])
    fiber_xy_2 = np.array([1290.2, 4569.8])
    fiber_xy_3 = np.array([5285.4, 3077.8])
    

    pixel_size = 0.006 # 6 micron pixels.

    # Convert the fiber xy positions into millimeters also. The resulting scale will
    # give us 1 / magnification.
    fiber_xy_1 *= pixel_size
    fiber_xy_2 *= pixel_size
    fiber_xy_3 *= pixel_size

    scale_12, angle_12 = two_point_solve(fiber_en_1[0], fiber_en_1[1], fiber_xy_1[0], fiber_xy_1[1],
                                         fiber_en_2[0], fiber_en_2[1], fiber_xy_2[0], fiber_xy_2[1])
    scale_13, angle_13 = two_point_solve(fiber_en_1[0], fiber_en_1[1], fiber_xy_1[0], fiber_xy_1[1],
                                         fiber_en_3[0], fiber_en_3[1], fiber_xy_3[0], fiber_xy_3[1])
    scale_23, angle_23 = two_point_solve(fiber_en_2[0], fiber_en_2[1], fiber_xy_2[0], fiber_xy_2[1],
                                         fiber_en_3[0], fiber_en_3[1], fiber_xy_3[0], fiber_xy_3[1])

    mag_12 = 1.0 / scale_12
    mag_13 = 1.0 / scale_13
    mag_23 = 1.0 / scale_23

    ang_12 = math.degrees(angle_12)
    ang_13 = math.degrees(angle_13)
    ang_23 = math.degrees(angle_23)

    # Mean values, errors
    mag = np.array([mag_12, mag_13, mag_23])
    ang = np.array([ang_12, ang_13, ang_23])
    pos = ['12', '13', '23']

    print('')
    print('Magnification  = {0:.2f} +/- {1:.2f}'.format(mag.mean(), mag.std()))
    print('Position Angle = {0:.2f} +/- {1:.2f} deg'.format(ang.mean(), ang.std()))
    print('')
    print('Individual Measurements:')
    print(' {0:6s}  {1:6s}  {2:6s}'.format('Pos', '  Mag', 'PA_deg'))
    for ii in range(len(mag)):
        print(' {0:6s}  {1:6.2f}  {2:6.1f}'.format(pos[ii], mag[ii], ang[ii]))
    
    
