import numpy as np
from astropy.table import Table
from jlu.astrometry import Transform2D as T2D


def calc_pos():

    #star positions in arcminutes (RA, Dec) from Zemax 
    #order is same as Toomey spreadsheet, WFS3, WFS2, WFS1, WFS0, star 5 -> star 5 position is commented out
    star_ra = np.array([-5.2866, 6.3783, -6.36, 5.4133]) # 3.01166])
    star_dec = np.array([-5.12833, -1.13822, 4.73166, 6.856667]) # -7.79833])

    #want coordinates in polar form -> calulate scale change from radial part, calculate theta from offset
    #this assumes field is correctly centered on the pointing (between stars and cal plate)

    #ferrule x/y positions (inches) based on manufactured specifications
    fer_x =  np.array([-1.384, 1.667, -1.665, 1.414])#, .783])
    fer_y = np.array([-1.342, -.298, 1.238, 1.793])#, -2.039])

    #offsets between star center and ferrule center in arcseconds

    #nov 18th data
    #fer_off_ra = np.array([10, -1, -1, -8])
    #fer_off_dec = np.array([10, 16, 1, 8])

    #nov 14th data
    #fer_off_ra =  np.array([0.5, -5.5, 10.5, 0.5])
    #fer_off_dec = np.array([5, -6, -3, -12.5])

    #Jan 9th data WFS3, WFS 2, WFS 1, WFS0
    #fer_off_ra = np.array([3.0, 1., -1. , -2.5])
    #fer_off_dec = np.array([-14, -9., -14.0, -8.5])

    #Jan 9th data after rotation
    #fer_off_ra = np.array([3.5,0.5,0.75,-1])
    #fer_off_dec = np.array([-2.5,2,-2.5,2.5])

    #Jan 9th data after second rotation
    fer_off_ra = np.array([0.5,-.5,2.,0.])
    fer_off_dec = np.array([-2,-2.,-1.0,-.5])
    '''
    need to fit rotation/scale between the stars and the ferrules....  Want to know 2 things
    1) what roation angle should the instrument be put at
    2)what scale should be used in the generation of the cal plate
    -this requires fitting for the linear offsets as well
    '''

    #calulate corrected star positions (in arcseconds)
    #possible sign ambiguity (subtract or add offsets?)
    star_ra_cor = star_ra * 60.0 - fer_off_ra
    star_dec_cor = star_dec * 60.0 - fer_off_dec

    #now I have corrected stellar positions and the machined ferrule positions
    #Want to calculate the scale, rotation and offset between the two data sets

    trans = T2D.four_param(fer_x, fer_y, star_ra_cor, star_dec_cor)

    print('scale is '+str(trans.params['scale'] / 25.4)+' "/mm') #convert scale from "/inch to "/mm
    print('rotation is '+str(trans.params['angle'])+ ' degrees')
    print('x offset is '+str(trans.params['transX'])+ ' arcseconds')
    print('y offset is '+str(trans.params['transY'] )+ ' arcseconds')
    #import pdb;pdb.set_trace() 
