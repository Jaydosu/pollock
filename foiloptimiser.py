"""
Using xfoil to optimise a pollock aerofoil. Aerofoil profile to be used for a foil for an ozgoose sailboat.
Pollock aerofoil: https://www.desmos.com/calculator/vmy0t1wg7x


"""
# importing libraries
import os
import numpy as np
import matplotlib.pyplot as plt

def create_pollock_aerofoil(thickness, chordlength, x_le, x_te, S):
  """
  Function to create a pollock aerofoil using the piecewise functions defined in the function. Returns a dat file of the aerofoil.
  Example dat file located at "/joukowsk.dat".

  To be used with xfoil to optimise the aerofoil.

  Parameters:
  thickness: float
    Thickness of the aerofoil.
  chordlength: float
    Chord length of the aerofoil.
  x_le: float
    Leading edge thickness.
  x_te: float
    Trailing edge thickness.
  S: float
    Shape factor.

  Returns:
  dat file of the aerofoil.
  """

  # Defining piecewise functions for the pollock aerofoil.
  y1 = lambda x, x_le: (thickness / 2) * \
                    (((8*np.sqrt(x))/\
                    (3*np.sqrt(x_le*thickness))) - \
                    ((2*x)/(x_le*thickness)) + ((x**2)/\
                    (3*(x_le*thickness)**2)))

  y2 = lambda x, x_te, S: (thickness / 2) * \
                      (1\
                      + ((2 * S) - 3) * \
                      (((x - chordlength + x_te*thickness)**2)/((x_te*thickness)**2))\
                      + (2 * (1 - S)) * \
                      (((x - chordlength + x_te*thickness)**3)/((x_te*thickness)**3)))

  y3 = lambda x: thickness/2 + 0*x

  divisions = 500
  x_values = np.linspace(0, chordlength, 3 * divisions)
  x1 = np.linspace(0, x_te * thickness, divisions)
  x2 = np.linspace(chordlength - x_te*thickness, chordlength, divisions)
  x3 = np.linspace(x_le*thickness, chordlength - x_te*thickness, divisions)

  # Calculate y values using the y1 function
  fx = y1(x1, x_le)
  gx = y2(x2, x_te, S)
  hx = y3(x3)

  y_values = np.concatenate((fx, hx, gx))

  # scale the aerofoil to a unit chord length
  x_values = x_values / chordlength
  y_values = y_values / chordlength

  # Date file Format: 
  # name
  # (six spaces) number of points (float) (six spaces) number of points (float)
  # x1 y1 (float) (float) (7 significant figures)
  # x2 y2 (float) (float) (7 significant figures) etc. for all points
  # single line break
  # x1 y1 (float) (float) (7 significant figures) (bottom side)
  # single line break

  aerofoil_name = f'pollock_{thickness}_{chordlength}_{x_le}_{x_te}_{S}'
  # lets have 50 points on the top and bottom surfaces.
  # lets make sure there are points at 0 and 1, and the bounds of piecewise functions
  # x1 = np.linspace(0, x_te * thickness, divisions)
  # x2 = np.linspace(chordlength - x_te*thickness, chordlength, divisions)
  # x3 = np.linspace(x_le*thickness, chordlength - x_te*thickness, divisions)
  # [0, etc, (x_le*thickness)/chordlength, etc, (chordlength - x_te*thickness)/chordlength, etc, 1]
  # we'll have 3 points between (x_le*thickness)/chordlength and (chordlength - x_te*thickness)/chordlength
  # because they are all the same y value
  points = 50
  xlist = []
  points = points - len(xlist)

  points_div_2 = points // 2 - 4
  v1 = (x_le*thickness)/chordlength
  v2 = (chordlength - x_te*thickness)/chordlength

  iota_1 = v1 / points_div_2 # this is essentially the spacing between the points
  xbtw_1 = np.linspace(0, v1 - iota_1, points_div_2)

  iota_2 = (1 - v2) / points_div_2 # this is essentially the spacing between the points
  xbtw_2 = np.linspace(v2 + iota_2, 1, points - 4 - points_div_2)

  xlist = [v1, v1 + (v2-v1)/3, v1 + 2*(v2-v1)/3, v2]
  xlist = xlist + list(xbtw_1) + list(xbtw_2)
  xlist = sorted(xlist) # just in case
  print(len(xlist))

  # get scaled y values for each point in xlist
  ylist = []
  for x in xlist:
    if x <= (x_le*thickness)/chordlength:
      ylist.append(y1(x * chordlength, x_le) / chordlength)
    elif x >= (chordlength - x_te*thickness)/chordlength:
      ylist.append(y2(x * chordlength, x_te, S) / chordlength)
    else:
      ylist.append(y3(x * chordlength) / chordlength)

  with open(f'{aerofoil_name}.dat', 'w') as f:
    f.write(f'{aerofoil_name}\n\n')
    
    for x in xlist[::-1]:
      f.write(f'{x:.7f} {ylist[xlist.index(x)]:.7f}\n')

    f.write('\n')

    # do it again for the bottom side (y values are simply inverted)
    # where y values are negative, drop the leading 0
    for x in xlist:
      f.write(f'{x:.7f} -{ylist[xlist.index(x)]:.7f}\n')

    return f'{aerofoil_name}.dat'

thickness, chordlength, x_le, x_te, S = 20, 235, 2, 2, 0.25
aerofoil = create_pollock_aerofoil(thickness, chordlength, x_le, x_te, S)

example_aerofoil = 'pollock_20_235_2_2_0.25.dat'

# Run xfoil with the generated aerofoil dat file
def run_xfoil(aerofoil_dat):
  
  xfoilpath = 'xfoil.exe'
  commands = f"""
  LOAD {aerofoil_dat}
  MDES
  FILT
  EXEC

  PANE
  GDES
  CADD




  
  OPER
  ITER 200
  VISC 1000000
  PACC
  {aerofoil_dat.replace('.dat', '.pol')}
  
  ASEQ 0 5 0.1
  PACC
  QUIT
  """
  
  with open('xfoil_input.txt', 'w') as f:
    f.write(commands)
  
  os.system(f'{xfoilpath} < xfoil_input.txt')

# Run xfoil with the generated aerofoil dat file
run_xfoil(example_aerofoil)

