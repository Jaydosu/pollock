"""
Using xfoil to optimise a pollock aerofoil. Aerofoil profile to be used for a foil for an ozgoose sailboat.
Pollock aerofoil: https://www.desmos.com/calculator/vmy0t1wg7x


"""
# importing libraries
import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, differential_evolution

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

  aerofoil_name = f'p_{x_le}_{x_te}_{S}'
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

  xlist = [v1, v1 + (v2-v1)/3, v1 + 2*(v2-v1)/3, v2] # add the bounds
  xlist = xlist + list(xbtw_1) + list(xbtw_2)
  xlist = sorted(xlist) # just in case

  # get scaled y values for each point in xlist
  ylist = []
  for x in xlist:
    if x <= (x_le*thickness)/chordlength:
      ylist.append(y1(x * chordlength, x_le) / chordlength)
    elif x >= (chordlength - x_te*thickness)/chordlength:
      ylist.append(y2(x * chordlength, x_te, S) / chordlength)
    else:
      ylist.append(y3(x * chordlength) / chordlength)

  with open(os.getcwd()+f'\dat\{aerofoil_name}.dat', 'w') as f:
    f.write(f'{aerofoil_name}\n\n')
    
    for x in xlist[::-1]:
      f.write(f'{x:.7f} {ylist[xlist.index(x)]:.7f}\n')

    f.write('\n')

    # do it again for the bottom side (y values are simply inverted)
    # where y values are negative, drop the leading 0
    for x in xlist:
      f.write(f'{x:.7f} -{ylist[xlist.index(x)]:.7f}\n')

    return f'{aerofoil_name}.dat'

def run_xfoil(aerofoil_dat):
  """
  Function to run xfoil with the generated aerofoil dat file. Returns the cl_max and aoa_max of the aerofoil.

  Parameters:
  aerofoil_dat: str
    Dat file of the aerofoil.
  """
  # ensure that there is no polar file for the aerofoil
  if os.path.exists(os.getcwd() + '\pol\\' + aerofoil_dat.replace('.dat', '.pol')):
    os.remove(os.getcwd() + '\pol\\' + aerofoil_dat.replace('.dat', '.pol'))

  xfoilpath = 'xfoil.exe'
  commands = f"""
  LOAD dat\{aerofoil_dat}

  MDES
  FILT
  EXEC

  PANE
  GDES
  CADD




  
  OPER
  ITER 200
  VPAR
  N 13
  
  VISC 1439701
  PACC
  pol\{aerofoil_dat.replace('.dat', '.pol')}
  
  ASEQ 0 6 0.1
  PACC
  QUIT
  """
  
  with open('xfoil_input.txt', 'w') as f:
    f.write(commands)
  
  with open(os.devnull, 'w') as devnull:
    # Set up the startup info to prevent a console window from appearing
    startupinfo = subprocess.STARTUPINFO()
    startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
    startupinfo.wShowWindow = subprocess.SW_HIDE
    CREATE_NO_WINDOW = 0x08000000

    subprocess.run([xfoilpath], shell=True, stdin=open('xfoil_input.txt'), stdout=devnull, stderr=devnull, startupinfo=startupinfo, creationflags=subprocess.CREATE_NO_WINDOW)

def read_pol(aerofoil_pol):
  with open(os.getcwd() + '\pol\\' + aerofoil_pol, 'r') as f:
    lines = f.readlines()
    
    # find max cl
    cl_max = 0
    aoa_max = 0
    for line in lines[13:]:
      cl = float(line.split()[1])
      aoa = float(line.split()[0])
      if cl > cl_max:
        cl_max = cl
        aoa_max = aoa

  return cl_max, aoa_max

def objective_function(param):
  # Clean up params to remove floating point errors
  param = [round(p, 3) for p in param]

  # Generate the aerofoil data file for the current parameter
  aerofoil_dat = create_pollock_aerofoil(thickness, chordlength, x_le, param[0], param[1])

  print("Current parameters: ", param)
  print("Current aerofoil: ", aerofoil_dat)
  
  # Run xfoil with the generated aerofoil dat file
  run_xfoil(aerofoil_dat)
  
  # Read the results
  cl_max, aoa_max = read_pol(aerofoil_dat.replace('.dat', '.pol'))

  cl_max = float(cl_max)
  
  # We want to maximize cl_max, so we return its negative for minimization
  print("Current cl_max: ", cl_max)
  return -cl_max, aoa_max, aerofoil_dat


# Baseline parameters
# immutable parameters
thickness = 22
chordlength = 235
x_le = 2.468

baseline = [2.143, 0.803]

# Define the bounds for the parameters
bounds = [(2, 2.3), (.78, 0.85)] # x_le, x_te, S

# Perform the optimization
result = minimize(lambda param: objective_function(param)[0], baseline, bounds=bounds, method='L-BFGS-B', options={'ftol': 1e-9, 'eps': 1e-3, 'maxiter':100})
#result = differential_evolution(lambda param: objective_function(param)[0], bounds, maxiter=10, popsize=10, disp=True)

# Extract the optimal parameters
optimal_params = result.x
optimal_aerofoil_dat = objective_function(optimal_params)[2]
optimal_cl_max = objective_function(optimal_params)[0]
optimal_aoa_max = objective_function(optimal_params)[1]

list_of_x_te = []
list_of_S = []
list_of_cl = []

for pol in os.listdir(os.getcwd() + '\pol'):
  list_of_x_te.append(float(pol.split('_')[2]))
  list_of_S.append(float(pol.split('_')[3].split('.')[0]))
  list_of_cl.append(read_pol(pol)[0])

# Plot the results and colour by cl
plt.scatter(list_of_x_te, list_of_S, c=list_of_cl, cmap='viridis')
plt.colorbar()
plt.xlabel('x_te')
plt.ylabel('S')
plt.title('Optimisation of Ozgoose Foil')
plt.show()

'''
ozgoose = [4, 6, 0.75]
ozgoose_dat = create_pollock_aerofoil(thickness, chordlength, ozgoose[0], ozgoose[1], ozgoose[2])
run_xfoil(ozgoose_dat)
ozgoose_cl_max, ozgoose_aoa_max = read_pol(ozgoose_dat.replace('.dat', '.pol'))
print(f'Ozgoose cl_max: {ozgoose_cl_max}')'''