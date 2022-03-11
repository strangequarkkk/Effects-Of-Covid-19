
################################################################
#  Data and useful code for the Computing Extended Coursework  #
################################################################


# You can copy the data lists from this file and paste them into your Jupyter Notebook
# to create the numpy arrays you will need for the extended coursework question.
# You can also copy the sample code at the bottom to use and modify for your own model of the pandemic.


#########################################################################################################

# The following data list contains the number of weekly deaths (with a positive covid test within 28 days) in the UK, from early March until just before Christmas:

data_list = [  9.,  150.,  722., 2898., 5909., 6338., 5773., 4881., 3638.,
       2671., 2075., 1890., 1207.,  937.,  591.,  480.,  362.,  254.,
        164.,  107.,  113.,   89.,   89.,   56.,   74.,   50.,   81.,
         97.,  197.,  300.,  390.,  701., 1054., 1608., 2165., 2808.,
       2847., 3258., 3085., 2986., 2970.]

#########################################################################################################

# The following data list contains the cumulative sum of positive covid tests, week by week from early March until just before Christmas:

cases_list = [ 797.0, 3983.0, 14548.0, 39282.0, 73758.0, 108692.0, 143464.0,
       177454.0, 211364.0, 236711.0, 254195.0, 271222.0, 283311.0, 292950.0,
       301815.0, 309360.0, 284276.0, 288133.0, 293239.0, 297914.0, 303181.0,
       309005.0, 316367.0, 323313.0, 331644.0, 342351.0, 361677.0, 385936.0,
       423236.0, 467146.0, 575679.0, 689257.0, 830998.0, 989745.0, 1146484.0,
       1317496.0, 1473508.0, 1589301.0, 1690432.0, 1809455.0, 1977167.0]

########################################################################################################

# The following code is an implementation of the simplest SIR model for the dynamics of epidemics,
# under the assumption of a constant population (i.e. neglecting birth and death) 


# Maximum number of steps (weeks) to evolve in time
tmax=80
# Define empty arrays for the three variables of the SIR model
Iarray=np.zeros(tmax)
Sarray=np.zeros(tmax)
Rarray=np.zeros(tmax)
# Model parameters
R0=1.6  # Basic Reproduction number
gamma=1/2. # This is the inverse of the infectious period (in weeks)

N=70000000 # Population taken to be 70 million
Iarray[0]=100 # Start with 100 infectious individuals on week 1 
Rarray[0]=0
Sarray[0]=N-Iarray[0]-Rarray[0]
# Initialise I, S and R to be used in the loop 
I=Iarray[0]
S=Sarray[0]
R=Rarray[0]
# For-loop to update I, S, R and store in arrays Iarray, Sarray, Rarray 
for t in np.arange(1,tmax):
    DeltaI = I * ( R0*gamma*S/(S+I+R) - gamma ) # Read from eqn (3), Deltat=1 week 
    DeltaS = - I * R0*gamma*S/(S+I+R) # Read from eqn (4), Deltat=1 week 
    I += DeltaI # Update I 
    S += DeltaS # Update S
    R = N - I - S # No need to solve eqn (5). To find R use eqn (2) instead.
    Iarray[t] = I # Store updated values in corresponding array variable
    Sarray[t] = S
    Rarray[t] = R 
