# Imports 
import numpy as np

# Any initialization step can be done here
# This method is called during Packagename::initialize_impl
def init ():
    pass

# Probably what we would call during run_impl
def main (p_mid, T_mid, fname=None):
    # Print out max value of p_mid and T_mid
    if fname is not None:
        with open(fname, 'a') as log_file:
            log_file.write("\t Python - Max value for p_mid: ", np.max(p_mid))
            log_file.write("\t Python - Max value for T_mid: ", np.max(T_mid))

            # scale both by 0.5
            p_mid = p_mid * 0.5
            T_mid = T_mid * 0.5

            # log max values
            log_file.write("\t Python - Max value for p_mid scaled by 0.5: ", np.max(p_mid))
            log_file.write("\t Python - Max value for T_mid scaled by 0.5: ", np.max(T_mid))

            #scale both by 2.0
            p_mid = p_mid * 2.0
            T_mid = T_mid * 2.0

            # log max values
            log_file.write("\t Python - Max value for p_mid scaled by 2.0: ", np.max(p_mid))
            log_file.write("\t Python - Max value for T_mid scaled by 2.0: ", np.max(T_mid))

    else:
        # scale both by 0.5
        p_mid = p_mid * 0.5
        T_mid = T_mid * 0.5

        #scale both by 2.0
        p_mid = p_mid * 2.0
        T_mid = T_mid * 2.0

# Define any other needed functions