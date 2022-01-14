''' ################### SET/CHECK THESE PARAMETERS BEFORE RUNNING ################## '''
LOWPASS_FREQ = 500  # Hz
STIM_TIME = 520     # ms
POST_STIM = 250     # ms, amount of time after stimulus to look for max value
TP_START = 5        # ms, time of start of test pulse
VM_JUMP = 10        # mV, test pulse voltage jump
PRE_TP = 3          # ms, amount of time before test pulse start to get baseline
UNIT_SCALER = -12   # unitless, scaler to get back to A, from pA
AMP_FACTOR = 1      # scaler for making plots in pA
FS = 25             # kHz, the sampling frequency