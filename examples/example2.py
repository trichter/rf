import obspy.iris
import rf

# Create getwaveform function which has to be passed to rf function
client = obspy.iris.Client()
def getwaveform(station, t1, t2):
    return client.getWaveform('TA', station, '', 'BH?', t1, t2)

# Set output path here or in config file
rf.set_paths('~/obspyDMT-data/client_test')

# Create event file from given events
rf.create_rfeventsfile('./events.xml')

#Calculate receiver functions. Station coordinates are given in stations.txt
rf.rf('client', getwaveform, './stations.txt', deconvolve='freq')
