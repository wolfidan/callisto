
import datetime

def convert_decimal_hour_timedelta(time_df):
    '''
    This function converts the times of the log file of the motors to timedelta
    '''
    hhmmss = []
    for t in time_df:
        
        # Extract hours
        hh = int(t)
        # Extract minutes
        mm = int((t - hh) * 60)
        # Extract seconds
        ss = int(((t - hh) * 60 - mm) * 60)
        # Extract milliseconds
        mil = int(((t - hh) * 60 - mm - ss/60) * 60000)
        # Format time 
        hhmmss.append(datetime.timedelta(hours = hh, minutes = mm, seconds = ss,
                                         milliseconds = mil))
    
    return np.array(hhmmss)


#**************************************************************************

def average_around_time(times, data, ref_time, avg_period = 1.5):
    dic = {"timestamp": times, "value": data}
    df = pd.DataFrame(dic)
    
    # Ensure timestamp is a datetime type and sorted
    df['timestamp'] = pd.to_datetime(df['timestamp'])
    df = df.sort_values('timestamp')

    # Set timestamp as index for efficient slicing
    df = df.set_index('timestamp')

    # Duration of the averaging window
    delta = pd.Timedelta(seconds=avg_period)

    # Function to compute mean in the time window around a given time
    def average_around_date(date):
        window = df.loc[date - delta : date + delta]
        return window['value'].mean()

    # Apply the function to each date
    averages = [average_around_date(i) for i in ref_time]

    return np.array(averages)

