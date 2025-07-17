import straxion
st = straxion.qualiphide()
r = st.get_array(
    "timeS429", # run id
    "records", # data type
    config=dict(
        daq_input_dir="/Users/lanqingyuan/Documents/GitHub/axioph/QUALIPHIDE/debug_straxion/timeS429", 
        record_length=5_000_000, 
        fs=500_000, 
        iq_finescan_dir="/Users/lanqingyuan/Documents/GitHub/axioph/QUALIPHIDE/debug_straxion/finescan"
    )
)