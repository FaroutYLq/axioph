import straxion
st = straxion.qualiphide()
rr = st.make("timeS429", "raw_records", config=dict(daq_input_dir="/Users/lanqingyuan/Documents/GitHub/axioph/QUALIPHIDE/debug_straxion/timeS429/", record_length=5000000, fs=500_000))