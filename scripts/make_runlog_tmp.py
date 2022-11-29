of = open('runlog.txt', 'w')
of.write("run_number , run_description , start_time , exposure_sec , GEM3_V , GEM2_V , GEM1_V , T2_V , T1_V , DRIFT_V , OFFSET_V , PMT1_V , PMT2_V , PMT3_V , PMT4_V , HV_STATE , stop_time , number_of_events\n")

pattFe = "%d  S002:DATA:Fe, with pmt Trg: LVL2, 5mV thr, veto 10 us, Fe 25cm far from GEMs; gas flux 20 l/h , 2022-11-24 05:16:46 , 0.3 , 420 , 420 , 420 , 500 , 500 , 960 , 0 , 813 , 836 , 774 , 770 , 0 , 2022-11-24 05:17:42 , 200\n"
pattPed = pattFe.replace("DATA:Fe","PED:Fe")
ped = 6121

of.write(pattPed % ped)
for r in range(6122,6291):
    of.write(pattFe % r)

of.close()

    
