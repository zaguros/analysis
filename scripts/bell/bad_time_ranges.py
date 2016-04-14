import datetime
bad_time_ranges = [[datetime.datetime(2015,3, 4,10, 0,33),      datetime.datetime(2015,3, 4,10,14, 4)], #EWI lights on
                   [datetime.datetime(2015,3, 9,15,50,49),      datetime.datetime(2015,3, 9,16, 4, 3)], #EWI lights on
                   [datetime.datetime(2015,3,10,11,13, 6),      datetime.datetime(2015,3,10,11,19,22)], #EWI lights on
                   [datetime.datetime(2015,3,11, 8,13,31),      datetime.datetime(2015,3,11, 8,18,57)], #EWI lights on
                   [datetime.datetime(2015,3,11,17,32,42),      datetime.datetime(2015,3,11,17,56, 4)], #EWI lights on
                   [datetime.datetime(2015,3,12,12,55,24),      datetime.datetime(2015,3,12,12,57,22)], #EWI lights on
                   [datetime.datetime(2015,3,19,12,31,54),      datetime.datetime(2015,3,19,12,43,26)], #EWI lights on
                   [datetime.datetime(2015,3,25, 5,51, 0),      datetime.datetime(2015,3,25, 7, 0, 0)], #EWI lights on
                   [datetime.datetime(2015,3,26, 7,48,54),      datetime.datetime(2015,3,26, 7,49,13)], #EWI lights on
                   [datetime.datetime(2015,3,26,14,56,58),      datetime.datetime(2015,3,26,15, 8,43)], #EWI lights on
                   [datetime.datetime(2015,4,23,14,31, 8),      datetime.datetime(2015,4,23,14,36,22)], #EWI lights on
                   [datetime.datetime(2015,4,24,11,14,54),      datetime.datetime(2015,4,24,11,23,31)], #EWI lights on
                   [datetime.datetime(2015,4,28, 1,30, 0),      datetime.datetime(2015,4,28,20, 0, 0)], #ch1-ch0 offset 4th
                   [datetime.datetime(2015,4,29, 6,40, 0),      datetime.datetime(2015,4,29, 7,40, 0)], #ch1-ch0 offset 4th
                   [datetime.datetime(2015,4,29,16,13, 0),      datetime.datetime(2015,4,29,23,30, 0)], #ch1-ch0 offset 4th
                   [datetime.datetime(2015,5,15, 8,50,00),      datetime.datetime(2015,5,15,10, 0, 0)], #ch1-ch0 offset ZZ
                   [datetime.datetime(2015,5,19,23,22,13),      datetime.datetime(2015,5,20,10,25,10)], #lasermeister crash ZZ
                   [datetime.datetime(2015,6,12,21,11,00),      datetime.datetime(2015,6,13,04,25,10)], #ch1-ch0 offset XX
                   [datetime.datetime(2015,6,13,15,31,44),      datetime.datetime(2015,6,13,18,22,23)], # Bell_optimizer not running on LT3 XX
                   [datetime.datetime(2015,6,13,18,22,23),      datetime.datetime(2015,6,13,20, 8,23)], # Bell_optimizer not ruuning on LT3 XX
                   [datetime.datetime(2015,6,14, 4,34,07),      datetime.datetime(2015,6,14, 9, 0, 0)], #trouble with yellow laser or WM from LT4 XX
                   [datetime.datetime(2015,6,24,14,39,03),      datetime.datetime(2015,6,24,15, 0,38)], #EWI lights on
                   [datetime.datetime(2015,6,24,19, 0, 0),      datetime.datetime(2015,6,25,02,0, 0)],#datetime.datetime(2015,6,24,20,40, 0)], #ch1=ch0 offset XX
                   [datetime.datetime(2015,6,25,19, 0, 0),      datetime.datetime(2015,6,25,23, 0, 0)], #ch1=ch0 offset XX
                   [datetime.datetime(2015,6,29,10,24,41),      datetime.datetime(2015,6,29,12,44, 0)], #EWI lights on
                   [datetime.datetime(2015,7,01, 0, 0, 0),      datetime.datetime(2015,7,01,01, 0, 0)], # qt lab not restarted properly
                   ]