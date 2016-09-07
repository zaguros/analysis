"""
script that takes an array of times, NV_nrs, and prints a table.
"""

import numpy as np

# time_NVs = {}

# time_NVs['20160719011038']=(1,8,10)
# time_NVs['20160719011938'] = (5,14)
# time_NVs['20160719012810'] = (9)

# time_NVs['20160719050837']=() #scan2d_Sophiex=-10,y=70,z=62.4 + 2.5
# time_NVs['20160719050006']=() #scan2d_Sophiex=-10,y=70,z=62.4 + 1.5
# time_NVs['20160719045134']=() #scan2d_Sophiex=-10,y=70,z=62.4 + 0.5
# time_NVs['20160719044233']=() #scan2d_Sophiebleaching_x=-10,y=70,z=62.4
# time_NVs['20160719043318']=() #scan2d_Sophiex=-30,y=70,z=63.3 + 2.5
# time_NVs['20160719042446']=() #scan2d_Sophiex=-30,y=70,z=63.3 + 1.5
# time_NVs['20160719041615']=() #scan2d_Sophiex=-30,y=70,z=63.3 + 0.5
# time_NVs['20160719040714']=() #scan2d_Sophiebleaching_x=-30,y=70,z=63.3
# time_NVs['20160719035759']=() #scan2d_Sophiex=-50,y=70,z=63.8 + 2.5
# time_NVs['20160719034926']=() #scan2d_Sophiex=-50,y=70,z=63.8 + 1.5
# time_NVs['20160719034055']=() #scan2d_Sophiex=-50,y=70,z=63.8 + 0.5
# time_NVs['20160719033155']=() #scan2d_Sophiebleaching_x=-50,y=70,z=63.8
# time_NVs['20160719032239']=() #scan2d_Sophiex=-70,y=70,z=65.0 + 2.5
# time_NVs['20160719031408']=() #scan2d_Sophiex=-70,y=70,z=65.0 + 1.5
# time_NVs['20160719030537']=() #scan2d_Sophiex=-70,y=70,z=65.0 + 0.5
# time_NVs['20160719025637']=() #scan2d_Sophiebleaching_x=-70,y=70,z=65.0
# time_NVs['20160719024722']=() #scan2d_Sophiex=-90,y=70,z=61.3 + 2.5
# time_NVs['20160719023851']=() #scan2d_Sophiex=-90,y=70,z=61.3 + 1.5
# time_NVs['20160719023019']=() #scan2d_Sophiex=-90,y=70,z=61.3 + 0.5
# time_NVs['20160719022119']=() #scan2d_Sophiebleaching_x=-90,y=70,z=61.3
# time_NVs['20160719021205']=() #scan2d_Sophiex=-10,y=55,z=62.3 + 2.5
# time_NVs['20160719020333']=() #scan2d_Sophiex=-10,y=55,z=62.3 + 1.5
# time_NVs['20160719015502']=() #scan2d_Sophiex=-10,y=55,z=62.3 + 0.5
# time_NVs['20160719014556']=() #scan2d_Sophiebleaching_x=-10,y=55,z=62.3
# time_NVs['20160719013641']=() #scan2d_Sophiex=-30,y=55,z=62.6 + 2.5
# time_NVs['20160719012810']=() #scan2d_Sophiex=-30,y=55,z=62.6 + 1.5
# time_NVs['20160719011938']=() #scan2d_Sophiex=-30,y=55,z=62.6 + 0.5
# time_NVs['20160719011038']=() #scan2d_Sophiebleaching_x=-30,y=55,z=62.6
# time_NVs['20160719010122']=() #scan2d_Sophiex=-50,y=55,z=63.4 + 2.5
# time_NVs['20160719005250']=() #scan2d_Sophiex=-50,y=55,z=63.4 + 1.5
# time_NVs['20160719004419']=() #scan2d_Sophiex=-50,y=55,z=63.4 + 0.5
# time_NVs['20160719003519']=() #scan2d_Sophiebleaching_x=-50,y=55,z=63.4
# time_NVs['20160719002603']=() #scan2d_Sophiex=-70,y=55,z=64.2 + 2.5
# time_NVs['20160719001732']=() #scan2d_Sophiex=-70,y=55,z=64.2 + 1.5
# time_NVs['20160719000900']=() #scan2d_Sophiex=-70,y=55,z=64.2 + 0.5
# time_NVs['20160719000000']=() #scan2d_Sophiebleaching_x=-70,y=55,z=64.2




class NV(m2.M2Analysis):







print time_NVs