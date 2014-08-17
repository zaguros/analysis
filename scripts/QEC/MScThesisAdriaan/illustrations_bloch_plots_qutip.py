import qutip as qutip
import numpy as np
import os

### Defining states and vectors to plot.
def figure(save=True) :

    up = qutip.basis(2,0)
    n0 = [1,0,0]
    n1 = [-1,0,0]


    ### Plotting

    b=qutip.Bloch()
    # b.vector_style = 'simple'
    b.add_states(up)
    # b.vector_style = '-|>'
    b.add_vectors([n0,n1])

    b.vector_color = ['r','g','g']



    savename = 'bloch_entangling_gate.svg'
    if save == True:
        savedir = r'/Users/'+os.getlogin()+r'/Documents/MasterThesis/Img/'
        b.save(name = savedir +savename)
        print savedir
    else:
        b.show()
    return b



    ## Saving figure
