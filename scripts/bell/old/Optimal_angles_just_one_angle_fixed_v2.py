########## DENSITY MATRIX CALCULATIONS ##########

## IMPORT MODULES ##
import numpy as np
import scipy.optimize as opt
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.colors import BoundaryNorm
import time
from analysis.scripts.bell import Functions_optimal_angle_simulation as Func_sim
######################################################################
"""
If save is set to true all plots are saved with datestamps and names
corresponding to settings in the folder in which the simulation runs.
"""
save = False
plot_others=False

## Define Qubit state and calculate density matrix ##
# Give your wave function in basis |00>, |01>, |10>, |11> #
# Make sure to normalize #
Psi = np.array([[0.0],[1.0],[1.0],[0.0]]) * 1/np.sqrt(2)
Rho_pure = np.dot(Psi, Psi.transpose())

## Define matrices representing most important errors (Visisbility and dark counts (Error_dist and Error_dark)) ##
Error_dist = 1./2. * np.array([[0.0, 0.0, 0.0, 0.0],[0.0, 1.0, 0.0, 0.0],
                       [0.0, 0.0, 1.0, 0.0],[0.0, 0.0, 0.0, 0.0]])

Error_dark = 1./4. * np.array([[1.0, 0.0, 0.0, 0.0],[0.0, 1.0, 0.0, 0.0],
                       [0.0, 0.0, 1.0, 0.0],[0.0, 0.0, 0.0, 1.0]])

## Define Readout Fidelities ##
"""
Give the readout fidelity of LT1 up state (ms = 0) in F[0] and 
the readout fidelity of the LT1 down state (ms = 1) in F[1], 
do the same for LT3 in for G[0] and G[1].
"""
F = np.array([1,1])
G = np.array([1.,1])

"""
Set 1- visibility as D and the darkcount percentage as dc. 
"""
D = 20./100
dc = 3.3/100

## Construct finale state used for angle optimization ##
Rho = Rho_pure
Rho_error = Func_sim.Include_errors(Rho_pure,Error_dist,Error_dark,D,dc)

## 2D plot visibility, dark counts versus CHSH ##
def vis_dark_CHSH(Rho_pure,Error_dist,Error_dark,F,G, num_points_2D, vis_min):
    """
    Makes a 2D plot of with the visibility on the x-axis and the dark
    count on the y-axis. The values in plot are the CHSH values. 
    The values are calculated for the nnumber of points entered in the
    function. The readout fidelity is entered in the function through F
    ang G.
    """
    N = num_points_2D
    Max_D = 1.0
    Max_dc = 0.03

    Ang_a = np.zeros((N,N))
    Ang_b = np.zeros((N,N))
    Ang_a_prime = np.zeros((N,N))
    Ang_b_prime = np.zeros((N,N))
    CHSH_val = np.zeros((N,N))
    CHSH_val_non_op = np.zeros((N,N))
    Dif_a_b = np.zeros((N,N))
    Dif_a_b_prime = np.zeros((N,N))
    Dif_a_prime_b = np.zeros((N,N))
    Dif_a_prime_b_prime = np.zeros((N,N))
    Dif_b_b_prime = np.zeros((N,N))
    Dif_a_a_prime = np.zeros((N,N))
    
    # set range for D
    x = np.linspace(vis_min,Max_D,N)
    # set range for dc
    y = np.linspace(0.,Max_dc,N)

    X,Y = np.meshgrid(x,y)


                
    print "Calculating data for fidelity graph"
    for i in range(N):    #x-axis D
            for j in range(N): #y-axis dc
            # Define readout fidelities, all states equal #
                    D = Max_D/(N-1) *i
                    dc = Max_dc/(N-1) *j
                    Rho = Func_sim.Include_errors(Rho_pure,Error_dist,Error_dark,D,dc)
                    Minimization = opt.minimize(Func_sim.Invert, np.array([-22.5,45.,67.5]),
                       args = tuple([Rho,F,G]), method = 'SLSQP',
                       bounds = ((-90.,90.),(-90.,90.),(-90.,90.)),
                       tol=1e-12, options={'disp': False})
                    Ang = Minimization.x
                    Ang_b[j,i] = Ang[0]
                    Ang_a_prime[j,i] = Ang[1]
                    Ang_b_prime[j,i] = Ang[2]
                    Dif_a_b[j,i] = abs(0. - Ang[0])
                    Dif_a_b_prime[j,i] = abs(0. - Ang[2])
                    Dif_a_prime_b[j,i] = abs(Ang[1]- Ang[0])
                    Dif_a_prime_b_prime[j,i] = abs(Ang[1] - Ang[2])
                    Dif_b_b_prime[j,i] = abs(Ang[0] - Ang[2])
                    Dif_a_a_prime[j,i] = abs(0. - Ang[1]) 



                    CHSH_val[j,i] = Func_sim.CHSH(Ang,Rho,F,G)
                    Ang_non_op = [Ang_b[0,0],Ang_a_prime[0,0],Ang_b_prime[0,0]]
                    CHSH_val_non_op[j,i] = Func_sim.CHSH(Ang_non_op,Rho,F,G)
    print "All necessary data calculated"

    Time = time.strftime("%Y-%m-%d_t%H-%M-%S")
    F0 = str(F[0])
    F0 = F0.replace(".","dot")
    F1 = str(F[1])
    F1 = F1.replace(".","dot")
  
    plt.figure(1)
    levels = np.arange(1.2,3.,0.1)

    plot2 = plt.contourf(X, Y, CHSH_val, levels)
    cbar = plt.colorbar(plot2)
    plt.contour(plot2,levels = [2.], colors = 'r',hold='on' )
    plt.axis([X.min(), X.max(), Y.min(), Y.max()])
    plt.title('With optimized angles')
    
    cbar.ax.set_ylabel('CHSH value')
    plt.xlabel('1-Visibility')
    plt.ylabel('Probability fake entanglement')
    
    if save:
        filename = Time + '_CHSH_val_vs_D_dc' + '_F0G0_eq' + F0 + '_F1G1_eq' + F1
        plt.savefig(filename + ".png")
        plt.savefig(filename + ".eps")

    plt.figure(15)
    levels = np.arange(1.2,3.,0.1)
    plot1 = plt.contourf(X, Y, CHSH_val_non_op,levels)
    cbar = plt.colorbar(plot1)
    plt.contour(plot1,levels=[2], colors = 'r',hold='on' )
    plt.axis([X.min(), X.max(), Y.min(), Y.max()])
    plt.title('With standard angles')

    cbar.ax.set_ylabel('CHSH value')
    plt.xlabel('1-Visibility')
    plt.ylabel('Probability fake entanglement')

    if save:
        filename = Time + '_CHSH_val_non_op_vs_D_dc' + '_F0G0_eq' + F0 + '_F1G1_eq' + F1
        plt.savefig(filename + ".png")
        plt.savefig(filename + ".eps")

    if plot_others:
        print 'dfssdfsd'
        plt.figure(2)
        plot2 = plt.pcolor(X, Y, Ang_a, vmin = -2, vmax = 2)
        plt.axis([X.min(), X.max(), Y.min(), Y.max()])
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('Angle a')
        plt.xlabel('1-Visibility')
        plt.ylabel('Probability fake entanglement')
        
        if save:
            filename = Time + '_Angle_a_vs_D_dc' + '_F0G0_eq' + F0 + '_F1G1_eq' + F1
            plt.savefig(filename + ".png")
            plt.savefig(filename + ".eps")
       

        plt.figure(3)
        plot3 = plt.pcolor(X, Y, Ang_b)
        plt.axis([X.min(), X.max(), Y.min(), Y.max()])
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('Angle b')
        plt.xlabel('1-Visibility')
        plt.ylabel('Probability fake entanglement')
        
        if save:
            filename = Time + '_Angle_b_vs_D_dc' + '_F0G0_eq' + F0 + '_F1G1_eq' + F1
            plt.savefig(filename + ".png")
            plt.savefig(filename + ".eps")
       

        plt.figure(4)
        plot4 = plt.pcolor(X, Y, Ang_a_prime, vmin = -90, vmax = 90)
        plt.axis([X.min(), X.max(), Y.min(), Y.max()])
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('Angle a_prime')
        plt.xlabel('1-Visibility')
        plt.ylabel('Probability fake entanglement')
        
        if save:
            filename = Time + '_Angle_a_prime_vs_D_dc' + '_F0G0_eq' + F0 + '_F1G1_eq' + F1
            plt.savefig(filename + ".png")
            plt.savefig(filename + ".eps")
       

        plt.figure(5)
        plot5 = plt.pcolor(X, Y, Ang_b_prime)
        plt.axis([X.min(), X.max(), Y.min(), Y.max()])
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('Angle b_prime')
        plt.xlabel('1-Visibility')
        plt.ylabel('Probability fake entanglement')
        
        if save:
            filename = Time + '_Angle_b_prime_vs_D_dc' + '_F0G0_eq' + F0 + '_F1G1_eq' + F1
            plt.savefig(filename + ".png")
            plt.savefig(filename + ".eps")

        plt.figure(6)
        plot6 = plt.pcolor(X, Y, Dif_a_b)
        plt.axis([X.min(), X.max(), Y.min(), Y.max()])
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('Abs(a-b)')
        plt.xlabel('1-Visibility')
        plt.ylabel('Probability fake entanglement')
        
        if save:
            filename = Time + '_Dif_a_b_vs_D_dc' + '_F0G0_eq' + F0 + '_F1G1_eq' + F1
            plt.savefig(filename + ".png")
            plt.savefig(filename + ".eps")

        plt.figure(7)
        plot7 = plt.pcolor(X, Y, Dif_a_b_prime)
        plt.axis([X.min(), X.max(), Y.min(), Y.max()])
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('Abs(a-b_prime)')
        plt.xlabel('1-Visibility')
        plt.ylabel('Probability fake entanglement')
        
        if save:
            filename = Time + '_Dif_a_b_prime_vs_D_dc' + '_F0G0_eq' + F0 + '_F1G1_eq' + F1
            plt.savefig(filename + ".png")
            plt.savefig(filename + ".eps")

        plt.figure(8)
        plot8 = plt.pcolor(X, Y, Dif_a_prime_b)
        plt.axis([X.min(), X.max(), Y.min(), Y.max()])
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('Abs(a_prime-b)')
        plt.xlabel('1-Visibility')
        plt.ylabel('Probability fake entanglement')
        
        if save:
            filename = Time + '_Dif_a_prime_b_vs_D_dc' + '_F0G0_eq' + F0 + '_F1G1_eq' + F1
            plt.savefig(filename + ".png")
            plt.savefig(filename + ".eps")
       
        plt.figure(9)
        plot9 = plt.pcolor(X, Y, Dif_a_prime_b_prime)
        plt.axis([X.min(), X.max(), Y.min(), Y.max()])
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('Abs(a_prime-b_prime)')
        plt.xlabel('1-Visibility')
        plt.ylabel('Probability fake entanglement')
        
        if save:
            filename = Time + '_Dif_a_prime_b_prime_vs_D_dc' + '_F0G0_eq' + F0 + '_F1G1_eq' + F1
            plt.savefig(filename + ".png")
            plt.savefig(filename + ".eps")

        plt.figure(17)
        plot10 = plt.pcolor(X, Y, Dif_b_b_prime)
        plt.axis([X.min(), X.max(), Y.min(), Y.max()])
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('Abs(b-b_prime)')
        plt.xlabel('1-Visibility')
        plt.ylabel('Probability fake entanglement')
        
        if save:
            filename = Time + '_Dif_b_b_prime_vs_D_dc' + '_F0G0_eq' + F0 + '_F1G1_eq' + F1
            plt.savefig(filename + ".png")
            plt.savefig(filename + ".eps")

        plt.figure(18)
        plot11 = plt.pcolor(X, Y, Dif_a_a_prime)
        plt.axis([X.min(), X.max(), Y.min(), Y.max()])
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('Abs(a-a_prime)')
        plt.xlabel('1-Visibility')
        plt.ylabel('Probability fake entanglement')
        
        if save:
            filename = Time + '_Dif_a_a_prime_vs_D_dc' + '_F0G0_eq' + F0 + '_F1G1_eq' + F1
            plt.savefig(filename + ".png")
            plt.savefig(filename + ".eps")

## Plots CHSH value vs. ROF for both pure and non pure states ##
def Fidelity(Rho1,Rho2, num_points_1D, min_ROF, max_ROF):
        """
        Makes a plot of the CHSH value vs. the Readout fidelity with the 
        limits as set above. Now it sweeps the up state (ms =  0) but this
        can be changed easily.
        """
        Ang_a = []
        Ang_b = []
        Ang_a_prime = []
        Ang_b_prime = []
        CHSH_val = []
        min_ROF = 0.8
        max_ROF = 1
        N = num_points_1D
        print dc
        
        print "Calculating data for fidelity graph"
        for i in np.linspace(min_ROF,max_ROF,N):
                # Define readout fidelities, all states equal #
                F = np.array([i, 1])
                G = np.array([i, 1])
                Minimization = opt.minimize(Func_sim.Invert, np.array([-22.5,45.,67.5]),
                           args = tuple([Rho1,F,G]), method = 'SLSQP',
                           bounds = ((-90.,90.),(-90.,90.),(-90.,90.)),
                           tol=1e-12, options={'disp': False})
                Ang = Minimization.x
                Ang_a.append(0.)
                Ang_b.append(Ang[0])
                Ang_a_prime.append(Ang[1])
                Ang_b_prime.append(Ang[2])
                CHSH_val.append(Func_sim.CHSH(Ang,Rho1,F,G))
        print "All necessary data calculated"
        print "The maximal CHSH value is"
        print CHSH_val[N-1]
        print "The maximal violation of the CHSH value is"
        print 2 * np.sqrt(2)

        ## Plotting ##

        Time = time.strftime("%d-%m-%Y_%H-%M-%S")
        D_str = str(D)
        D_str = D_str.replace(".","dot")
        dc_str = str(dc)
        dc_str = dc_str.replace(".","dot")
        F1 = str(F[1])
        F1 = F1.replace(".","dot")
        
        plt.figure(10)
        host = host_subplot(111, axes_class=AA.Axes)
        plt.subplots_adjust(right=0.90)

        par1 = host.twinx()
        host.set_xlim(min_ROF, max_ROF)
        host.set_ylim(-91.0, 91.)
        host.set_xlabel("Read out Fidelity LT1 = ms0,LT3 = ms0")
        host.set_ylabel("Angle($^\circ$)")
        par1.set_ylabel("CHSH value")

        p1, = host.plot(np.linspace(min_ROF, max_ROF,N),Ang_a, label="a")
        p2, = host.plot(np.linspace(min_ROF, max_ROF,N),Ang_b, label="b")
        p3, = host.plot(np.linspace(min_ROF, max_ROF,N),Ang_a_prime, label="a'")
        p4, = host.plot(np.linspace(min_ROF, max_ROF,N),Ang_b_prime, label="b'")
        p5, = par1.plot(np.linspace(min_ROF, max_ROF,N),CHSH_val, label="CHSH value pure state")
        p6, = par1.plot(np.linspace(min_ROF, max_ROF,N), 2*np.ones(N), '--g', label = "Maximal classical CHSH value")
        p7, = par1.plot(np.linspace(min_ROF, max_ROF,N), 2*np.sqrt(2)*np.ones(N), '--k',label = "Maximal CHSH violation")

        par1.set_ylim(1, 3)

        legend = host.legend(loc='upper left', shadow=False)
        #host.legend()
        host.axis["left"].label.set_color(p1.get_color())
        par1.axis["right"].label.set_color(p5.get_color())

        plt.draw()
        
        if save:
            filename = Time + '_Ang_pure_vs_F0G0' + '_D_eq' + D_str + '_dc_eq' + dc_str + '_F1G1_eq'+ F1
            plt.savefig(filename + ".png")
            plt.savefig(filename + ".eps")

        ## End Plotting #

        # In order to compare two different density matrices uncomment the bottum part and add Rho2 in the function
        
        Ang_a2 = []
        Ang_b2 = []
        Ang_a2_prime = []
        Ang_b2_prime = []
        CHSH_val_2 = []

        print "Calculating data for fidelity graph"
        for i in np.linspace(min_ROF, max_ROF,N):
                # Define readout fidelities, all states equal #
                F = np.array([i, 1])
                G = np.array([i, 1])
                Minimization = opt.minimize(Func_sim.Invert, np.array([-22.5,45.,67.5]),
                           args = tuple([Rho2,F,G]), method = 'SLSQP',
                           bounds = ((-90.,90.),(-90.,90.),(-90.,90.)),
                           tol=1e-12, options={'disp': False})
                Ang = Minimization.x
                Ang_a2.append(0.)
                Ang_b2.append(Ang[0])
                Ang_a2_prime.append(Ang[1])
                Ang_b2_prime.append(Ang[2])
                CHSH_val_2.append(Func_sim.CHSH(Ang,Rho2,F,G))
        print "All necessary data calculated"

        ## Plotting ##

        plt.figure(11)
        host = host_subplot(111, axes_class=AA.Axes)
        plt.subplots_adjust(right=0.90)

        par1 = host.twinx()
        host.set_xlim(min_ROF, max_ROF)
        host.set_ylim(-91.0, 91.)
        host.set_xlabel("Read out Fidelity LT1 = ms0, LT3 = ms0")
        host.set_ylabel("Angle($^\circ$)")
        par1.set_ylabel("CHSH value")

        p1, = host.plot(np.linspace(min_ROF, max_ROF,N),Ang_a2, label="a")
        p2, = host.plot(np.linspace(min_ROF, max_ROF,N),Ang_b2, label="b")
        p3, = host.plot(np.linspace(min_ROF, max_ROF,N),Ang_a2_prime, label="a'")
        p4, = host.plot(np.linspace(min_ROF, max_ROF,N),Ang_b2_prime, label="b'")
        p5, = par1.plot(np.linspace(min_ROF, max_ROF,N),CHSH_val_2, label="CHSH value state including errors")
        p6, = par1.plot(np.linspace(min_ROF, max_ROF,N), 2*np.ones(N), '--g', label = "Maximal classical CHSH value")
        p7, = par1.plot(np.linspace(min_ROF, max_ROF,N), 2*np.sqrt(2)*np.ones(N), '--k',label = "Maximal CHSH violation")

        par1.set_ylim(1, 3)

        legend = host.legend(loc='upper left', shadow=False)
        #host.legend()
        host.axis["left"].label.set_color(p1.get_color())
        par1.axis["right"].label.set_color(p5.get_color())

        plt.draw()
        
        if save:
            filename = Time + '_Ang_inclerrors_vs_F0G0' + '_D_eq' + D_str + '_dc_eq' + dc_str + '_F1G1_eq' + F1
            plt.savefig(filename + ".png")
            plt.savefig(filename + ".eps")

        plt.figure(12)
        host = host_subplot(111,axes_class = AA.Axes)
        plt.subplots_adjust(right = 0.90)

        host.set_xlim(min_ROF, max_ROF)
        host.set_ylim(0.95,3.)
        host.set_xlabel("Read out Fidelity F0,G0")
        host.set_ylabel("CHSH value")

        host.plot(np.linspace(min_ROF, max_ROF,N), CHSH_val, label="CHSH value pure state")
        host.plot(np.linspace(min_ROF, max_ROF,N), CHSH_val_2, label = "CHSH value including errors")
        host.plot(np.linspace(min_ROF, max_ROF,N), 2*np.ones(N), '--r', label = "Maximal classical CHSH value")
        host.plot(np.linspace(min_ROF, max_ROF,N), 2*np.sqrt(2)*np.ones(N), '--c',label = "Maximal CHSH violation")

        legend = host.legend(loc='upper left', shadow=False)

        plt.draw()
        
        if save:
            filename = Time + '_CHSH_vs_F0G0' + '_D_eq' + D_str + '_dc_eq' + dc_str + '_F1G1_eq' + F1
            plt.savefig(filename + ".png")
            plt.savefig(filename + ".eps")

        ## End Plotting ##

## 2D plot visibility, readout fidelity 
def vis_fid_CHSH(Rho_pure,Error_dist,Error_dark, num_points_2D, min_ROF_LT3_ms0, Max_D):
        """
        Makes a 2D plot with 1 - the Visibility on the x-axis and the readout fidelity
        of the ms = 0  state of LT3 on the y-axis.  The CHSH value is plotted in colors.
        The number of points used for the grid is given by num_points_2D
        """
        N = num_points_2D
        dc = 0.01
        Min_G_ms0 = min_ROF_LT3_ms0

        Ang_a = np.zeros((N,N))
        Ang_b = np.zeros((N,N))
        Ang_a_prime = np.zeros((N,N))
        Ang_b_prime = np.zeros((N,N))
        CHSH_val = np.zeros((N,N))
        CHSH_val_non_op = np.zeros((N,N))
        Dif_a_b = np.zeros((N,N))
        Dif_a_b_prime = np.zeros((N,N))
        Dif_a_prime_b = np.zeros((N,N))
        Dif_a_prime_b_prime = np.zeros((N,N))
        Dif_b_b_prime = np.zeros((N,N))
        Dif_a_a_prime = np.zeros((N,N))
        
        # set range for D
        x = np.linspace(0.,Max_D,N)
        # set range for dc
        y = np.linspace(Min_G_ms0,1.,N)

        X,Y = np.meshgrid(x,y)


                
        print "Calculating data for fidelity graph"
        for i in range(N):    #x-axis D (Visibility)
                for j in range(N): #y-axis Read oud Fidelity of G up state (ms = 0)
                # Define Func_sim.readout fidelities, all states equal #
                        D = Max_D/(N-1) *i
                        F = np.array([0.927, 0.99])
                        G = np.array([Min_G_ms0 + (j*(1-Min_G_ms0)/(N-1)),  0.99])
                        Rho = Func_sim.Include_errors(Rho_pure,Error_dist,Error_dark,D,dc)
                        test = np.array([-22.5,45.,67.5])
                        Minimization = opt.minimize(Func_sim.Invert, np.array([-22.5,45.,67.5]),
                           args = tuple([Rho,F,G]), method = 'SLSQP',
                           bounds = ((-90.,90.),(-90.,90.),(-90.,90.)),
                           tol=1e-12, options={'disp': False})
                        Ang = Minimization.x
                        Ang_b[j,i] = Ang[0]
                        Ang_a_prime[j,i] = Ang[1]
                        Ang_b_prime[j,i] = Ang[2]
                        Dif_a_b[j,i] = abs(0. - Ang[0])
                        Dif_a_b_prime[j,i] = abs(0. - Ang[2])
                        Dif_a_prime_b[j,i] = abs(Ang[1]- Ang[0])
                        Dif_a_prime_b_prime[j,i] = abs(Ang[1] - Ang[2])
                        Dif_b_b_prime[j,i] = abs(Ang[0] - Ang[2])
                        Dif_a_a_prime[j,i] = abs(0. - Ang[1]) 



                        CHSH_val[j,i] = Func_sim.CHSH(Ang,Rho,F,G)
                        Ang_non_op = [Ang_b[0,0],Ang_a_prime[0,0],Ang_b_prime[0,0]]
                        CHSH_val_non_op[j,i] = Func_sim.CHSH(Ang_non_op,Rho,F,G)
        print "All necessary data calculated"


        Time = time.strftime("%Y-%m-%d_t%H-%M-%S")
        F0 = str(F[0])
        F0 = F0.replace(".","dot")
        F1 = str(F[1])
        F1 = F1.replace(".","dot")
      
        plt.figure(1)
        levels = np.arange(1.2,3.,0.1)

        plot2 = plt.contourf(X, Y, CHSH_val, levels)
        cbar = plt.colorbar(plot2)
        plt.contour(plot2,levels = [2.], colors = 'r',hold='on' )
        plt.axis([X.min(), X.max(), Y.min(), Y.max()])
        plt.title('CHSH Optimized angles')
        cbar.ax.set_ylabel('CHSH value')
        plt.xlabel('1-Visibility')
        plt.ylabel('Read out fidelity m_s = 0 LT3')
        
        if save:        
            filename = Time + '_CHSH_val_vs_D_G0' + '_F0_eq' + F0 + '_F1G1_eq' + F1
            plt.savefig(filename + ".png")
            plt.savefig(filename + ".eps")

        plt.figure(15)
        levels = np.arange(1.,3.,0.1)
        plot1 = plt.contourf(X, Y, CHSH_val_non_op,levels)
        cbar = plt.colorbar(plot1)
        plt.contour(plot1,levels=[2], colors = 'r',hold='on' )
        plt.axis([X.min(), X.max(), Y.min(), Y.max()])
        plt.title('CHSH Standard angles')

        cbar.ax.set_ylabel('CHSH value')
        plt.xlabel('1-Visibility')
        plt.ylabel('Read out fidelity m_s = 0 LT3')
        
        if save:
            filename = Time + '_CHSH_val_non_op_vs_D_G0' + '_F0_eq' + F0 + '_F1G1_eq' + F1
            plt.savefig(filename + ".png")
            plt.savefig(filename + ".eps")
        if plot_others:
            plt.figure(2)
            plot2 = plt.pcolor(X, Y, Ang_a, vmin = -2, vmax = 2)
            plt.axis([X.min(), X.max(), Y.min(), Y.max()])
            cbar = plt.colorbar()
            cbar.ax.set_ylabel('Angle a')
            plt.xlabel('1-Visibility')
            plt.ylabel('Read out fidelity m_s = 0 LT3')
            
            if save:
                filename = Time + '_Angle_a_vs_D_dc' + '_F0G0_eq' + F0 + '_F1G1_eq' + F1
                plt.savefig(filename + ".png")
                plt.savefig(filename + ".eps")
            

            plt.figure(3)
            plot3 = plt.pcolor(X, Y, Ang_b)
            plt.axis([X.min(), X.max(), Y.min(), Y.max()])
            cbar = plt.colorbar()
            cbar.ax.set_ylabel('Angle b')
            plt.xlabel('1-Visibility')
            plt.ylabel('Read out fidelity m_s = 0 LT3')
            
            if save:
                filename = Time + '_Angle_b_vs_D_dc' + '_F0G0_eq' + F0 + '_F1G1_eq' + F1
                plt.savefig(filename + ".png")
                plt.savefig(filename + ".eps")
            

            plt.figure(4)
            plot4 = plt.pcolor(X, Y, Ang_a_prime, vmin = -90, vmax = 90)
            plt.axis([X.min(), X.max(), Y.min(), Y.max()])
            cbar = plt.colorbar()
            cbar.ax.set_ylabel('Angle a_prime')
            plt.xlabel('1-Visibility')
            plt.ylabel('Read out fidelity m_s = 0 LT3')
            
            if save:
                filename = Time + '_Angle_a_prime_vs_D_dc' + '_F0G0_eq' + F0 + '_F1G1_eq' + F1
                plt.savefig(filename + ".png")
                plt.savefig(filename + ".eps")
            

            plt.figure(5)
            plot5 = plt.pcolor(X, Y, Ang_b_prime)
            plt.axis([X.min(), X.max(), Y.min(), Y.max()])
            cbar = plt.colorbar()
            cbar.ax.set_ylabel('Angle b_prime')
            plt.xlabel('1-Visibility')
            plt.ylabel('Read out fidelity m_s = 0 LT3')       

            if save:
                filename = Time + '_Angle_b_prime_vs_D_dc' + '_F0G0_eq' + F0 + '_F1G1_eq' + F1
                plt.savefig(filename + ".png")
                plt.savefig(filename + ".eps")

            plt.figure(6)
            plot6 = plt.pcolor(X, Y, Dif_a_b)
            plt.axis([X.min(), X.max(), Y.min(), Y.max()])
            cbar = plt.colorbar()
            cbar.ax.set_ylabel('Abs(a-b)')
            plt.xlabel('1-Visibility')
            plt.ylabel('Read out fidelity m_s = 0 LT3')
            
            if save:
                filename = Time + '_Dif_a_b_vs_D_dc' + '_F0G0_eq' + F0 + '_F1G1_eq' + F1
                plt.savefig(filename + ".png")
                plt.savefig(filename + ".eps")

            plt.figure(7)
            plot7 = plt.pcolor(X, Y, Dif_a_b_prime)
            plt.axis([X.min(), X.max(), Y.min(), Y.max()])
            cbar = plt.colorbar()
            cbar.ax.set_ylabel('Abs(a-b_prime)')
            plt.xlabel('1-Visibility')
            plt.ylabel('Read out fidelity m_s = 0 LT3')
            
            if save:
                filename = Time + '_Dif_a_b_prime_vs_D_dc' + '_F0G0_eq' + F0 + '_F1G1_eq' + F1
                plt.savefig(filename + ".png")
                plt.savefig(filename + ".eps")

            plt.figure(8)
            plot8 = plt.pcolor(X, Y, Dif_a_prime_b)
            plt.axis([X.min(), X.max(), Y.min(), Y.max()])
            cbar = plt.colorbar()
            cbar.ax.set_ylabel('Abs(a_prime-b)')
            plt.xlabel('1-Visibility')
            plt.ylabel('Read out fidelity m_s = 0 LT3')
            
            if save:
                filename = Time + '_Dif_a_prime_b_vs_D_dc' + '_F0G0_eq' + F0 + '_F1G1_eq' + F1
                plt.savefig(filename + ".png")
                plt.savefig(filename + ".eps")
            
            plt.figure(9)
            plot9 = plt.pcolor(X, Y, Dif_a_prime_b_prime)
            plt.axis([X.min(), X.max(), Y.min(), Y.max()])
            cbar = plt.colorbar()
            cbar.ax.set_ylabel('Abs(a_prime-b_prime)')
            plt.xlabel('1-Visibility')
            plt.ylabel('Read out fidelity m_s = 0 LT3')
            
            if save:
                filename = Time + '_Dif_a_prime_b_prime_vs_D_dc' + '_F0G0_eq' + F0 + '_F1G1_eq' + F1
                plt.savefig(filename + ".png")
                plt.savefig(filename + ".eps")

            plt.figure(17)
            plot10 = plt.pcolor(X, Y, Dif_b_b_prime)
            plt.axis([X.min(), X.max(), Y.min(), Y.max()])
            cbar = plt.colorbar()
            cbar.ax.set_ylabel('Abs(b-b_prime)')
            plt.xlabel('1-Visibility')
            plt.ylabel('Read out fidelity m_s = 0 LT3')
            
            if save:
                filename = Time + '_Dif_b_b_prime_vs_D_dc' + '_F0G0_eq' + F0 + '_F1G1_eq' + F1
                plt.savefig(filename + ".png")
                plt.ylabel('Read out fidelity m_s = 0 LT3')

            plt.figure(18)
            plot11 = plt.pcolor(X, Y, Dif_a_a_prime)
            plt.axis([X.min(), X.max(), Y.min(), Y.max()])
            cbar = plt.colorbar()
            cbar.ax.set_ylabel('Abs(a-a_prime)')
            plt.xlabel('1-Visibility')
            plt.ylabel('Read out fidelity m_s = 0 LT3')
            
            if save:
                filename = Time + '_Dif_a_a_prime_vs_D_dc' + '_F0G0_eq' + F0 + '_F1G1_eq' + F1
                plt.savefig(filename + ".png")
                plt.savefig(filename + ".eps")

## Call Final Functions ##
""" 
In order to get nice plots only call one function at a time
"""

"""
Determines how many point are taken
"""
num_points_2D = 10
num_points_1D = 10

min_ROF = 0.8
max_ROF = 1.0

#Fidelity(Rho,Rho_error, num_points_1D, min_ROF, max_ROF)

vis_min = 0.8

#vis_dark_CHSH(Rho_pure,Error_dist,Error_dark,F,G, num_points_2D, vis_min)

min_ROF_LT3_ms0 = 0.8

"""
The Maximum of 1- Visibility is given by Max_D
"""
Max_D = 0.4

vis_fid_CHSH(Rho_pure,Error_dist,Error_dark, num_points_2D, min_ROF_LT3_ms0, Max_D)

plt.show()
