import numpy as np 
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import sys
import math
from scipy import linalg
from tqdm import trange


class system(object):

    def __init__(self,*args):
        self.E_c = args[0]
        self.E_j = args[1]
        self.alpha = args[1]/args[0]
        self.ncs = args[2]
        self.Omega_t = np.sqrt(8*self.E_c*self.E_j)
        self.eta = 1/np.sqrt(2)*(8*self.E_c/self.E_j)**0.25
        
class state(object):

    def __init__(self,*args):
        sys = args[0]
        self.p = args[1]
        self.z = args[2]
        self.ncs = sys.ncs
        self.E_c = sys.E_c
        self.alpha = sys.alpha
        self.E_j = sys.E_j
        self.Omega_t = sys.Omega_t
        self.eta = sys.eta

    def exp_Et(self,p,z):
        f = 0
        for i in range(self.ncs):
            for j in range(self.ncs):
                f += (self.Omega_t*np.conj(p[i]*z[i])*p[j]*z[j]-0.5*self.E_j*(np.exp(-0.5*self.eta**2)*(np.exp(1j*self.eta*(np.conj(z[i])+z[j]))
                    +np.exp(-1j*self.eta*(np.conj(z[i])+z[j])))+self.eta**2*(1+np.conj(z[i])**2+z[j]**2
                    +2*np.conj(z[i])*z[j]))*np.conj(p[i])*p[j])*np.exp(-0.5*(np.abs(z[i])**2+np.abs(z[j])**2-2*np.conj(z[i])*z[j]))
        return f.real
    
    def norm(self,p,z):
        f = 0
        for i in range(self.ncs):
            for j in range(self.ncs):
                f+= np.conj(p[i])*p[j]*np.exp(-0.5*(np.abs(z[i])**2+np.abs(z[j])**2-2*np.conj(z[i])*z[j]))
        return f.real
    


    def z_dot(self,p,z):       # X is the output containing the first order time deriviative of real and imaginary parts of p and z
        M = np.zeros((self.ncs,self.ncs),dtype = 'complex64')
        for i in range(self.ncs):
            for j in range(self.ncs):
                # M[i,j] = np.dot(np.conj(coherent(self.ncs,z[i])),coherent(self.ncs,z[j]))
                M[i,j] = np.exp(-0.5*(np.abs(z[i])**2+np.abs(z[j])**2-2*np.conj(z[i])*z[j]))

        A = np.zeros((4*self.ncs,4*self.ncs),dtype = 'float64')
        A[:self.ncs,:self.ncs] = M.real         #1               #numbering based on moving from left to right in an eq. and then going to other eq.
        A[:self.ncs,self.ncs:2*self.ncs] = -M.imag     #2
        A[self.ncs:2*self.ncs,:self.ncs] = M.imag      #5
        A[self.ncs:2*self.ncs,self.ncs:2*self.ncs] = M.real   #6

        for i in range (self.ncs):
            for j in range (self.ncs):
                A[:self.ncs,2*self.ncs:3*self.ncs][i,j] = -0.5*((np.conj(z[j]-2*z[i])+z[j])*p[j]*M[i,j]).real    #3
                A[:self.ncs,3*self.ncs:4*self.ncs][i,j] = 0.5*((np.conj(z[j]-2*z[i])-z[j])*p[j]*M[i,j]).imag      #4
                A[self.ncs:2*self.ncs,2*self.ncs:3*self.ncs][i,j] = -0.5*((np.conj(z[j]-2*z[i])+z[j])*p[j]*M[i,j]).imag   #7
                A[self.ncs:2*self.ncs,3*self.ncs:4*self.ncs][i,j] = -0.5*((np.conj(z[j]-2*z[i])-z[j])*p[j]*M[i,j]).real    #8
                A[2*self.ncs:3*self.ncs,:self.ncs][i,j] = (z[j]*M[i,j]).real    #9
                A[2*self.ncs:3*self.ncs,self.ncs:2*self.ncs][i,j] = -(z[j]*M[i,j]).imag    #10
                A[2*self.ncs:3*self.ncs,2*self.ncs:3*self.ncs][i,j] = ((1-0.5*z[j]*np.conj(z[j]-2*z[i])-0.5*z[j]**2)*p[j]*M[i,j]).real    #11
                A[2*self.ncs:3*self.ncs,3*self.ncs:4*self.ncs][i,j] = -((1-0.5*z[j]*np.conj(z[j]-2*z[i])+0.5*z[j]**2)*p[j]*M[i,j]).imag   #12
                A[3*self.ncs:4*self.ncs,:self.ncs][i,j] = (z[j]*M[i,j]).imag    #13
                A[3*self.ncs:4*self.ncs,self.ncs:2*self.ncs][i,j] = (z[j]*M[i,j]).real   #14
                A[3*self.ncs:4*self.ncs,2*self.ncs:3*self.ncs][i,j] = ((1-0.5*z[j]*np.conj(z[j]-2*z[i])-0.5*z[j]**2)*p[j]*M[i,j]).imag    #15
                A[3*self.ncs:4*self.ncs,3*self.ncs:4*self.ncs][i,j] = ((1-0.5*z[j]*np.conj(z[j]-2*z[i])+0.5*z[j]**2)*p[j]*M[i,j]).real    #16

        B = np.zeros((4*self.ncs),dtype = 'float64')

        for i in range(self.ncs):
            f0 = 0
            f1 = 0
            for j in range(self.ncs):
                # f0 += p[j]*z[j]*np.conj(z[i])*M[i,j]            # for cavity
                f0 += (self.Omega_t*np.conj(z[i])*z[j]-0.5*self.E_j*(np.exp(-0.5*self.eta**2)*(np.exp(1j*self.eta*(np.conj(z[i])+z[j]))
                        +np.exp(-1j*self.eta*(np.conj(z[i])+z[j])))+self.eta**2*(1+np.conj(z[i])**2+z[j]**2+2*np.conj(z[i])*z[j])))*p[j]*M[i,j]
                # f1 += (1+np.conj(z[i])*z[j])*p[j]*z[j]*M[i,j]   # for cavity
                f1 += self.Omega_t*(1+np.conj(z[i])*z[j])*p[j]*z[j]*M[i,j]-0.5*self.E_j*( ( (1.5*z[i]-2*z[j])*(np.exp(-0.5*self.eta**2)*(np.exp(1j*self.eta*(np.conj(z[i])+z[j]))
                        +np.exp(-1j*self.eta*(np.conj(z[i])+z[j])))+self.eta**2*(1+np.conj(z[i])**2+z[j]**2+2*np.conj(z[i])*z[j]))
                        +self.eta*(1j*np.exp(-0.5*self.eta**2)*(np.exp(1j*self.eta*(np.conj(z[i])+z[j]))-np.exp(-1j*self.eta*(np.conj(z[i])+z[j])))
                        +2*self.eta*(np.conj(z[i])+z[j])) )*p[j]*M[i,j] + 1.5*(np.exp(-0.5*self.eta**2)*(np.exp(1j*self.eta*(np.conj(z[j])+z[i]))
                        +np.exp(-1j*self.eta*(np.conj(z[j])+z[i])))+self.eta**2*(1+np.conj(z[j])**2+z[i]**2+2*np.conj(z[j])*z[i]))*np.conj(p[j]/p[i])*p[i]*z[i]*M[j,i] )
            # B[:self.ncs][i] = (-1j*self.Omega_t*f0).real
            B[:self.ncs][i] = (-1j*f0).real
            B[self.ncs:2*self.ncs][i] = (-1j*f0).imag
            B[2*self.ncs:3*self.ncs][i] = (-1j*f1).real
            B[3*self.ncs:4*self.ncs][i] = (-1j*f1).imag

        X = linalg.solve(A,B)
        p_real_dot = X[:self.ncs]
        p_imag_dot = X[self.ncs:2*self.ncs]
        p_dot = p_real_dot + 1j*p_imag_dot
        z_real_dot = X[2*self.ncs:3*self.ncs]
        z_imag_dot = X[3*self.ncs:4*self.ncs]
        z_dot = z_real_dot + 1j*z_imag_dot
        return p_dot, z_dot
    

    def runge_kutta4(self,t_f,num_steps):
        dt = t_f/num_steps*np.ones((num_steps))
        p_val = []
        z_val = []
        p_val.append(self.p)
        z_val.append(self.z)
        for i in trange (num_steps):
            k1, l1 = self.z_dot(p_val[-1],z_val[-1])
            k2, l2 = self.z_dot(p_val[-1]+0.5*dt[i]*k1,z_val[-1]+0.5*dt[i]*l1)
            k3, l3 = self.z_dot(p_val[-1]+0.5*dt[i]*k2,z_val[-1]+0.5*dt[i]*l2)
            k4, l4 = self.z_dot(p_val[-1]+dt[i]*k3,z_val[-1]+dt[i]*l3)
            p_next = p_val[-1] + dt[i]*1/6*(k1+2*k2+2*k3+k4)
            z_next = z_val[-1] + dt[i]*1/6*(l1+2*l2+2*l3+l4)

            p_val.append( p_next )
            z_val.append( z_next )
        return p_val, z_val
    
    def result(self,z_vals):
        z_r = []
        z_i = []
        for i in range(self.ncs):
            for j in range(len(z_vals)):
                z_r.append(np.real(z_vals[j][i]))
                z_i.append(np.imag(z_vals[j][i]))
        return np.array(z_r), np.array(z_i)
    
    def length(self,p_vals,z_vals):
        length = []
        for i in range (len(z_vals)-1):
            length.append(self.norm(p_vals[i],z_vals[i]))
        return np.array(length)

    def E_t_tr(self,p_vals,z_vals):
        E_t_tr = []
        for i in range(len(z_vals)-1):
            E_t_tr.append(self.exp_Et(p_vals[i],z_vals[i]))
        return np.array(E_t_tr)
