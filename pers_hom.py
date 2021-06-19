import numpy as np
import os
import scipy.stats as st
import matplotlib.pyplot as plt
import scipy.spatial.distance as ssd 
import itertools as it

class ph:
    def __init__(self,my_data,my_hom_dim,my_thresh,my_dist_type=None,my_dist_mat = None,my_dist_max=None):
        self.data= my_data
        self.hom_dim = my_hom_dim
        self.thresh = my_thresh
        self.dist_type = my_dist_type
        self.dist_mat = my_dist_mat
        self.pp = None
        self.dist_max = my_dist_max
        self.birth=None
        self.death=None
        self.dims = None

    def build_distance(self,p=2):
        method = getattr(self,self.dist_type,lambda:"invalid distance type")
        self.dist_mat = method(p)
        self.dist_max = np.max(self.dist_mat)
    
    def pnorm(self,my_p):
        return ssd.squareform(ssd.pdist(self.data,metric='minkowski', p = my_p))

    def spherical_dist(self,p):
        diffMat = lambda x:x[:,np.newaxis]-x
        multMat = lambda x: x[:,np.newaxis]*x
        hs = lambda x: (1-np.cos(x))/2
        halpha = lambda x,y,z: x + y*z
        dis = lambda x: np.arcsin(np.sqrt(x))
        costheta = self.data[:,0]
        pol = np.arccos(costheta)	    
        az = self.data[:,1]

        hpol,haz = map(hs,list(map(diffMat,[pol,az])))
        cosmult = multMat(costheta)
        ha=halpha(hpol,cosmult, haz)
        my_dist_mat = dis(ha)
        return my_dist_mat

    def unit_circle_dist(self,p):
        return np.mod(np.abs(self.data[:,np.newaxis] - self.data),np.pi)

    def mat_pert_dist(self,data):
        dist_mat = np.zeros_like(data)
        for i in np.arange(dist_mat.shape[0]):
            for j in np.arange(i+1,dist_mat.shape[1]):
                dist_mat[j,i] = np.abs(data[i,j]/(data[i,i]-data[j,j]))
        return dist_mat + dist_mat.T


    def betti(self,radii_list):
        tmp = []
        curr_dim = 0
        for i in radii_list:
            tmp.append(curr_dim)
            if i == True:
                curr_dim +=1
                tmp.pop()
        return np.array(tmp)

    def pers_pairs(self,dat):
        radii_list = dat[:,1] == -1.0
        my_dims = self.betti(radii_list)
        birth_tmp = np.array(dat[:,0]) 
        my_birth = np.delete(birth_tmp, np.where(birth_tmp==-1.0))
        death_tmp = np.array(dat[:,1])
        death_tmp[death_tmp > 1000000] = self.dist_max
        my_death = np.delete(death_tmp, np.where(death_tmp==-1.0))
        self.birth, self.death, self.dims = my_birth, my_death,my_dims
    
    def run_ripser(self,input_str,output_str):
        np.savetxt(input_str, self.dist_mat, delimiter=",")
        if self.thresh > 0:
            os.system("ripser.exe {} --format distance --dim {} --threshold {} > {}".format(input_str,self.hom_dim,self.thresh,output_str))
            ripser_output=np.loadtxt(output_str,delimiter=",",skiprows= 1)
        else:
            os.system("ripser.exe {} --format distance --dim {} > {}".format(input_str,self.hom_dim,output_str))
            ripser_output=np.loadtxt(output_str,delimiter=",")
        self.pers_pairs(ripser_output)


    def plot_pp(self,title_str):
        tmp = np.unique(self.dims)
        for dim in tmp:
            b = self.birth[self.dims == dim]
            d = self.death[self.dims==dim]
            plt.scatter(b,d,label = "dim {}".format(dim))

        plt.xlabel("birth radius")
        plt.ylabel("death radius")
        x = np.arange(np.max(self.death))
        plt.plot(x,x,color='r',linestyle="dashed")
        plt.title(title_str)
        plt.legend()
        plt.ylim([0,5])
        plt.xlim([-.01,np.max(self.birth)])
        plt.show()

    def hist_pp(self,title_str):
        tmp = np.unique(self.dims)
        for dim in tmp:
            b = self.birth[self.dims == dim]
            d = self.death[self.dims==dim]
            plt.hist(d-b,label = "dim {}".format(dim),alpha=.3,density=True,bins=30)
        plt.xlabel("lifetime sum")
        plt.ylabel("density")
        plt.title(title_str)
        plt.legend()
        plt.show()

    def fractal_dim(self, alpha):
        un_dims = np.unique(self.dims)
        result = []
        for dim in un_dims:
            b= self.birth[self.dims == dim]
            d = self.death[self.dims == dim]
            result.append(np.sum((d-b)**alpha))
        return np.array(result)



    def pers_entropy(self,alpha):
        un_dims = np.unique(self.dims)
        result = []
        for dim in un_dims:
            b= self.birth[self.dims == dim]
            d = self.death[self.dims == dim]
            s = d-b
            prob = s/np.sum(s)
            nor = np.log(len(b))
            if alpha == 1:
                entropy = -np.sum(prob*np.log(prob))/nor 
            else:
                entropy = np.log(np.sum(prob**alpha))/(1-alpha)
            result.append(entropy)
        return np.array(result)


    '''
    given persistence diagram P, the rank function r(a,b) is the sum of persistence points to the north-west of (a,b); i.e., the number of homology groups born before a and die after b. We have to decide a consistent gridding, or do a linear interpolation.
    '''
    def rank_pp(self):
        return 0

    '''
    calculate euler integral from the distance matrix. If hom_dim = i, then we sum up through i+1 simplices.
    '''
    def euler_integral(self):
        result = 0
        result += self.dist_mat.shape[0]*self.dist_max 
        result -= np.sum(self.dist_mat)/2.0
        if self.hom_dim == 0:
            return result
        count = 2
        return 0
        # come back, implement a max over combinations lambda function
        # for np.arange(2,self.hom_dim+2)

    def inter_pcd(self, ph1):
        new_ph = ph(None,self.hom_dim,0,self.dist_type)
        if self.data.ndim == 1:
            new_ph.data = np.concatenate((self.data,ph1.data))
            tmp = len(self.data)
        else:
            new_ph.data = np.vstack((self.dist_mat,ph1.dist_mat))
            tmp = self.data.shape[0]

        new_ph.build_distance()
        new_ph.dist_max  = np.max(new_ph.dist_mat[:tmp,tmp:])
        new_ph.dist_mat[:tmp,:tmp] = new_ph.dist_max + 1 
        new_ph.dist_mat[tmp:, tmp:] = new_ph.dist_max + 1
        new_ph.thresh = new_ph.dist_max 

        return new_ph

    def get_simplicity(self, clusters):
        n = clusters.size
        den = np.sum(np.log(np.arange(1,n+1)))
        num = 0.0
        for cluster in clusters:
            num += np.sum(np.log(np.arange(1,cluster.size+1)))
        return num/den

    def simplicities(self):
        result = []
        #den = np.sum(np.log(np.arange(1,len(self.dims)+1)))
        # set birth and death radii for all connected components
        b = self.birth[self.dims == 0]
        d = self.death[self.dims == 0]
        pairs = np.array([b,d]).T
        pairs[pairs[:,1].argsort()]
        clusters = np.arange(1,len(pairs)+1).reshape([len(pairs),1])
        for i, pair in enumerate(pairs):
            # locate cluster containing pair[0]
            pair0InCluster = np.array(list(map(lambda y: np.isin(pair[0],y),clusters))).tolist()
            # locate cluster containing pair[1]
            pair1InCluster = np.array(list(map(lambda y: np.isin(pair[1],y),clusters))).tolist()
            # merge both clusters
            merged = np.concatenate(( clusters[pair0InCluster], clusters[pair1InCluster] ))
            # remove both unmerged clusters from list and append single merged cluster
            clusters = np.append( clusters[  (np.logical_not(pair0InCluster))
                                           & (np.logical_not(pair1InCluster))], merged )
            # compute "simplicity" for this configuration and append to results
            # (together with merge step and death radius)
            result.append([i, pair[1], self.get_simplicity(clusters)])
        return np.array(result)
            







