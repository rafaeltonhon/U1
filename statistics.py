import numpy as np
from scipy.optimize import curve_fit

def v(r,sigma,b,v0):
    return v0+sigma*r+b/r

# function that does the linear regression for a function y=ax+b
def linearregression(x,y):
    sumx=sumy=sumx2=sumxy=0
    n=len(x)
    for i in range(0,n):
        sumx+=x[i]
        sumx2+=x[i]**2
        sumy+=y[i]
        sumxy+=x[i]*y[i]

    a=-(sumx*sumy-n*sumxy)/(sumx**2-n*sumx2)
    b=(sumx*sumxy-sumx2*sumy)/(sumx**2-n*sumx2)
    return a, b


# function that does the linear regression for a function y=ax+b
# the error on the y variables are considered
def linearregression_yerr(x,y,sigmay):
    sumx=sumy=sumx2=sumxy=sumsigmay=0
    n=len(x)
    for i in range(0,n):
        sumx+=x[i]/sigmay[i]**2
        sumx2+=(x[i]/sigmay[i])**2
        sumy+=y[i]/sigmay[i]**2
        sumxy+=x[i]*y[i]/sigmay[i]**2
        sumsigmay+=1/sigmay[i]**2

    a=(sumsigmay*sumxy-sumy*sumx)/(sumsigmay*sumx**2-sumx2*sumx)
    return a
    

# function that computes the potential
def pot(a,b,filein,fileout):
    #file=open('mean-wilson.dat','r')
    file=open(filein,'r')
    w=list()
    v=list()
    # read the file
    for line in file.readlines():
        w.append([float(x) for x in line.split()])
    file.close()
    r=[]
    for i in range(0,b):
        r.append(i+1)

    for i in range(0,a):
        l=w[i][:]

        sigma, b=linearregression(r,l)
        v.append([sigma,jackknifefit(r,l)])
        # we make the fit
        #for j in range(1,b):
        #    if(w[i][j]/w[i][j-1]>0): l.append(-np.log(w[i][j]/w[i][j-1]))
        #v.append([np.mean(l),jackknife(l)])
        #v.append([np.mean(l),np.std(l)/np.sqrt(len(l))])


    #file=open('potqq.dat','w')
    file=open(fileout,'w')
    for i in range(0,a):
        file.write(" ".join([str(x) for x in [i+1,v[i][0],v[i][1],"\n"]]))
    file.close()


# function that computes the mean and the the std of W(C)
def meanwilson(a,b,p,fileout1,fileout2):
    # make the mean value of the wilson loops
    w=list()
    s=list()
    for i in range(1,a+1):
        # open the file
        if(p=='full'): file=open(f"fort.{100+i}",'r') # full loops
        if(p=='proj'): file=open(f"fort.{200+i}",'r') # vortex projected loops  
        # read the file inputs
        wi=[]
        for line in file.readlines():
            l=[]
            l=[float(x) for x in line.split()]
            wi.append(l[1:]) # append the usefull values of the array
            
        file.close()
        # now we make the mean of each column
        wm=[]
        sm=[]
        for column in range(0,b):
            l=[]
            for i in range(0,len(wi)):
                l.append(wi[i][column])
            wm.append(np.mean(l))
            sm.append(np.std(l)/np.sqrt(len(l)))
        w.append(wm)
        s.append(sm)
    # print w on the file
    #file=open('mean-wilson.dat','w')
    #file1=open('meanwilsonerror.dat','w')
    file=open(fileout1,'w')
    file1=open(fileout2,'w')        
    for i in range(0,a):
        file.write(" ".join(str(x) for x in w[i]))
        file.write('\n')
        for j in range(0,b):
            file1.write(" ".join(str(x) for x in [i+1,j+1,w[i][j],s[i][j],"\n"]))
        file1.write("\n")
    file.close()
    file1.close()


# function that computes the creutz ratios
def creutz(a,filein,fileout):
    file=open(filein,'r')
    w=list()
    chi=list()
    # read the file
    for line in file.readlines():
        w.append([float(x) for x in line.split()])
    file.close()
    # compute the ratios
    chi.append(-np.log(w[0][0]))
    for i in range(1,a):
        chi.append(-np.log((w[i][i]*w[i-1][i-1])/(w[i][i-1]*w[i-1][i])))
    # get out with data
    file=open(fileout,'w')
    for i in range(0,a):
        file.write(" ".join(str(x) for x in [i+1,chi[i],'\n']))


# function jackknefi
def jackknife(thetaN):
    thetan=list()
    thetahat=np.mean(thetaN)
    N=len(thetaN)
    for i in range(0,N):
        ln=[]
        ln=thetaN[:]
        del ln[i]
        thetan.append(np.mean(ln))
    sigmatheta=0
    for i in range(0,N):
        sigmatheta+=(thetan[i]-thetahat)**2
    sigmatheta*=((N-1)/N)
    return np.sqrt(sigmatheta)


# function jackknefi
def jackknifefit(xN,thetaN):
    thetan=list()
    a, b=linearregression(xN,thetaN)
    thetahat=a
    N=len(thetaN)
    for i in range(0,N):
        ln=[]
        LN=xN[:]
        ln=thetaN[:]
        del ln[i]
        del LN[i]
        a, b=linearregression(LN,ln)
        thetan.append(a)
    sigmatheta=0
    for i in range(0,N):
        sigmatheta+=(thetan[i]-thetahat)**2
    sigmatheta*=((N-1)/N)
    return np.sqrt(sigmatheta)

#def pot
#=============================================================================================================#
#=============================================================================================================#
#=============================================================================================================#
#               MAIN PROGRAM
#=============================================================================================================#
#=============================================================================================================#
#=============================================================================================================#
a=int(input("a value: "))
b=int(input("b value: "))
# for the full loops
meanwilson(a,b,'full','mean-wilson.dat','meanwilsonerror.dat')
pot(a,b,'mean-wilson.dat','potqq.dat')
#creutz(a,'mean-wilson.dat','creutz.dat')

# for the projected loops
proj=str(input("projetados [s/n] "))
if(proj[0]=='s'):
	meanwilson(a,b,'proj','mean-wilson-proj.dat','meanwilsonerror-proj.dat')
	pot(a,b,'mean-wilson-proj.dat','potqq-proj.dat')	
	creutz(a,'mean-wilson-proj.dat','creutz-proj.dat')
	
