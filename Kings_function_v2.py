
# coding: utf-8

# In[1]:

##What to do
#Include error bars - sqrt of counts
from math import e as e_constant
def Kings(pars, x):
    """ This should plot the Windowed Timming Mode Swift Point spread function
    as defined by Moretti et al. 2005. the Kings Function takes the form
    
    KING(r) = (1 +(r/rc)^(2))^(-beta)
    rc = core radius 
    beta = slope
    
    rc and beta are defined in the following way:
    coeff_ai + coeff_bi*theta + coeff_ci*Energy + coeffi_di*Energy*theta
    
    These coeffitions are defined in the header of the Swift psf function v5. 
    
    """
    rc =pars[0] + (pars[1]*pars[8]) + (pars[2]*pars[9])+(pars[3]*pars[9]*pars[8])
    beta = pars[4] + (pars[5]*pars[8]) + (pars[6]*pars[9])+(pars[7]*pars[9]*pars[8])
    return  pars[10]*((1+(x/rc)**2)**(-beta))


# In[7]:

def PSF_Ground_Simplify(pars, x):
    """ This should plot the Windowed Timming Mode Swift Point spread function
    as defined by 2011 CALDB psf file (swxpsf200010101v005.fits) in WT mode. This is function
    is also referenced in 2005 Moretti, psf ground based analysis. 
    
    PSF = P0*Gauss(P1) + (1-P0)*King(P2,P3)
    
    P_i = ai + (bi*theta) + (ci*energy)+ (di*energy*theta)
    KING(r) = (1 +(r/rc)^(2))^(-beta)
    rc and beta are defined in the following way:
    
    (P0,P1,P2,P3)===>(W,sigma,rc,beta)
    
    These coeffitions are defined in the header of the Swift psf function v5. 
    
    """    
    
    W = (pars[0] + (pars[1]*pars[16]) + (pars[2]*pars[17])+(pars[3]*pars[17]*pars[16]))
    sigma = (pars[4] + (pars[5]*pars[16]) + (pars[6]*pars[17])+(pars[7]*pars[17]*pars[16]))
    rc = (pars[8] + (pars[9]*pars[16]) + (pars[10]*pars[17])+(pars[11]*pars[17]*pars[16]))
    beta= (pars[12] + (pars[13]*pars[16]) + (pars[14]*pars[17])+(pars[15]*pars[17]*pars[16]))
    GaussianProfile = (e_constant**(-1*((x**2)/(2*(sigma)**2))))
    Kings = ((1+((x-pars[19])/rc)**2)**(-beta))
    return pars[18]*(W*GaussianProfile + (1-W)*Kings)


# In[6]:

def PSF_Ground_Simplify_version_3(pars, x):
    """ This should plot the Windowed Timming Mode Swift Point spread function
    as defined by 2004 CALDB psf file (swxpsf200010101v003.fits) in WT mode. 
    
    PSF = P0*Gauss(P1) + King(P2,P3)
    
    P_i = ai + (bi*theta) + (ci*energy)+ (di*energy*theta)
    KING(r) = (1 +(r/rc)^(2))^(-beta)
    rc and beta are defined in the following way:
    
    (P0,P1,P2,P3)===>(W,sigma,rc,beta)
    
    These coeffitions are defined in the header of the Swift psf function v3. 
    
    """    
    W = (pars[0] + (pars[1]*pars[16]) + (pars[2]*pars[17])+(pars[3]*pars[17]*pars[16]))
    sigma = (pars[4] + (pars[5]*pars[16]) + (pars[6]*pars[17])+(pars[7]*pars[17]*pars[16]))
    rc = (pars[8] + (pars[9]*pars[16]) + (pars[10]*pars[17])+(pars[11]*pars[17]*pars[16]))
    beta= pars[12] + (pars[13]*pars[16]) + (pars[14]*pars[17])+(pars[15]*pars[17]*pars[16])
    GaussianProfile = e_constant**(-1*((x**2)/(2*(sigma)**2)))
    Kings = (1+(x/rc)**2)**(-beta)
    return pars[18]*(W*GaussianProfile + Kings)

