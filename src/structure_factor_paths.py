import numpy as np

# Class for k-space paths through the first BZ used for structure factor
class Path(object):
  # init method for Path class
  # ARGS:
  #   xs <- list of k_x coordinates of path corners
  #   ys <- list of k_y coordinates of path corners
  def __init__(self,xs,ys):
    self.xs = xs
    self.ys = ys
    segments = [np.sqrt(np.square(xs[i+1]-xs[i])+np.square(ys[i+1]-ys[i]))\
                for i in range(len(xs)-1)]
    path_length = sum(segments)
    # for parametric evaluation, calculate what fraction of total path length
    # each segment is
    self.ts = [sum(segments[:i])/path_length for i in range(len(segments)+1)]
  
  # calculate (k_x,k_y) coordinate along path parameterized by 0 <= t <= 1
  # ARGS:
  #   t <- parametric evaluation parameter
  def coord(self,t):
    # find segment that t belongs to
    for i in range(1,len(self.ts)):
      if self.ts[i-1] <= t <= self.ts[i]:
        t_eff = (t-self.ts[i-1])/(self.ts[i]-self.ts[i-1])
        # interpolate
        x = self.xs[i]*t_eff + self.xs[i-1]*(1-t_eff)
        y = self.ys[i]*t_eff + self.ys[i-1]*(1-t_eff)
        return np.array([x,y])
    raise ValueError

# hardcoded paths. I know this is pretty gross...
pth1 = Path([0, 2*np.pi, 2*np.pi,            0,                  0],\
            [0,       0, 2*np.pi/np.sqrt(3), 4*np.pi/np.sqrt(3), 0])

pth2 = Path([0, -2*np.pi, -2*np.pi,            0,                  0],\
            [0,        0,  2*np.pi/np.sqrt(3), 4*np.pi/np.sqrt(3), 0])

pth3 = Path([            2*np.pi,           2*np.pi,                 0,            2*np.pi],\
            [-2*np.pi/np.sqrt(3),2*np.pi/np.sqrt(3),4*np.pi/np.sqrt(3),-2*np.pi/np.sqrt(3)])

pth4 = Path([           -2*np.pi,          -2*np.pi,                 0,           -2*np.pi],\
            [-2*np.pi/np.sqrt(3),2*np.pi/np.sqrt(3),4*np.pi/np.sqrt(3),-2*np.pi/np.sqrt(3)])

pth5 = Path([ -2*np.pi, 2*np.pi, 0, -2*np.pi],\
            [2*np.pi/np.sqrt(3),2*np.pi/np.sqrt(3),4*np.pi/np.sqrt(3),2*np.pi/np.sqrt(3)])

PATHS = [pth1,pth2,pth3,pth4,pth5]
