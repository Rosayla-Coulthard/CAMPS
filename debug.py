import numpy as np
from scipy.spatial import KDTree

pts = np.array([(2, 1), (3, 1), (4, 1), (1, 2), (2, 2), (3, 2), (4, 2), (1, 3), (2, 3), (3, 3), (4, 3), (1, 4), (2, 4), (3, 4), (4, 4)])

T = KDTree(pts)
idx = T.query_ball_point([1,1],r=1.5)
print (pts[idx])

