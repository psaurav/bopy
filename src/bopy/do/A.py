import sys

class A:
    def __init__(self, n_in=1.33, n_out=1.):
        #n_in = 1.33 #water
        #n_out = 1.485 #acrylic plastic
        #n_out = 1. #air
        #c_light = 299792458000 # mm/s
        #v = c_light/n_in
        self.n = float(n_in)/n_out
        print >>sys.stderr, "Message:", __name__, " -- n_in/n_out =", self.n

    def getContiniA(self):
        n = self.n
        n2 = n*n
        n3 = n2*n
        n4 = n2*n2 
        if n < 1:                      # Eq (A2)
            A =  3.084635 - 6.531194*n + 8.357854*n2 - 5.082751*n3 + 1.171382*n4
        elif n > 1:                    # Eq (A3)
            n5    = n2*n3
            n6    = n3*n3
            n7    = n3*n4
            A =  504.332889 - 2641.00214*n + 5923.699064*n2 - 7376.355814*n3 + 5507.53041*n4 - 2463.357945*n5 + 610.956547*n6 - 64.8047*n7
        else:
            A = 1.0
        return A

