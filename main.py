from km15gen import genOneEvent
import argparse
import numpy as np
import time

M  = 0.938272081
x1 = 1/2/M/8.604
x2 = 1/(5-M**2)
x3 = (10.604/8.604-1)/M*10.604* (1-np.cos(np.radians(35)))
x4 = (1-(4-M**2)/2/10.604/M)/(1+(4-M**2)/2/10.604**2/(1-np.cos(np.radians(35))))

y1 = 1
y2 = 1.456
y3 = 2.510
y4 = 4.326
y5 = 7.671

c0 = y2/2/M/8.604
d0 = 1/(1+(4-M*M)/y2)
c1 = np.sqrt(y2*y3)/2/M/8.604
d1 = 1/(1+(4-M*M)/np.sqrt(y2*y3))
c2 = y3/2/M/8.604
d2 =  1/(1+(4-M*M)/y3)
c3 = np.sqrt(y3*y4)/2/M/8.604
d3 = 1/(1+(4-M*M)/np.sqrt(y3*y4))
c4 = y4/2/M/8.604
d4 = 1/(1+(4-M*M)/y4)
c5 = np.sqrt(y4*y5)/2/M/8.604
d5 = 1/(1+(4-M*M)/np.sqrt(y4*y5))
c6 = y5/2/M/8.604
d6 = 1/(1+(4-M*M)/y5)

newxBbins = [x1, c0, c1, c2, c3, c4, d2]
newQ2bins = [y1, y2, np.sqrt(y2*y3), y3, np.sqrt(y3*y4), y4, np.sqrt(y4*y5)]
newxBbins2 = [x1, c0, c1, c2, c3, c4,                             c5, d2, d4]
newQ2bins2 = [y1, 1.2, y2, np.sqrt(y2*y3), y3, np.sqrt(y3*y4)    , y4, np.sqrt(y4*y5), 7]
newtbins = [0.11, 0.15, 0.25, 0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 1.79]

def main(args):

  M = 0.938272081 # target mass
  Ed    = args.Ed
  xBmin = args.xBmin
  xBmax = args.xBmax
  Q2min = args.Q2min
  Q2max = args.Q2max
  tmin  = args.tmin
  tmax  = args.tmax
  ymin  = args.ymin
  ymax  = args.ymax
  w2min  = args.w2min


  now = time.time()
  trig = args.trig
  num  = 0
  while num <trig:
    result   = genOneEvent(xBmin, xBmax, Q2min, Q2max, tmin, tmax, ymin, ymax, w2min, 0,  rad = 1, filename = args.fname)
    num      = num + result
  later = time.time()
  print(later-now)


  # if args.find_max:
  #   xs_array = []
  #   while len(xs_array) < args.trig:
  #     xs = bhdvcs(xBmin, xBmax, Q2min, Q2max, tmin, tmax, ymin, ymax, xs_only = 1)
  #     if xs:
  #       xs_array.append(xs)

  #   print(np.max(xs_array))
if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument("-Ed", "--Ed", type = float, default = 10.604)
  parser.add_argument("-trig", "--trig", type = int, default = 1)
  parser.add_argument("-fname", "--fname", type = str, default = "km15gen")
  parser.add_argument("-bin", "--bin", type = float, default = 0.05)
  parser.add_argument("-xBmin", "--xBmin", type = float, default = 0.05)
  parser.add_argument("-xBmax", "--xBmax", type = float, default = 0.75)
  parser.add_argument("-Q2min", "--Q2min", type = float, default = 0.9)
  parser.add_argument("-Q2max", "--Q2max", type = float, default = 11)
  parser.add_argument("-tmin", "--tmin", type = float, default = 0.085)
  parser.add_argument("-tmax", "--tmax", type = float, default = 1.79)
  parser.add_argument("-ymin", "--ymin", type = float, default = 0.19)
  parser.add_argument("-ymax", "--ymax", type = float, default = 0.85)
  parser.add_argument("-w2min", "--w2min", type = float, default = 3.61)
  args = parser.parse_args()

  main(args)