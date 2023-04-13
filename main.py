from km15gen import genOneEvent
import argparse
import numpy as np
import time

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument("-Ed", "--Ed", type = float, default = 10.604)
  parser.add_argument("-trig", "--trig", type = int, default = 1)
  parser.add_argument("-fname", "--fname", type = str, default = "km15gen")
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